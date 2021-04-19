
setwd("C:\\4_R\\ToxValue\\Prediction\\Prediction_Script")
args <- commandArgs(TRUE)	# this line will take the command line parameter into args.
process_id <- args[1]			# this line will send args[1] to the infile

#process_id <- "pi_438159"
model_input <- paste("../Prediction_temp_files/", process_id,"_input.txt",sep="")

# input file should only have 1 SMILES in it
# if there are more than one, only the first will be used
newchem.smiles<- scan(model_input,what=character())[1]

message("input read")
# load libraries
# install.packages("randomForest")  # If the library is not installed
# install.packages("rcdk")  		# If the library is not installed
# install.packages("rJava")  		# If the library is not installed
# .libPaths("C:/Users/wchiu/Documents/R/win-library/3.3") # For testing
library("randomForest")
library("rcdk") # Note this requres JAVA that is that same build (e.g., 64-bit) as R

## Set up -- only do this once
message("libraries loaded")
load("ToxValRFModels.Rdata")
tv.residuals <- read.csv("Prediction_Residuals.csv",row.names=1)
row.names(tv.residuals) <- ToxValuesNames
dc <- get.desc.categories()
dnames<-c(get.desc.names(dc[1]),get.desc.names(dc[2]),get.desc.names(dc[3]),get.desc.names(dc[4]),get.desc.names(dc[5]))
dnames<-unique(dnames)
load_time<-Sys.time()
message("config loaded")
## Get SMILES and check at least 1 carbon and no metals/metaloids

mol <- parse.smiles(newchem.smiles)[[1]]
atomicnumbers <- unique(unlist(lapply(get.atoms(mol),get.atomic.number)))

nocarbon <- !(6 %in% atomicnumbers)
metals_etc <- c(5,13:14,21:33,39:52,57:71,89:103,72:84,104:112)
hasmetals <- sum(metals_etc %in% atomicnumbers) > 0
message(" smiles obtained")
## Get descriptors from smiles - do this for each chemical

mol.desc<-cbind(data.frame(smiles=newchem.smiles,stringsAsFactors = FALSE),eval.desc(mol,dnames))

message(" descriptor obtained")
## Clean up descriptors
mol.desc.clean <- mol.desc[,col.keep]
frac.imputed <- sum(is.na(mol.desc.clean))/length(mol.desc.clean)

# Impute NA with median of training set (if necessary)
mol.desc.clean.imputed <- mol.desc.clean
for (i in 1:ncol(mol.desc.clean.imputed)) {
  mol.desc.clean.imputed[is.na(mol.desc.clean.imputed[,i]),i]<-col.med[i]
}

## Warnings

if (nocarbon | hasmetals | (frac.imputed > 0.1)) {
  warnfile <- paste("..\\Prediction_temp_files\\", process_id,"_warn.txt",sep="")
  fileConn<-file(warnfile,open="wt")
  if (nocarbon) writeLines("Warning: Inorganic (no carbons) - model not intended for this use",fileConn)
  if (hasmetals) writeLines("Warning: Contains metals/metalloids - model not intended for this use",fileConn)
  if (frac.imputed > 0.1) writeLines("Warning: >10% of molecular descriptors imputed",fileConn)
  close(fileConn)
}

message("warnings done")
## Make predictions
x.new<-mol.desc.clean.imputed[,-1]
y.pred<-data.frame(tv=sub(".train","",ToxValuesNames),row.names=ToxValuesNames)
y.pred$prediction<-0
y.pred$lower95<-0
y.pred$upper95<-0
y.pred$appl.domain<-0





for (tvname in ToxValuesNames) {

  y.pred[tvname,"prediction"] <- predict(tv.model[[tvname]],x.new)
  y.pred[tvname,"lower95"] <- y.pred[tvname,"prediction"] + tv.residuals[tvname,"CDK.ci.lb"]
  y.pred[tvname,"upper95"] <- y.pred[tvname,"prediction"] + tv.residuals[tvname,"CDK.ci.ub"]
}

## Check applicability domain Z-score

for (tvname in ToxValuesNames) {
 message(tvname)
 
 #cat(tvname, file="output.text", append=TRUE)
  ## Scale based on tox value descriptors
  cas.train <- tv.train[[tvname]][,1]
  desc.train <- tv.train[[tvname]][,-(1:2)]
  lo<-apply(desc.train,2,min,na.rm=TRUE)
  hi<-apply(desc.train,2,max,na.rm=TRUE)
  # only keep descriptors that hi-lo > 0
  notconstant <- (hi!=lo)
  desc.train.scale <- as.data.frame(scale(desc.train[,notconstant],
                                          center=lo[notconstant],
                                          scale=(hi-lo)[notconstant]))

  ## Write ".x" file for tox values
  tv_xfile<-paste("../Prediction_temp_files/", process_id, "_intermediate_",tvname,".x",sep="")
  write.table(as.data.frame(t(dim(desc.train.scale))),file=tv_xfile,
              row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t")
  write.table(as.data.frame(t(names(desc.train.scale))),file=tv_xfile,
              row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t",append=TRUE)
  write.table(cbind(cas.train,desc.train.scale),file=tv_xfile,
              row.names=TRUE,col.names=FALSE,quote=FALSE,sep="\t",append=TRUE)
  
  ## Write ".x" file for new chemical
  newchem_xfile<- paste("..\\Prediction_temp_files\\", process_id, "_intermediate", "_newchem.x",sep="")	# file name here
  x.new.scale <- as.data.frame(scale(x.new[,notconstant],
                                     center=lo[notconstant],
                                     scale=(hi-lo)[notconstant]))
  write.table(as.data.frame(t(dim(x.new.scale))),file=newchem_xfile,
              row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t")
  write.table(as.data.frame(t(names(x.new.scale))),file=newchem_xfile,
              row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t",append=TRUE)
  write.table(cbind(newchem.smiles,x.new.scale),file=newchem_xfile,
              row.names=TRUE,col.names=FALSE,quote=FALSE,sep="\t",append=TRUE)
  #message("predicition 1 part done")
  ## Calculate z-score
  gad_outfile <- paste("..\\Prediction_temp_files\\", process_id, "_intermediate", "_newchem_by_",tvname,".gad",sep="");
  system(paste("get_ad.exe ",tv_xfile," -4PRED=",newchem_xfile," -Z=3 -OUT=",gad_outfile,sep="")); 
  report.dat<-read.table(gad_outfile,header=TRUE,skip=3);
  y.pred[tvname,"appl.domain"]<-report.dat$Z.score
  #message ("z-score done")
}


##
outfile <- paste("..\\Prediction_temp_files\\", process_id,"_output.csv",sep="")
write.csv(y.pred,file=outfile,row.names=FALSE)

message("prediction done")
# remove intermediate files
setwd("C:/4_R/ToxValue/Prediction/Prediction_temp_files/")
intermediate_file_pattern <- paste(process_id, "_intermediate_.*..*..*",sep="")
junk <- dir(path="./", pattern= intermediate_file_pattern) 
file.remove(junk)

message("junk deleted")
#close(fileConn)