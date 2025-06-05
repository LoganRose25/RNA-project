library(tdatools)
library(tidyverse)
library(plotly)
library(umap)
library(Rtsne)
library(Matrix)

setwd('/Users/matthewwheeler/Dropbox (UFL)/Research/Biomedical Research/RNA Structure Profililng')
source("./struc_prof_tools_v2.0.R")

#setwd("/Users/matthewwheeler/Dropbox (UFL)/TDA_tools/SetA_no_shape/")

dx<-0.01
hom <- 1
save_dir <- paste0("./output041524/B/",paste0("hom",hom))

#### Load the Data ######### 

#### Files TestSet A no Shape #########

filesTA <- list.files(path="./TDA_tools/TestSetB_no_shape", pattern="*no_shape.csv", full.names=TRUE, recursive=FALSE)
files<-filesTA[file.size(filesTA) > 0]
filesTA<-files
NTA<-length(files)
K <- 100
L <- 0
matNorms <- c()
for(i in 0:4){
  start <- i*K+1
  end <- min((i+1)*K,NTA)
  files <- filesTA[start:end]
  M <- max_length(files,thres=5,dx=dx)
  matTA <- get_matrix_landscapes_shifted(files,homology=hom,M,thres=5,dx=dx)
  label <- replicate(nrow(matTA),"Test_B_No_Shape")
  matNorm <- cbind(label,matTA[,1])
  matNorms <- rbind(matNorms,matNorm)
  print(paste0("Finished matrix group ",i))
}
matNormsTA <- as.data.frame(matNorms)
matNormsTA.1 <- cbind("RNA_Label"=as.numeric(gsub("_no_shape.csv","",row.names(matNormsTA))),matNormsTA)
matNormsTA.1 <- matNormsTA.1[order(matNormsTA.1$RNA_Label),]

#### Files TestSet A with Shape #########

filesTAS <- list.files(path="./TDA_tools/TestSetB_with_shape33/", pattern="*with_shape.csv", full.names=TRUE, recursive=FALSE)
files<-filesTAS[file.size(filesTAS) > 0]
filesTAS<-files
NTAS<-length(filesTAS)
K <- 100
L <- 0
matNorms <- c()
for(i in 0:4){
  start <- i*K+1
  end <- min((i+1)*K,NTAS)
  files <- filesTAS[start:end]
  M <- max_length(files,thres=5,dx=dx)
  matTA <- get_matrix_landscapes_shifted(files,homology=hom,M,thres=5,dx=dx)
  label <- replicate(nrow(matTA),"Test_B_With_Shape")
  matNorm <- cbind(label,matTA[,1])
  matNorms <- rbind(matNorms,matNorm)
  print(paste0("Finished matrix group ",i))
}
matNormsTAS <- as.data.frame(matNorms)
matNormsTAS.1 <- cbind("RNA_Label"=as.numeric(gsub("_with_shape.csv","",row.names(matNormsTAS))),matNormsTAS)
matNormsTAS.1 <- matNormsTAS.1[order(matNormsTAS.1$RNA_Label),]
#### Files TrainSet A no Shape #########

filesTrA <- list.files(path="./TDA_tools/TrainSetB_no_shape/", pattern="*no_shape.csv", full.names=TRUE, recursive=FALSE)
files<-filesTrA[file.size(filesTrA) > 0]
filesTrA<-files
NTrA<-length(filesTrA)
K <- 100
L <- 0
matNorms <- c()
for(i in 0:10){
  start <- i*K+1
  end <- min((i+1)*K,NTrA)
  files <- filesTrA[start:end]
  M <- max_length(files,thres=5,dx=dx)
  matTA <- get_matrix_landscapes_shifted(files,homology=hom,M,thres=5,dx=dx)
  label <- replicate(nrow(matTA),"Train_B_No_Shape")
  matNorm <- cbind(label,matTA[,1])
  matNorms <- rbind(matNorms,matNorm)
  print(paste0("Finished matrix group ",i))
}
matNormsTrA <- as.data.frame(matNorms)
matNormsTrA.1 <- cbind("RNA_Label"=as.numeric(gsub("_no_shape.csv","",row.names(matNormsTrA))),matNormsTrA)
matNormsTrA.1 <- matNormsTrA.1[order(matNormsTrA.1$RNA_Label),]

#### Files TrainSet A with Shape #########

filesTrAS <- list.files(path="./TDA_tools/TrainSetB_with_shape33", pattern="*with_shape.csv", full.names=TRUE, recursive=FALSE)
files<-filesTrAS[file.size(filesTrAS) > 0]
filesTrAS <- files
NTrAS<-length(filesTrAS)
K <- 100
L <- 0
matNorms <- c()
for(i in 0:10){
  start <- i*K+1
  end <- min((i+1)*K,NTrAS)
  files <- filesTrAS[start:end]
  M <- max_length(files,thres=5,dx=dx)
  matTA <- get_matrix_landscapes_shifted(files,homology=hom,M,thres=5,dx=dx)
  label <- replicate(nrow(matTA),"Train_B_With_Shape33")
  matNorm <- cbind(label,matTA[,1])
  matNorms <- rbind(matNorms,matNorm)
  print(paste0("Finished matrix group ",i))
}

matNormsTrAS <- as.data.frame(matNorms)
matNormsTrAS.1 <- cbind("RNA_Label"=as.numeric(gsub("_with_shape.csv","",row.names(matNormsTrAS))),matNormsTrAS)
matNormsTrAS.1 <- matNormsTrAS.1[order(matNormsTrAS.1$RNA_Label),]


write.csv(matNormsTA.1,file.path(save_dir,"Norms_TestB.csv"))
write.csv(matNormsTAS.1,file.path(save_dir,"Norms_TestBS.csv"))
write.csv(matNormsTrA.1,file.path(save_dir,"Norms_TrainB.csv"))
write.csv(matNormsTrAS.1,file.path(save_dir,"Norms_TrainBS.csv"))

