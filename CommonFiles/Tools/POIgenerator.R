library(data.table)
library(ggplot2)
library("doParallel")
library("parallel")
library("foreach")
library(doSNOW)
library(progress)
library(digest)
library(writexl)
library(seqinr)
path=commandArgs(TRUE)
co.n<-0.15
pathf<-path[1]
co.n<-as.numeric(path[2])
n.start<-as.numeric(path[3])
n.end<- as.numeric(path[4])
poicurrent<- path[5]
if(poicurrent!="0" & poicurrent!="auto") poicurrent<-as.numeric(gsub(" ","",unlist(base::strsplit(poicurrent, ","))))


noise.path<-list.files(pathf, pattern = ".*noise.tsv",full.names = TRUE, recursive = TRUE)

if(exists("noise.table")) rm(noise.table)
for (i in 1:length(noise.path)) {
  if(exists("noise.table.d"))rm(noise.table.d)
  try(noise.table.d<-read.csv(noise.path[i], sep = "\t", header = FALSE), silent = TRUE)
  if(exists("noise.table.d")){
  noise.table.d$Sample<-gsub(".*/","",gsub(".noise.tsv","",noise.path[i]))
  if(!exists("noise.table")){
    noise.table<-noise.table.d
  }else{
    noise.table<-rbind(noise.table,noise.table.d)
  }
  }
}



noise.table.agg<-aggregate(V2~V1, noise.table, mean)
depth.table.agg<-aggregate(V3~V1, noise.table, mean)

poi.co<-  length(depth.table.agg$V3>20)*0.025

poit<-noise.table.agg[which(noise.table.agg$V2>co.n & noise.table.agg$V1>n.start & noise.table.agg$V1 < n.end & depth.table.agg$V3>20),]
poi<-poit$V1[order(poit$V2, decreasing=TRUE)][1:min(round(poi.co),nrow(poit))]
if(poicurrent[1]!="0")poi<-unique(c(poi,poicurrent))
if(length(which(poi>n.end))>0) poi<-poi[-which(poi>n.end)]
if(length(which(poi<n.start))>0) poi<-poi[-which(poi<n.start)]

if(poicurrent=="auto"){
  for (j in 1:length(noise.path)) {
    if(exists("df2")) rm(df2)
    try(df2<-read.csv(noise.path[j], sep = "\t", header = FALSE) , silent = TRUE)
   if(exists("df2")){
    
     den<-density(df2$V2)
     bins<-den$y
     border2<-bins[c(2*c(1:250))]
     border1<-bins[c(2*c(1:250))-1]
     border1<-as.data.frame(border1)
     border1$border2<-border2
     border1$mean<-apply(border1, 1, mean)
     border1$sd<-apply(border1, 1, sd)
     border1$cov<-border1$sd/border1$mean
     borderco<-den$x[min(which(border1$cov>0.5 & as.numeric(rownames(border1))>15))]
     poi<-c(poi, df2$V1[which(df2$V2 > borderco & df2$V3>50)])
    }
  }
}

poi<-unique(poi)

write.csv(poi, "poi.temp",row.names = FALSE)
