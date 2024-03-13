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
co.n<-as.numeric(path[2])
path<-path[1]
noise.path<-list.files(path, pattern = ".*noise.tsv",full.names = TRUE)
sampleid<-gsub("\\.noise.tsv","",gsub(".*/","",noise.path))
bam.path<-list.files(path, pattern = ".*bam$",full.names = TRUE)
ref.path<-paste(path,"spike.cons.aligned.fa",sep = "")

allpos<-TRUE

noise.table<-read.csv(noise.path, sep = "\t", header = FALSE)

ggplot(noise.table)+
  geom_line(aes(V1,V2))+
  xlab("Position")+
  ylab("Noise")+
  geom_hline(yintercept=co.n, linetype='dotted', col = 'red')+
  theme_minimal()
ggsave(paste(path, "Noise.pdf",sep = ""))

ggplot(noise.table)+
  geom_line(aes(V1,V3))+
  xlab("Position")+
  ylab("Depth")+
  theme_minimal()
ggsave(paste(path, "Coverage.pdf",sep = ""))
