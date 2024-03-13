library(data.table)
library(seqinr)
library("doParallel")
library("parallel")
library("foreach")
library(doSNOW)
library(progress)

positions<-c(   "T22032C", "C22033A", "A22034G")

positions<-toupper(positions)

positions<-gsub("A|T|C|G", "", positions)
positions<-as.numeric(positions)
if(length(which(is.na(positions)))>0){
  positions<-positions[-which(is.na(positions))]
}

positions<-positions[order(positions)]

if(positions[1]>21563 & positions[1]<25384 ){
  positions<- positions-21563+1+12
}else if(positions[1]>1 & positions[1] <38321){
  positions<-positions
}else{
  positions<-"auto"
}



files<-list.files("/media/nacho/Data/wastewater/Kmersearchsource/", pattern = "*.tsv$", full.names = TRUE)
ref<-read.fasta("/media/nacho/Data/DockerImages/Wastewater_SARS-CoV-2/CommonFiles/reference/SpikeRef.fa")
ref<-unlist(ref)

#23012-21563+1

kmer<-c("TAAGCATAGTG")
kmer.names<-c("BA.2.86")


counter<-0
pb<-txtProgressBar(min = 1, initial = 1, max = length(kmer)*length(files))
for (k in 1:length(kmer)) {


  df<-as.data.frame(files)
  kmer<-toupper(kmer)
  df$kmer<-kmer[k]
  df$countQuery<-NA
  df$countRef<-NA
  df$RatioQuery<-NA
  df$kmerName<-kmer.names[k]

  for (i in 1:nrow(df)) {
     gc()
     counter<-counter+1
     setTxtProgressBar(pb,counter)
     file.to.check<-fread(df$files[i])
     file.to.check<-as.data.frame(file.to.check)
     
     file.to.check$`#query_msa`<-toupper(file.to.check$`#query_msa`)
     file.to.check$ref_msa<-toupper(file.to.check$ref_msa)
     index.query<-grep(kmer[k],file.to.check$`#query_msa`)
     df$countQuery[i]<-length(index.query)
     df$countRef[i]<-length(grep(kmer[k],file.to.check$ref_msa))
     df$RatioQuery[i]<-df$countQuery[i]/nrow(file.to.check)
     df$RatioReference[i]<-df$countRef[i]/nrow(file.to.check)
     if(length(index.query)>0 & kmer.names[k]=="BA.2.86") write.table(file.to.check[index.query,], gsub("tsv","BA.2.86.tsv",df$files[i]),sep = "\t", quote = FALSE, row.names = FALSE)
     
  }
  
  if(!exists("df.out")){
    df.out<-df
  }else{
    df.out<-rbind(df.out, df)
  }

}

library(ggplot2)


write.csv(df.out, "/media/nacho/Data/wastewater/BA.2.86.results.csv")

df.out<-read.csv("/media/nacho/Data/wastewater/BA.2.86.results.csv")
df.out$Sample<-gsub(".*/","",df.out$files)
df.out$Sample<-gsub(".sorted.*","",df.out$Sample)
df.out<-df.out[-which(df.out$kmerName=="WT"),]
df.out$Location<-gsub(".*_","",df.out$Sample)
df.out$Week<-gsub("_.*","",df.out$Sample)
df.out<-df.out[-grep("Pos|Neg",df.out$Sample),]
df.out$Run<-gsub("\\..*","",df.out$Week)
df.out$Week<-as.numeric(gsub(".*\\.","",df.out$Week))

ggplot(df.out)+
  geom_bar(aes(Location,RatioQuery ),stat = "identity", fill="red")+
  theme_minimal()+
  facet_wrap(~Run)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ggtitle("Ratio G22895C + T22896A + G22898A")+
  ylab("Ratio")
