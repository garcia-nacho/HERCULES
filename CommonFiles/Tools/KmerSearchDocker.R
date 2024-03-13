library(data.table)
library(seqinr)
library("doParallel")
library("parallel")
library("foreach")
library(doSNOW)
library(progress)
library(ggplot2)
library(writexl)

arg=commandArgs(TRUE)

compressed<-list.files(full.names = TRUE, pattern = "b2f.tsv.gz",recursive = TRUE)
#ref<-read.fasta("/media/nacho/Data/DockerImages/Wastewater_SARS-CoV-2/CommonFiles/reference/SpikeRef.fa")
#kmers<-c("TAAGCATAGTG")


kmers<-arg[1]
ref<-read.fasta("/home/docker/CommonFiles/reference/SpikeRef.fa")
ref<-unlist(ref)

dir.create("kmeruncompressed")

for (i in 1:length(compressed)) {
  target.file<-gsub(".*/","kmeruncompressed/",gsub(".gz","",compressed[i]))
  if(file.exists(target.file)) file.remove(target.file)
  system(paste("gunzip -c ", compressed[i], " > ",target.file,sep = ""))
}

kmer<-base::strsplit(kmers[1],",")
files<-list.files("kmeruncompressed",full.names = TRUE )
counter<-0
pb<-txtProgressBar(min = 1, initial = 1, max = length(kmer)*length(files))

for (k in 1:length(kmer)) {

  df<-as.data.frame(files)
  kmer<-toupper(kmer)
  df$kmer<-kmer[k]
  df$countQuery<-NA
  df$countRef<-NA
  df$RatioQuery<-NA

  for (i in 1:nrow(df)) {
     gc()
    if(!dir.exists("ExtractedKmer")) dir.create("ExtractedKmer")
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
     if(length(index.query)>0 ) write.table(file.to.check[index.query,], 
                                            paste(gsub(".tsv",paste("_",kmer[k],".tsv",sep=""),gsub(".*/","ExtractedKmer/",df$files[i])),sep = ""),sep = "\t", quote = FALSE, row.names = FALSE)
     
  }
  df$files<-gsub(".*/","",df$files)
  
  if(!exists("df.out")){
    df.out<-df
  }else{
    df.out<-rbind(df.out, df)
  }

}


df.out$Sample<-gsub(".*/","",df.out$files)
df.out$Sample<-gsub(".sorted.*","",df.out$Sample)

write_xlsx(df.out, "KmerSearchResults.xlsx")

df.out$Location<-gsub(".*_","",df.out$Sample)
df.out$Week<-gsub("_.*","",df.out$Sample)
df.out<-df.out[-grep("Pos|Neg",df.out$Sample),]
df.out$Run<-gsub("\\..*","",df.out$Week)
df.out$Week<-as.numeric(gsub(".*\\.","",df.out$Week))


ggplot(df.out)+
  geom_bar(aes(Location,RatioQuery ),stat = "identity", fill="red")+
  theme_minimal()+
  facet_wrap(~Run+kmer)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ylab("Ratio")

ggsave("KmerSearchResult.pdf")
