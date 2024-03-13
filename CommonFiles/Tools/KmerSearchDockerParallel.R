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

kmer<-unlist(base::strsplit(kmers,","))
files<-list.files("kmeruncompressed",full.names = TRUE )

for (k in 1:length(kmer)) {
  if(!dir.exists("ExtractedKmer")) dir.create("ExtractedKmer")
  df<-as.data.frame(files)
  kmer<-toupper(kmer)
  df$kmer<-kmer[k]
  df$countQuery<-NA
  df$countRef<-NA
  df$RatioQuery<-NA
  df$RatioReference<-NA

  samples.to.analyze<-gsub(".*/","",gsub("\\.sorted.*","",compressed))
  
  pb <- progress_bar$new(
    format = "Sample: :samp.pb [:bar] :elapsed | eta: :eta",
    total = length(samples.to.analyze),    # 100 
    width = 80)
  samp <- samples.to.analyze
  
  progress <- function(n){
    pb$tick(tokens = list(samp.pb = samp[n]))
  } 
  
  opts <- list(progress = progress)
  
  gc()
  
  cores<-as.numeric(detectCores())-2
  cluster.cores<-makeCluster(cores)
  registerDoSNOW(cluster.cores)
  
  out.par<-foreach(i=1:length(compressed), .verbose=FALSE, .packages = c("data.table","seqinr"), .options.snow = opts) %dopar%{
     df2<-df
     gc()

     file.to.check<-fread(df$files[i])
     file.to.check<-as.data.frame(file.to.check)
     
     file.to.check$`#query_msa`<-toupper(file.to.check$`#query_msa`)
     file.to.check$ref_msa<-toupper(file.to.check$ref_msa)
     index.query<-grep(kmer[k],file.to.check$`#query_msa`)
     
     df2$countQuery[i]<-length(index.query)
     df2$countRef[i]<-length(grep(kmer[k],file.to.check$ref_msa))
     df2$RatioQuery[i]<-df2$countQuery[i]/nrow(file.to.check)
     df2$RatioReference[i]<-df2$countRef[i]/nrow(file.to.check)
     if(length(index.query)>0 ) write.table(file.to.check[index.query,], 
                                            paste(gsub(".tsv",paste("_",kmer[k],".tsv",sep=""),gsub(".*/","ExtractedKmer/",df2$files[i])),sep = ""),sep = "\t", quote = FALSE, row.names = FALSE)
     rm(file.to.check)
     return(df2[i,])
     
  }
  stopCluster(cluster.cores)
  dfunlist<-do.call(rbind, out.par)
  
  if(!exists("df.out")){
    df.out<-dfunlist
  }else{
    df.out<-rbind(df.out, dfunlist)
  }

}

file.remove(list.files("kmeruncompressed/",full.names = TRUE) )

df.out$Sample<-gsub(".*/","",df.out$files)
df.out$Sample<-gsub(".sorted.*","",df.out$Sample)

write_xlsx(df.out, "KmerSearchResults.xlsx")

df.out$Location<-gsub(".*_","",df.out$Sample)
df.out$Week<-gsub("_.*","",df.out$Sample)
df.out<-df.out[-grep("Pos|Neg",df.out$Sample),]
df.out$Run<-gsub("\\..*","",df.out$Week)
df.out$Week<-as.numeric(gsub(".*\\.","",df.out$Week))


ggplot(df.out)+
  geom_bar(aes(Run,RatioQuery ),stat = "identity", fill="red")+
  theme_minimal()+
  facet_wrap(~kmer+Location)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ylab("Ratio")
ggsave("KmerSearchBarplot.pdf" ,width = 20, height = 16)


