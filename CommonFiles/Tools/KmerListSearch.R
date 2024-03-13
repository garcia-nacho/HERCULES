library(data.table)
library(seqinr)
library("doParallel")
library("parallel")
library("foreach")
library(doSNOW)
library(progress)
library(ggplot2)
library(writexl)


compressed<-list.files(full.names = TRUE, pattern = "b2f.tsv.gz",recursive = TRUE)

#kmerlist<-read.csv("/media/nacho/Data/DockerImages/HERCULES/CommonFiles/reference/kmerlist.csv")
#ref<-read.fasta("/media/nacho/Data/DockerImages/HERCULES/CommonFiles/reference/SpikeRef.fa")

kmerlist<-read.csv("/home/docker/CommonFiles/reference/kmerlist.csv")
ref<-read.fasta("/home/docker/CommonFiles/reference/SpikeRef.fa")

ref<-unlist(ref)

#Uncompressing and cleaning

dir.create("kmeruncompressed")
pb<-txtProgressBar(max = length(compressed))
for (i in 1:length(compressed)) {
  setTxtProgressBar(pb,i)
  target.file<-gsub(".*/","kmeruncompressed/",gsub(".gz","",compressed[i]))
  if(file.exists(target.file)) file.remove(target.file)
  system(paste("gunzip -c ", compressed[i], " > ",target.file,sep = ""))
  
  file.to.check<-fread(target.file)
  file.to.check<-as.data.frame(file.to.check)
  for (j in 1:nrow(file.to.check)) {
    if(length(grep("-",file.to.check$ref_msa[j]))==1){
      ref<-unlist(strsplit(file.to.check$ref_msa[j],""))
      qry<-unlist(strsplit(file.to.check$`#query_msa`[j],""))
      qry<-qry[-which(ref=="-")]
      ref<-ref[-which(ref=="-")]
      file.to.check$ref_msa[j]<-paste(ref, collapse = "")
      file.to.check$`#query_msa`[j]<-paste(qry, collapse="")
    }
  }
  write.table(file.to.check,target.file,sep = "\t", quote = FALSE, row.names = FALSE)
  
}
#Recompress db and save

files<-list.files("kmeruncompressed",full.names = TRUE )


samples.to.analyze<-gsub(".*/","",gsub("\\.sorted.*","",files))


counter<-paste("kmerID:",c(1:nrow(kmerlist)),":Lineage:",kmerlist$Lineage,sep = "")

pb <- progress_bar$new(
  format = ":samp.pb [:bar] :elapsed | eta: :eta",
  total = length(counter),
  width = 80)

samp <- counter

progress <- function(n){
  pb$tick(tokens = list(samp.pb = samp[n]))
} 

opts <- list(progress = progress)

gc()

cores<-as.numeric(detectCores())-2
cluster.cores<-makeCluster(cores)
registerDoSNOW(cluster.cores)

out.par<-foreach(k=1:nrow(kmerlist), .verbose=FALSE, .packages = c("data.table"), .options.snow = opts) %dopar%{
  
  df<-as.data.frame(files)
  kmer<-toupper(kmerlist$Seq[k])
  df$kmer<-toupper(kmerlist$Seq[k])
  df$countQuery<-NA
  df$countRef<-NA
  df$RatioQuery<-NA
  df$RatioReference<-NA
  df$Lineage<-kmerlist$Lineage[k]
  
  for ( i in 1:length(files)) {
    file.to.check<-fread(df$files[i])
    #file.to.check<-as.data.frame(file.to.check)
    
    file.to.check$`#query_msa`<-toupper(file.to.check$`#query_msa`)
    file.to.check$ref_msa<-toupper(file.to.check$ref_msa)
    
    index.query<-grep(kmer,file.to.check$`#query_msa`)
    
    df$countQuery[i]<-length(index.query)
    df$countRef[i]<-length(grep(kmer,file.to.check$ref_msa))
    df$RatioQuery[i]<-df$countQuery[i]/nrow(file.to.check)
    df$RatioReference[i]<-df$countRef[i]/nrow(file.to.check)
  }
  
  return(df)  
}

stopCluster(cluster.cores)
dfunlist<-do.call(rbind, out.par)

df.out<-dfunlist[-which(dfunlist$RatioReference>0.98),]

df.out$Location<-gsub(".*_","",gsub(".*/","",gsub("_sorted.*","",df.out$files)))
df.out$Date<-as.Date(gsub(".*\\.","",gsub("_.*","",gsub(".*/","",gsub("_sorted.*","",df.out$files)))))
df.out$Run<-gsub("\\..*","",gsub(".*/","",df.out$files))


ggplot(df.out[which(df.out$Lineage=="BA.2.86"),])+
  geom_line(aes(Date, RatioQuery, group=kmer), alpha=0.2)+
  theme_minimal()+
  facet_wrap(~Location)

#file.remove(list.files("kmeruncompressed/",full.names = TRUE) )

df.out$Sample<-gsub(".*/","",df.out$files)
df.out$Sample<-gsub(".sorted.*","",df.out$Sample)

write_xlsx(df.out, "KmerSearchResults.xlsx")




