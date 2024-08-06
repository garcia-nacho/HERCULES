library(seqinr)

sam<-list.files(pattern ="*voc.sam" )
seq<-read.fasta("/Data/VariantSpike.fasta")

all<-read.csv(sam, sep = "\t", header = FALSE, skip = length(seq)+1,  quote = "")
if(length(which(all$V3==""))>0) all<-all[-which(all$V3==""),]

#if(length(grep("^SID",all$V3))>0) all$V3<-gsub(".*_","",all$V3)
results<-as.data.frame(table(all$V3))
 
write.csv(results, gsub("voc.sam","voc_depth.csv",sam),row.names = FALSE  )
