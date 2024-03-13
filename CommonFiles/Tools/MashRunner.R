library(phylotools)
library("doParallel")
library("parallel")
library("foreach")
library(doSNOW)
library(progress)

fasta.to.read<-list.files(pattern = "uncompressed.fasta")
reads<-read.fasta(fasta.to.read)
reads<-reads$seq.text
reads<-tolower(reads)
samples.to.analyze<-c(1:length(reads))

pb <- progress_bar$new(
  format = "Read: :samp.pb [:bar] :elapsed | eta: :eta",
  total = length(samples.to.analyze),    # 100 
  width = 60)

samp <- samples.to.analyze

progress <- function(n){
  pb$tick(tokens = list(samp.pb = samp[n]))
} 

opts <- list(progress = progress)

###
gc()

cluster.cores<-makeCluster(9)
registerDoSNOW(cluster.cores)
if(exists("out.par")) rm(out.par)
out.par<-foreach(i=1:length(samples.to.analyze), .verbose=FALSE, .packages = "seqinr", .options.snow = opts) %dopar%{
  
  templist<-strsplit(system(paste('printf ">Read\n',reads[i],'\n" | ./mash dist ./reference.msh -',sep = "" ), intern = TRUE),"\t") 
  
  templist<-as.data.frame(do.call("rbind", templist))
  newlin<-gsub(".fasta","",gsub(".*_","",templist$V1[which(as.numeric(templist$V3)==min(as.numeric(templist$V3)))]))
  
  if(length(newlin)==1){
    return(newlin) 
  }else{
    if(length(unique(newlin))==1){
      return(unique(newlin))
    }
    
    newlin.s<-strsplit(newlin, "\\.")
    shortest<-which(unlist(lapply(newlin.s, length))== min(unlist(lapply(newlin.s, length))))[1]
    
    for (nl in length(newlin.s[[shortest]]):1){
      temp.newlin<-paste(newlin.s[[shortest]][1:nl],collapse = ".")
      
      if(length(grep(temp.newlin,newlin))==length(newlin.s)){
        newlin<-temp.newlin
        newlin<-paste(newlin,".X",sep = "")
        break
      }
      
    }
    if(length(newlin)>1) newlin<-"Unclassified"
    return(newlin)
  }
  
}

stopCluster(cluster.cores)

out.vec<-unlist(out.par)
out.vec[which(out.vec==".X")]<-"Unclassified"

lineages<-as.data.frame(table(out.vec))
lineages$Count<-lineages$Freq
lineages$Freq<-lineages$Freq/sum(lineages$Freq)

write.csv(lineages, gsub(".uncompressed.fasta","_AF.lineages.csv",fasta.to.read), row.names = FALSE)


