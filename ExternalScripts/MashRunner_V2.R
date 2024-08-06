library(phylotools)
library("doParallel")
library("parallel")
library("foreach")
library(doSNOW)
library(progress)

fasta.to.read<-list.files(pattern = "uncompressed.fasta")
reads<-phylotools::read.fasta(fasta.to.read)

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
  coord<-templist$V3
  newlin<-gsub(".fasta","",gsub(".*_","",templist$V1[which(as.numeric(templist$V3)==min(as.numeric(templist$V3)))]))
  
  if(length(newlin)==1){
    return(list(newlin,coord)) 
  }else{
    if(length(unique(newlin))==1){
      return(list(unique(newlin),coord))
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
    return(list(newlin,coord))
  }
  
}

stopCluster(cluster.cores)


out.vec<-lapply(out.par, function(x)x[[1]])
out.vec<-unlist(out.vec)
if(length(which(out.vec==".X"))>0) out.vec[which(out.vec==".X")]<-"Unclassified"
if(length(which(is.na(out.vec)))>0) out.vec<-out.vec[-which(is.na(out.vec))]

out.umap<-lapply(out.par, function(x)as.numeric(x[[2]]))

empty<-unlist(lapply(out.umap, function(x) length(x) ))
empty<- which(empty==0)
if(length(empty)>0){
  out.vec<-out.vec[-empty]
  out.umap<-out.umap[-empty]
  
}

out.umap <- do.call("rbind", out.umap)

non.covid<-apply(out.umap, 1, function(x)length(which(x==1)))
NoSC2<-vector()
if(length(which(non.covid==ncol(out.umap)))>0) NoSC2<-c(NoSC2,which(non.covid==ncol(out.umap)))

mashdistvec<-apply(out.umap, 1, function(x)min(x))

if(length(which(mashdistvec>=0.2)) >0) NoSC2<-c(NoSC2,which(mashdistvec>=0.2))

if(length(NoSC2)>0){
  NoSC2<-unique(NoSC2)
  out.vec[NoSC2]<-"Non-SARS-CoV-2"
}

lineages<-as.data.frame(table(out.vec))
lineages$Count<-lineages$Freq
lineages$Freq<-lineages$Freq/sum(lineages$Freq)

out.umap<-as.data.frame(out.umap)
out.umap$Lineage<-out.vec

write.csv(lineages, gsub(".uncompressed.fasta","_AF.lineages.csv",fasta.to.read), row.names = FALSE)
 reads<-phylotools::read.fasta(fasta.to.read)
 if(length(empty)>0) reads<-reads[-empty,]
 out.umap$SeqID<-reads$seq.name

#write.csv(out.umap, gsub(".uncompressed.fasta","_MashMatrix.csv",fasta.to.read),row.names = FALSE)

