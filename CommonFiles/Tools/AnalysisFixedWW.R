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
path<-path[1]
noise.path<-list.files(path, pattern = ".*noise.tsv",full.names = TRUE)
bam.path<-list.files(path, pattern = ".*bam$",full.names = TRUE)
ref.path<-paste(path,"spike.cons.aligned.fa",sep = "")
cores<-8
read.co<-30
allpos<-TRUE

samples.to.analyze<-poi

pb <- progress_bar$new(
  format = "Position: :samp.pb [:bar] :elapsed | eta: :eta",
  total = length(samples.to.analyze),    # 100 
  width = 60)

samp <- samples.to.analyze

progress <- function(n){
  pb$tick(tokens = list(samp.pb = samp[n]))
} 

opts <- list(progress = progress)

###
gc()
cluster.cores<-makeCluster(cores)
registerDoSNOW(cluster.cores)

out.par<-foreach(k=1:length(poi), .verbose=FALSE, .options.snow = opts) %dopar%{
  system(paste(path,"bbasereader -p ",poi[k]  ,
               " -f ", bam.path, "> ", gsub("\\.sorted.bam","_P",bam.path), poi[k],".tsv",sep = ""))
}

stopCluster(cluster.cores)

position.files<-list.files(path, full.names = TRUE, pattern = "_P.*\\.tsv")

pb<-txtProgressBar(min = 0, max = length(position.files), initial = 1)
print("Getting the reads for the positions of interest")
try(rm(out))
for (i in 1:length(position.files)) {
  setTxtProgressBar(pb,i)
  dummy<-fread(position.files[i],sep = "\t", header = FALSE)
  
  if(nrow(dummy)>0){
    #dummy<-as.data.frame(dummy)
    colnames(dummy)<-c("ReadID","Position","Base")
    if(length( which(duplicated(dummy)))>0)dummy<-dummy[-which(duplicated(dummy)),]
    #if(length(which(dummy$Base=="I"))>0)dummy<-dummy[-which(dummy$Base=="I"),]
    which(duplicated(dummy$ReadID))
    colnames(dummy)[3]<-paste("P_",unique(dummy$Position),sep = "")
    dummy<-dummy[,c(1,3)]
    
    if(!exists("out")){
      out<-dummy
      
    }else{
      setkey(out, ReadID)
      setkey(dummy, ReadID)
      setkey(out, ReadID)
      out<-out[dummy, on="ReadID"]
      #out<-merge(out, dummy, by="ReadID",all=TRUE)
    }
  }
}
write.csv(out, "MutationMatrix.csv", row.names = FALSE)
close(pb)

file.remove(position.files)
out<-as.data.frame(out)

index.c<-vector()
for (c in 2:ncol(out)) {
  if(max(table(out[-which(out[,c] %in% c("D","I")),c]))> nrow(out)*0.85) index.c<-c(index.c,c)
}

if(length(index.c)>0) out<-out[,-index.c]
out<-out[-unique(which(is.na(out), arr.ind = TRUE)[,1]),]

print("Finding variants")

out$variant<-NA
out$varianthash<-NA
variants<-vector()

out$variant<-apply(out[,grep("P_", colnames(out))], 1, paste, collapse="/")
out$varianthash<-apply(out[,grep("P_", colnames(out))], 1, function(x) digest(paste(x, collapse = "/"), algo = "md5"))

variantable<-as.data.frame(table(out$variant))
variantable<-variantable[which(variantable$Freq>read.co),]

if(length(grep("//", variantable$Var1))) variantable<-variantable[-grep("//", variantable$Var1),]
if(allpos & length(grep("NA", variantable$Var1))>0)variantable<-variantable[-grep("NA", variantable$Var1),]

try(rm(variant.out))
for (i in 1:nrow(variantable)) {
  dum<-out[which(out$variant==variantable$Var1[i])[1],-1]
  dum$variantID<- paste(unlist(base::strsplit(out$varianthash[which(out$variant==variantable$Var1[i])[1]],""))[1:10] ,collapse = "")
  dum$count<- variantable$Freq[i]
  dum$ratio<-variantable$Freq[i]/sum(variantable$Freq)
  
  if(!exists("variant.out")){
    variant.out<-dum
  }else{
    variant.out<-rbind(variant.out, dum)
  }
}


pure.variants<-variant.out[-grep("D", variant.out$variant),]

for (i in 1:nrow(variant.out)) {
  if(!variant.out$variantID[i] %in% pure.variants$variantID){
    signature<- as.character(variant.out[i,grep("P_", colnames(variant.out))])
    dummy.pure.var<-pure.variants[,grep("P_", colnames(pure.variants))[-which(signature=="D")]]  
    dummy.pure.var$ID<- apply(dummy.pure.var, 1, function(x) paste(x,collapse = "/"))
    var.to.track<- which(dummy.pure.var$ID == paste(signature[-which(signature=="D")],collapse =   "/"))
    if(length(var.to.track)>0) variant.out$variantID<-pure.variants$variantID[var.to.track]
    pure.variants$count[var.to.track]<- variant.out$count[i]+pure.variants$count[var.to.track]
    pure.variants$ratio[var.to.track]<- variant.out$ratio[i]+pure.variants$ratio[var.to.track]
  }
}
variant.out<-variant.out[-which(variant.out$variantID %in% pure.variants$variantID),]

if(nrow(variant.out)>0) pure.variants<-rbind(pure.variants, variant.out)  

variant.out<-pure.variants

variant.out$variant<-NULL
variant.out$varianthash<-NULL
variant.out$ratio<- variant.out$count/sum(variant.out$count)
write_xlsx(variant.out, paste(path,"VariantResultsFixedPos.xlsx", sep = ""))


refspike<-read.fasta(ref.path)
reference<-as.character(refspike[[which(names(refspike)=="Spike")]])
consensus<-as.character(refspike[[-which(names(refspike)=="Spike")]])

consensus<-consensus[21563:25384]
reference<-reference[21563:25384]
consensus[which(consensus=="n")]<-reference[which(consensus=="n")]

pb<-txtProgressBar(min = 0, max = nrow(variant.out),initial = 1)

refspike<-refspike[-which(names(refspike)=="Spike")]
names(refspike)<-gsub("_consensus.*","",names(refspike))
refspike[[1]]<-consensus
for (i in 1:nrow(variant.out)) {
  setTxtProgressBar(pb,i)
  positions<-as.numeric(gsub("P_","",colnames(variant.out)[grep("P_", colnames(variant.out))]))
  dummy<-consensus
  dummy[positions]<-as.character(variant.out[i,grep("P_", colnames(variant.out))])  
  if(length(grep("D",toupper(dummy)))>0){
    dummy[grep("D",toupper(dummy))]<-as.character(refspike[[1]])[grep("D",toupper(dummy))]
  }
  if(length(which(dummy=="-"))>0)dummy<-dummy[-which(dummy=="-")]
  refspike[[1+i]]<-dummy
  names(refspike)[length(refspike)]<-paste("ID-",variant.out$variantID[i],sep = "")
}
if(length(which(refspike[[1]]=="-"))>0)refspike[[1]]<-refspike[[1]][-which(refspike[[1]]=="-")]
write.fasta(refspike, paste(path,"VariantsFixedPos.fa",sep = ""), names = names(refspike))