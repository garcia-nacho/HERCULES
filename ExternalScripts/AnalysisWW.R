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
#n.start<-1250
#n.end<-2250
n.start<-as.numeric(path[3])
n.end<- as.numeric(path[4])
path<-path[1]
noise.path<-list.files(path, pattern = ".*noise.tsv",full.names = TRUE)
sampleid<-gsub("\\.noise.tsv","",gsub(".*/","",noise.path))
bam.path<-list.files(path, pattern = ".*bam$",full.names = TRUE)
ref.path<-paste(path,"spike.cons.aligned.fa",sep = "")
cores<-as.numeric(detectCores())-2
poi.co<-  (n.end-n.start)*0.03
  

allpos<-TRUE

noise.table<-read.csv(noise.path, sep = "\t", header = FALSE)

ggplot(noise.table)+
  geom_line(aes(V1,V2))+
  xlab("Position")+
  ylab("Noise")+
  geom_hline(yintercept=co.n, linetype='dotted', col = 'red')+
  theme_minimal()+
  ggtitle(paste("Noise",sampleid))

ggsave(paste(path, "Noise.pdf",sep = ""))

ggplot(noise.table)+
  geom_line(aes(V1,V3))+
  xlab("Position")+
  ylab("Depth")+
  theme_minimal()+
  ggtitle(paste("Coverage",sampleid))

ggsave(paste(path, "Coverage.pdf",sep = ""))

poi<-noise.table$V1[which(noise.table$V2>co.n & noise.table$V1>n.start & noise.table$V1<n.end & noise.table$V3>20)]
poit<-noise.table[which(noise.table$V2>co.n & noise.table$V1>n.start & noise.table$V1<n.end & noise.table$V3>20),]
if(length(poi)>0){
  
  if(length(poi)>poi.co){
    poi<-poit$V1[order(poit$V2, decreasing=TRUE)][1:round(poi.co)]
  }
  
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
  
  top.noise<-grep(paste("_P",noise.table$V1[which(noise.table$V3>20 & noise.table$V2== max(noise.table$V2[which(noise.table$V1 %in% poi)]))],".tsv",sep = ""), position.files)
  
  position.files<-position.files[c(top.noise, c(1:length(position.files))[-top.noise])]
  
  print("Getting the reads for the positions of interest")
  pb2<-txtProgressBar(min = 0, max = length(position.files), initial = 1)
  
  if(exists("out")) try(rm(out))
  for (i in 1:length(position.files)) {
    setTxtProgressBar(pb2,i)
    dummy<-fread(position.files[i],sep = "\t", header = FALSE)
    
    if(nrow(dummy)>0){
      #dummy<-as.data.frame(dummy)
      colnames(dummy)<-c("ReadID","Position","Base")
      if(length( which(duplicated(dummy)))>0)dummy<-dummy[-which(duplicated(dummy)),]
      #if(length(which(dummy$Base=="I"))>0)dummy<-dummy[-which(dummy$Base=="I"),]
      
      colnames(dummy)[3]<-paste("P_",unique(dummy$Position),sep = "")
      dummy<-dummy[,c(1,3)]
      
      if(length(which(is.na(dummy[,2])))<nrow(dummy)*0.7){
        if(!exists("out")){
          out<-dummy
          
        }else{
          setkey(out, ReadID)
          setkey(dummy, ReadID)
          out<-out[dummy, on="ReadID"]
          #out<-merge(out, dummy, by="ReadID",all=TRUE)
        }
      }
    }
  }
  
  
  
  
  #Removing problematic reads
  out<-as.data.frame(out)
  col.to.remove<-which(apply(out[,-1],2,function(x)length(which(is.na(x))))> nrow(out)*0.7 )
  if(length(col.to.remove)>0) out<-out[,-which(colnames(out) %in% names(col.to.remove))]
  write.csv(out, "MutationMatrix.csv", row.names = FALSE)
  
  file.remove(position.files)
  
  
  index.c<-vector()
  for (c in 2:ncol(out)) {
    if(length(which(out[,c] %in% c("D","I")))>0){ 
      if(max(table(out[-which(out[,c] %in% c("D","I")),c]))> nrow(out)*0.85) index.c<-c(index.c,c)}
  }
  
  if(length(index.c)>0) out<-out[,-index.c]
  if(length(which(is.na(out)))>0) out<-out[-unique(which(is.na(out), arr.ind = TRUE)[,1]),]
  
  print("Finding variants")
  
  out$variant<-NA
  out$varianthash<-NA
  variants<-vector()
  
  out$variant<-apply(out[,grep("P_", colnames(out))], 1, paste, collapse="/")
  out$varianthash<-apply(out[,grep("P_", colnames(out))], 1, function(x) digest(paste(x, collapse = "/"), algo = "md5"))
  
  variantable<-as.data.frame(table(out$variant))
  
  read.co<-sum(variantable$Freq)/50
  
  variantable<-variantable[which(variantable$Freq>read.co),]
  
  if(length(grep("//", variantable$Var1))>0) variantable<-variantable[-grep("//", variantable$Var1),]
  if(allpos & length(grep("NA", variantable$Var1))>0)variantable<-variantable[-grep("NA", variantable$Var1),]
  
  if(exists("variant.out")) try(rm(variant.out))
  for (i in 1:nrow(variantable)) {
    dum<-out[which(out$variant==variantable$Var1[i])[1],-1]
    dum$variantID<- paste(sampleid,"-", paste(unlist(base::strsplit(out$varianthash[which(out$variant==variantable$Var1[i])[1]],""))[1:12],collapse = "") ,sep = "")
    dum$count<- variantable$Freq[i]
    dum$ratio<-variantable$Freq[i]/sum(variantable$Freq)
    
    if(!exists("variant.out")){
      variant.out<-dum
    }else{
      variant.out<-rbind(variant.out, dum)
    }
  }
  if(length(grep("D", variant.out$variant))>0){ 
    pure.variants<-variant.out[-grep("D", variant.out$variant),]
  }else{
    pure.variants<-variant.out
  }
  row.t.r<-vector()
  for (i in 1:nrow(variant.out)) {
    if(!variant.out$variantID[i] %in% pure.variants$variantID){
      signature<- as.character(variant.out[i,grep("P_", colnames(variant.out))])
      dummy.pure.var<-pure.variants[,grep("P_", colnames(pure.variants))[-which(signature=="D")], drop=FALSE]  
      dummy.pure.var$ID<- apply(dummy.pure.var, 1, function(x) paste(x,collapse = "/"))
      var.to.track<- which(dummy.pure.var$ID == paste(signature[-which(signature=="D")],collapse =   "/"))
      #Changed this to ==1 
      if(length(var.to.track)==1){ 
        variant.out$variantID[i]<-paste(pure.variants$variantID[var.to.track],collapse = "/")
        pure.variants$count[var.to.track]<- variant.out$count[i]+pure.variants$count[var.to.track]
        pure.variants$ratio[var.to.track]<- variant.out$ratio[i]+pure.variants$ratio[var.to.track]
      }
    }
  }
  
  if(length(which(variant.out$variantID %in% pure.variants$variantID))>0){ 
    variant.out<-variant.out[-which(variant.out$variantID %in% pure.variants$variantID),]
  }
  
  if(nrow(variant.out)>0){
    if(nrow(pure.variants)>0){ 
      pure.variants<-rbind(pure.variants, variant.out)
    }else{
      pure.variants<-variant.out
    }  
  }
  variant.out<-pure.variants
  
  variant.out$variant<-NULL
  variant.out$varianthash<-NULL
  variant.out$ratio<- variant.out$count/sum(variant.out$count)
  write_xlsx(variant.out, paste(path,"VariantResults.xlsx", sep = ""))
  
  
  refspike<-read.fasta(ref.path)
  reference<-as.character(refspike[[which(names(refspike)=="Spike")]])
  consensus<-as.character(refspike[[-which(names(refspike)=="Spike")]])
  
  consensus<-consensus[21563:25384]
  reference<-reference[21563:25384]
  consensus[which(consensus=="n")]<-reference[which(consensus=="n")]
  
  pb3<-txtProgressBar(min = 0, max = nrow(variant.out),initial = 1)
  
  refspike<-refspike[-which(names(refspike)=="Spike")]
  names(refspike)<-gsub("_consensus.*","",names(refspike))
  refspike[[1]]<-consensus
  for (i in 1:nrow(variant.out)) {
    setTxtProgressBar(pb3,i)
    positions<-as.numeric(gsub("P_","",colnames(variant.out)[grep("P_", colnames(variant.out))]))
    dummy<-consensus
    dummy[positions]<-as.character(variant.out[i,grep("P_", colnames(variant.out))])
    if(length(which(as.character(variant.out[i,grep("P_", colnames(variant.out))]) %in% c("D","I")))>0){
      dummy[which(dummy %in% c("D","I"))]<-consensus[which(dummy %in% c("D","I"))]
    }
    
    if(length(which(dummy=="-"))>0)dummy<-dummy[-which(dummy=="-")]
    refspike[[1+i]]<-dummy[c(n.start:n.end)]
    names(refspike)[length(refspike)]<-paste("ID-",variant.out$variantID[i],sep = "")
  }
  if(length(which(refspike[[1]]=="-"))>0)refspike[[1]]<-refspike[[1]][-which(refspike[[1]]=="-")]
  write.fasta(refspike, paste(path,"Variants.fa",sep = ""), names = names(refspike))
}
