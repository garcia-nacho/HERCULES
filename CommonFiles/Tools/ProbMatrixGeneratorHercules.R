library(seqinr)
library(digest)

arg=commandArgs(TRUE)
#1: Cutoff 0.1
#2: start 1250
#3: end 2250
#4: poicurrent "auto"

nstart<-as.numeric(arg[1])
nend<- as.numeric(arg[2])
input.path<-"/Reference"

# wu<-seqinr::read.fasta("/media/nacho/Data/DockerImages/Wastewater_SARS-CoV-2/CommonFiles/reference/SpikeRef.fa")
# nstart<-1250
# nend<-2250
# input.path<-"/media/nacho/Data/DockerImages/HERCULES/CommonFiles/reference/newref"


wu<-seqinr::read.fasta("/home/docker/CommonFiles/reference/SpikeRef.fa")

input.file <- list.files(input.path, pattern = "*\\.fa.*",full.names = TRUE)
wu<-wu[[1]][c(nstart:nend)]

#sqs<-read.fasta("/media/nacho/Data/wastewater/paper/data/new_ref/norwegian_ref/HERCULES_GISAID_Spike_Norway.fasta")


if(length(input.file)==1){
  sqs<-seqinr::read.fasta(input.file)  
  sqs<-lapply(sqs, function(x) x[c(nstart:nend)] )
  
hashes<-lapply(sqs, function(x) digest(paste(x,collapse = ""),algo="md5"))
uniqueregions<-unique(unlist(hashes))
df<-as.data.frame(uniqueregions)
colnames(df)<-"hash"
df$QC<-NA
df$Corrected<-NA
for (i in 1:length(uniqueregions)) {
  pangos<- gsub(".*_","",names(sqs)[ which(hashes== uniqueregions[i])])
  pangos.df<-as.data.frame(table(pangos))
  pangos.df$Freq<-pangos.df$Freq/sum(pangos.df$Freq)
  pangos.df$Hash<-paste(strsplit(uniqueregions[[i]],"")[[1]][1:10], collapse = "")
  
  if(length(unique(pangos))==1){
    df$QC[which(df$hash==uniqueregions[i])]<-"OK"
    df$Corrected[which(df$hash==uniqueregions[i])]<-pangos[1]

  }else{

    df$QC[which(df$hash==uniqueregions[i])]<-"Warning"
    
    lineages<- pangos
    lineages<-strsplit(lineages,"\\.")
    noe<-max(unlist(lapply(lineages, length)))
    maxnoe<-noe
    
    running<-TRUE
    while(running){
      new.lineages<- unlist(lapply(lineages, function(x) paste(x[(1:min(length(x),noe))], collapse = ".") ))  
      ratio<- max(table(new.lineages))/length(new.lineages)
      
      if(ratio>0.9){
        new.lin<-names(table(new.lineages)[which(table(new.lineages)==max(table(new.lineages)))])[1]
        #unique.index<-c(unique.index,which(df$hash==uniqueregions[i])[1])
        running<-FALSE
        
      }else{
        new.lin<-NA
        if(noe==1) running<-FALSE
      }
      noe<-noe-1
    }
    if(is.na(new.lin)){
      df$Corrected[which(df$hash==uniqueregions[i])]<-paste(new.lin,".",paste(strsplit(uniqueregions[[i]],"")[[1]][1:10], collapse = "") ,".X",sep = "")
    }else{
      if(maxnoe==noe+1){
        df$Corrected[which(df$hash==uniqueregions[i])]<-new.lin
      }else{
        df$Corrected[which(df$hash==uniqueregions[i])]<-paste(new.lin,".X",sep = "")
      }
    }
    
  }
  
  pangos.df$Lineage<-df$Corrected[which(df$hash==uniqueregions[i])][1]
  if(!(exists("pangos.out"))){
    pangos.out<-pangos.df
  }else{
    pangos.out<-rbind(pangos.out, pangos.df)
  }
  
}


df$seqName<-NA
hashes<-unlist(hashes)
for (i in 1:length(sqs)) {
  names(sqs)[i]<- paste(df$Corrected[which(df$hash==hashes[i])], "_SQ",i,sep = "")
}

if(length(which(duplicated(hashes)))>0){
  sqs<-sqs[-which(duplicated(hashes))]
}


write.fasta(sqs, "/home/docker/CommonFiles/reference/Custom_ReferencesAligned.fasta",names = names(sqs))
write.csv(pangos.out, "/home/docker/CommonFiles/reference/variant_hash.csv", row.names = FALSE)

spikes<-lapply(sqs,function(x) paste(toupper(x),collapse = ""))

spikes_tsv<-as.data.frame(unlist(spikes))

colnames(spikes_tsv)<-"#query_msa"
#spikes_tsv$`#query_msa`<-paste("A",spikes_tsv$`#query_msa`,sep="")
spikes_tsv$ref_msa<-tolower(paste(wu,collapse = ""))

write.table(spikes_tsv, "/home/docker/CommonFiles/reference/MSA_Refs.tsv", row.names = FALSE, quote = FALSE, sep = "\t")
write.csv(rownames(spikes_tsv),"/home/docker/CommonFiles/reference/MSA_RefsID.csv",row.names = FALSE)
try(system("rm /home/docker/CommonFiles/reference/MSA_Refs.tsv.gz"))
system("gzip -f /home/docker/CommonFiles/reference/MSA_Refs.tsv")

poi<-c(1:length(wu))

noise.ref<-do.call(rbind, sqs)
noisecalc<-function(x){
  if(length(which(x=="n"))>0)x<-x[-which(x=="n")]
  return(1-(as.numeric(table(x)[(which(table(x)==max(table(x))))]/sum(table(x))))) 
}

noiseunique<-apply(noise.ref,2,noisecalc)
noise.co<-which(noiseunique>0.05)

if(length(noise.co>0)){
  write.csv(noise.co,"/home/docker/CommonFiles/reference/NoiseRefs.csv",row.names = FALSE)
}else{
  write.csv(1841,"/home/docker/CommonFiles/reference/NoiseRefs.csv",row.names = FALSE)
}



lineages<-unique(gsub("_.*","", names(sqs)))

bases<-c("a","t","c","g","-")
poid<-poi+nstart-1
pb<-txtProgressBar(min = 1, max = length(poi),initial = 1)
for (i in 1:length(poi)) {
  setTxtProgressBar(pb,i)
  dummym<-matrix(0, ncol = length(lineages), nrow = 5 )
  
  base.m<-lapply(sqs, function(x)x[poi[i]])  
  names(base.m)<-gsub("_.*","", names(base.m))
  base.m<-unlist(base.m)
  base.m2<-match(tolower(base.m), c("a","t","c","g","-"))
  names(base.m2)<-names(base.m)
  
  for (j in 1:length(lineages)) {
    for (k in c(1:5)) {
      dummym[k,j]<-length(which(names(base.m2)==lineages[j] & base.m2==k))/ length(which(names(base.m2)==lineages[j]))
    }
  }
  
  dummym<-as.data.frame(dummym)
  colnames(dummym)<-lineages
  rownames(dummym)<-paste(poid[i],c("A","T","C","G","D"),sep = "")
  if(!exists("dummy.out")){ 
    dummy.out<-dummym
  }else{
    dummy.out<-rbind(dummy.out, dummym)
  }
}

print("Writing Probability matrix")
write.csv(dummy.out, "/home/docker/CommonFiles/reference/ProbMatrix.csv" )
}