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
library(reshape2)

path=commandArgs(TRUE)
co.n<-0.15
co.n<-as.numeric(path[2])
n.start<-1
n.end<- 3822
poicurrent<-path[3]
path<-path[1]

SpikeRef<-"/home/docker/CommonFiles/reference/SpikeRef.fa"
#SpikeRef<-"/media/nacho/Data/DockerImages/Wastewater_SARS-CoV-2/CommonFiles/reference/SpikeRef.fa"

if(poicurrent!="0" & poicurrent!="auto") poicurrent<-as.numeric(gsub(" ","",unlist(base::strsplit(poicurrent, ","))))

noise.path<-list.files(path, pattern = ".*noise.tsv",full.names = TRUE, recursive = TRUE)
if(exists("noise.table")) rm(noise.table)
for (i in 1:length(noise.path)) {
  if(exists("noise.table.d"))rm(noise.table.d)
  try(noise.table.d<-read.csv(noise.path[i], sep = "\t", header = FALSE), silent = TRUE)
  if(exists("noise.table.d")){
  noise.table.d$Sample<-gsub(".*/","",gsub(".noise.tsv","",noise.path[i]))
  if(!exists("noise.table")){
    noise.table<-noise.table.d
  }else{
    noise.table<-rbind(noise.table,noise.table.d)
  }
  }
}

#co.n<-co.n/length(noise.path)

noise.table.agg<-aggregate(V2~V1, noise.table, mean)
depth.table.agg<-aggregate(V3~V1, noise.table, mean)



poi.co<-  length(depth.table.agg$V3>20)*0.025
poit<-noise.table.agg[which(noise.table.agg$V2>co.n & noise.table.agg$V1>n.start & noise.table.agg$V1 < n.end & depth.table.agg$V3>20),]
poi<-poit$V1[order(poit$V2, decreasing=TRUE)][1:min(round(poi.co),nrow(poit))]

if(poicurrent!="0"& poicurrent!="auto") poi<-c(poi,poicurrent)

# if(poicurrent=="auto"){
#   
#   for (j in 1:length(noise.path)) {
#     if(exists("df2")) rm(df2)
#     try(df2<-read.csv(noise.path[j], sep = "\t", header = FALSE) , silent = TRUE)
#    if(exists("df2")){
#     
#     poi<-c(poi, df2$V1[which(((df2$V2-median(df2$V2))/sd(df2$V2))>2 & df2$V3>50)])
#     }
#   }
# }

if(poicurrent=="auto"){
  for (j in 1:length(noise.path)) {
    if(exists("df2")) rm(df2)
    try(df2<-read.csv(noise.path[j], sep = "\t", header = FALSE) , silent = TRUE)
    if(exists("df2")){
      
      den<-density(df2$V2)
      bins<-den$y
      border2<-bins[c(2*c(1:250))]
      border1<-bins[c(2*c(1:250))-1]
      border1<-as.data.frame(border1)
      border1$border2<-border2
      border1$mean<-apply(border1, 1, mean)
      border1$sd<-apply(border1, 1, sd)
      border1$cov<-border1$sd/border1$mean
      borderco<-den$x[min(which(border1$cov>0.5 & as.numeric(rownames(border1))>15))]
      poi<-c(poi, df2$V1[which(df2$V2 > borderco & df2$V3>50)])
    }
  }
}

poi<-unique(poi)


bam.path.total<-list.files(path, pattern = ".*bam$",full.names = TRUE, recursive = TRUE)

sizes<-file.info(bam.path.total)$size
if(length(which(sizes<5000))>0) bam.path.total<-bam.path.total[-which(sizes<5000)]


cores<-as.numeric(detectCores())-2

allpos<-TRUE

for (b in 1:length(bam.path.total)) {
  bam.path<-bam.path.total[b]

if(length(poi)>0){
  
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
                 " -f ", bam.path, "> ", gsub(".*/",path,gsub("\\.sorted.bam","_P",bam.path)), poi[k],".tsv",sep = ""))
  }
  
  stopCluster(cluster.cores)
  
}
  
  print(paste("Extracting Single Mutations for",gsub(".*/","",bam.path.total[b])) )
  position.files<-list.files(path, full.names = TRUE, pattern = "_P.*\\.tsv")
  basediff<-as.data.frame(matrix(NA, ncol = 7, nrow=length(position.files)))
  colnames(basediff)<-c("Base", "A","T","C","G", "Count","File")
  for (i in 1:length(position.files)) {
    dummy<-fread(position.files[i],sep = "\t", header = FALSE)
    
    if(nrow(dummy)>0){
      dummy<-as.data.frame(dummy[,-1])
      colnames(dummy)<-c("Position","Base")
      
      if(length(which(dummy$Base=="D"))>0) dummy<-dummy[-which(dummy$Base=="D"),]
      if(length(which(dummy$Base=="I"))>0) dummy<-dummy[-which(dummy$Base=="I"),]
    }
    if(nrow(dummy)>0){
      dummy<-as.data.frame(dummy)
      basediff$Base[i]<-unique(dummy$Position)[1]
      dummy<-as.data.frame(table(dummy$Base))
      
      basediff$Count[i]<-sum(dummy$Freq)
      
      if(length(which(dummy$Var1=="A"))>0){
        basediff$A[i]<-dummy$Freq[which(dummy$Var1=="A")]/sum(dummy$Freq)
      }else{
        basediff$A[i]<-0
      }
      
      if(length(which(dummy$Var1=="T"))>0){
        basediff$T[i]<-dummy$Freq[which(dummy$Var1=="T")]/sum(dummy$Freq)
      }else{
        basediff$T[i]<-0
      }
      
      if(length(which(dummy$Var1=="C"))>0){
        basediff$C[i]<-dummy$Freq[which(dummy$Var1=="C")]/sum(dummy$Freq)
      }else{
        basediff$C[i]<-0
      }
      
      if(length(which(dummy$Var1=="G"))>0){
        basediff$G[i]<-dummy$Freq[which(dummy$Var1=="G")]/sum(dummy$Freq)
      }else{
        basediff$G[i]<-0
      }
      basediff$File[i]<-gsub("_P.*","",gsub(".*/","",position.files[i]))
    }
    
    
  }
  
  if(length(which(is.na(basediff$File)))>0) basediff<-basediff[-which(is.na(basediff$File)),]
  
  if(length(which(is.na(basediff$Base)))>0)basediff<-basediff[-which(is.na(basediff$Base))]
  basediff.melted<-melt(basediff, id=c("Base","Count","File"))
  s<-read.fasta(SpikeRef)
  s<-as.character(s$Spike)
  
  basediff.melted$CI<-NA
  basediff.melted$AA<-NA
  
  basediff.melted$CI<-qt(0.975,df=basediff.melted$Count-1)*(sd(basediff.melted$value,na.rm = TRUE))/sqrt(basediff.melted$Count)    
  
  basediff.melted<-basediff.melted[which(basediff.melted$CI<basediff.melted$value),]
  
  bases.to.check<-unique(paste(basediff.melted$Base,basediff.melted$variable,sep = "-"))
  original<-translate(s)
  
  noise.table$Cover<-"High"
  noise.table$Cover[which(noise.table$V3<10)]<-"Low"
  noise.table<-noise.table[which(noise.table$V1 %in% as.numeric(gsub("-.*","",bases.to.check))) ,]
  
  cov.agg<-aggregate(V1~Sample+Cover, noise.table, length)
  cov.agg<-cov.agg[which(cov.agg$Cover=="Low"),]
  if(nrow(cov.agg)>0) cov.agg$Warning<-paste("Warn:",cov.agg$V1, "Positions with low coverage")
  
  for (i in 1:length(bases.to.check)){
    mutated<-s
    mutated[as.numeric(gsub("-.*","",bases.to.check[i]))]<-tolower(gsub(".*-","",bases.to.check[i]))
    mutated<-translate(mutated)
    if(length(which(original!=mutated))){
      basediff.melted$AA[which(paste(basediff.melted$Base,basediff.melted$variable,sep = "-") == bases.to.check[i])]<- 
        paste("S:",original[which(original!=mutated)], which(original!=mutated), mutated[which(original!=mutated)],sep = "")
    }else{
      basediff.melted$AA[which(paste(basediff.melted$Base,basediff.melted$variable,sep = "-") == bases.to.check[i])]<-"SynMut"
    }
    
  }
  
  if(!exists("basediff.melted.bk")){
    basediff.melted.bk<-basediff.melted
  }else{
    basediff.melted.bk<-rbind(basediff.melted.bk,basediff.melted)
  }
  
  file.remove(position.files)
  
  
}

basediff.melted<-basediff.melted.bk

  
if(length(which(basediff.melted$AA=="SynMut"))>0) basediff.melted<-basediff.melted[-which(basediff.melted$AA=="SynMut"),]
if(length(grep("\\*", basediff.melted$AA))>0) basediff.melted<-basediff.melted[-grep("\\*", basediff.melted$AA),]

basediff.melted$AApos<-gsub(".$","",basediff.melted$AA)
basediff.melted$AAnew<-gsub("S:.[0-9]","",basediff.melted$AA)
basediff.melted$AAnew<-gsub("^[0-9]","",basediff.melted$AAnew)
basediff.melted$AAnew<-gsub("^[0-9]","",basediff.melted$AAnew)
basediff.melted$AAnew<-gsub("^[0-9]","",basediff.melted$AAnew)


 

basediff.melted$AApos<-factor( basediff.melted$AApos, levels =unique(basediff.melted$AApos[order(as.numeric(gsub(".$","",gsub("S:.","",basediff.melted$AApos))))]) )


colnames(basediff.melted)[9]<-"Mutation"
if(nrow(cov.agg)>0){
  basediff.melted<-merge(basediff.melted, cov.agg, by.x="File", by.y="Sample", all.x=TRUE)
}else{
    basediff.melted$Warning<-NA
  }

basediff.melted$File[which(!is.na(basediff.melted$Warning))]<-paste(basediff.melted$File[which(!is.na(basediff.melted$Warning))], "/" ,basediff.melted$Warning[which(!is.na(basediff.melted$Warning))])

if(length(unique(basediff.melted$File))<10){
ggplot(basediff.melted)+
  geom_bar(aes(AApos,value, fill=Mutation), stat = "identity")+
  theme_minimal()+
  facet_wrap(~File,ncol = 1)+
    scale_fill_manual(values = rainbow(length(unique(basediff.melted$Mutation))))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  xlab("Position Spike")+
    ylab("Ratio mutant reads")+
    ylim(0,1)
  ggsave(paste(path,"SingleMutationsBarplots.pdf",sep = "") ,width =  8.27, height =  11.69)
}

if(length(unique(basediff.melted$File))>10){
  files.to.plot<-unique(basediff.melted$File)
printing<-TRUE
counter<-0
start<-1
stop<-10
  while(printing){
    counter<-counter+1
    if(stop>=length(files.to.plot)){
      stop<-length(files.to.plot)
      printing<-FALSE
    }
    ggplot(basediff.melted[which(basediff.melted$File %in% files.to.plot[start:stop]),])+
      geom_bar(aes(AApos,value, fill=Mutation), stat = "identity")+
      theme_minimal()+
      facet_wrap(~File,ncol = 1)+
      scale_fill_manual(values = rainbow(length(unique(basediff.melted$Mutation))))+
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
      xlab("Position Spike")+
      ylab("Ratio mutant reads")+
      ylim(0,1)
    
    ggsave(paste(path,"SingleMutationsBarplots_part",counter, ".pdf",sep = "") ,width =  8.27, height =  11.69)
    start<-stop+1
    stop<-start+9
    if(start>length(files.to.plot)) printing<-FALSE
}
  
}

file.remove(paste(path,"bbasereader",sep = ""))
write_xlsx(basediff.melted.bk, paste(path,"SingleMutationsResults.xlsx",sep = ""))     
