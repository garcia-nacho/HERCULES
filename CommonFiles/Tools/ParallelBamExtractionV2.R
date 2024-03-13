library(data.table)
library("doParallel")
library("parallel")
library("foreach")
library(doSNOW)
library(progress)
library(seqinr)
library(umap)
library(ggplot2)
library(digest)
library(plotly)
library(htmlwidgets)
library(writexl)
library(ggsankey)

#src/bam2msa /media/nacho/Data/DockerImages/Wastewater_SARS-CoV-2/CommonFiles/reference/SpikeRef.fa ../fullwaste18/results/bam/waste018.43_Oslo.sorted.bam Spike:1260-2000 bam2msa
arg=commandArgs(TRUE)
#1: Cutoff 0.1
#2: start 1250
#3: end 2250
#4: poicurrent "auto"
co.n<-as.numeric(arg[1])
n.start<-as.numeric(arg[2])
n.end<- as.numeric(arg[3])
poicurrent<- arg[4]

# co.n<-0.1
# n.start<-1250
# n.end<-2250
# poicurrent<-"G22895C, T22896A, G22898A, A22910G, C22916T,G23012A, C23013A, T23018C, T23019C, C23271T, C23423T, A23604G"
#1355	28548	waste011.40_Bergen	G	0,718649292	0,004937526	S:L452R
#1841	26601	waste011.40_Bergen	G	0,986955378	0,005115047	S:D614G
#2037	28145	waste011.40_Bergen	G	0,983798188	0,004972753	S:N679K
#1320	12470	waste011.40_Oslo	G	0,973857257	0,007451126	S:N440K
#1513	11999	waste016.34_Oslo	C	0,875239603	0,00759624	S:Y505H

# pmatrixfile<-"/media/nacho/Data/DockerImages/Wastewater_SARS-CoV-2/CommonFiles/reference/ProbMatrix.csv"
# reference<-read.fasta("/media/nacho/Data/DockerImages/Wastewater_SARS-CoV-2/CommonFiles/reference/SpikeRef.fa")
# file.copy("/media/nacho/Data/DockerImages/Wastewater_SARS-CoV-2/CommonFiles/Tools/bam2msa","bam2msa")
# file.copy("/media/nacho/Data/DockerImages/Wastewater_SARS-CoV-2/CommonFiles/reference/SpikeRef.fa","SpikeRef.fa")
# refmsa<-"/media/nacho/Data/DockerImages/Wastewater_SARS-CoV-2/CommonFiles/reference/MSA_Refs.tsv.gz"
# refmsaid<-"/media/nacho/Data/DockerImages/Wastewater_SARS-CoV-2/CommonFiles/reference/MSA_RefsID.csv"
# hasshes<-"/media/nacho/Data/DockerImages/Wastewater_SARS-CoV-2/CommonFiles/reference/variant_hash.csv"

pmatrixfile<-"/home/docker/CommonFiles/reference/ProbMatrix.csv"
reference<-read.fasta("/home/docker/CommonFiles/reference/SpikeRef.fa")
file.copy("/home/docker/CommonFiles/Tools/bam2msa","bam2msa")
file.copy("/home/docker/CommonFiles/reference/SpikeRef.fa","SpikeRef.fa")
refmsa<-"/home/docker/CommonFiles/reference/MSA_Refs.tsv.gz"
refmsaid<-"/home/docker/CommonFiles/reference/MSA_RefsID.csv"
hasshes<-"/home/docker/CommonFiles/reference/variant_hash.csv"

reference<-unlist(reference)

# Parallel generation of fastas -------------------------------------------

if(!dir.exists("MSAFastas")){
  print("Extracting bams")
bam.files<-list.files("bam",pattern = ".bam$",full.names = TRUE)
dir.create("MSAFastas")
samples.to.analyze<-c(1:length(bam.files))

pb <- progress_bar$new(
  format = "File: :samp.pb [:bar] :elapsed | eta: :eta",
  total = length(samples.to.analyze),    # 100 
  width = 60)
samp <- samples.to.analyze

progress <- function(n){
  pb$tick(tokens = list(samp.pb = samp[n]))
} 

opts <- list(progress = progress)

gc()
cores.n<-detectCores()
cores<- cores.n -2

cluster.cores<-makeCluster(cores)
registerDoSNOW(cluster.cores)

out.par<-foreach(i=1:length(samples.to.analyze), .verbose=FALSE, .options.snow = opts) %dopar%{
  system(paste("./bam2msa ", "./SpikeRef.fa ",bam.files[i]," Spike:1250-2250",
             " | cut -f 1-2 > ", gsub(".*/","MSAFastas/",gsub(".bam","_b2f.tsv",bam.files[i])), sep = ""))
  system(paste("gzip ",gsub(".*/","MSAFastas/",gsub(".bam","_b2f.tsv",bam.files[i])), sep = ""))
  
}
stopCluster(cluster.cores)
}

# Noise co calculation -------------------------------------------------------
print("Calculating noise")
file.copy(refmsa, "MSAFastas/Reference.b2f.tsv.gz")

if(poicurrent!="0" & poicurrent!="auto"){
  poicurrent<-(gsub(" ","",unlist(base::strsplit(poicurrent, ","))))
  positions<-toupper(poicurrent)
  
  positions<-gsub("A|T|C|G", "", positions)
  positions<-as.numeric(positions)
  if(length(which(is.na(positions)))>0){
    positions<-positions[-which(is.na(positions))]
  }
  
  positions<-positions[order(positions)]
  
  if(positions[1]>21563 & positions[1]<25384 ){
    positions<- positions-21563+1
  }else if(positions[1]>1 & positions[1] <38321){
    positions<-positions
  }else{
    positions<-"auto"
  }
  poicurrent<-positions
} 

noise.path<-list.files(pattern = ".*noise.tsv",full.names = TRUE, recursive = TRUE)

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


noise.table.agg<-aggregate(V2~V1, noise.table, mean)
depth.table.agg<-aggregate(V3~V1, noise.table, mean)

poi.co<-length(depth.table.agg$V3>20)*0.025

poit<-noise.table.agg[which(noise.table.agg$V2>co.n & noise.table.agg$V1>n.start & noise.table.agg$V1 < n.end & depth.table.agg$V3>20),]
poi<-poit$V1[order(poit$V2, decreasing=TRUE)][1:min(round(poi.co),nrow(poit))]
if(poicurrent[1]!="0")poi<-unique(c(poi,poicurrent))
if(length(which(poi>n.end))>0) poi<-poi[-which(poi>n.end)]
if(length(which(poi<n.start))>0) poi<-poi[-which(poi<n.start)]

if(poicurrent[1]=="auto"){
  poit<-noise.table.agg[which(noise.table.agg$V2>0 & noise.table.agg$V1>n.start & noise.table.agg$V1 < n.end & depth.table.agg$V3>20),]
  poi<- c(poi, poit$V1[order(poit$V2, decreasing=TRUE)][1:min(round(poi.co),nrow(poit))])
}

poi<-unique(c(poi,1841))

compressed<-list.files(full.names = TRUE, pattern = "b2f.tsv.gz",recursive = TRUE)
poi.full<-as.numeric(poi)
poi<-as.numeric(poi)-n.start
poi<-poi[order(poi)]

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

# Parallel running --------------------------------------------------------

cores<-as.numeric(detectCores())-2
cluster.cores<-makeCluster(cores)
registerDoSNOW(cluster.cores)

out.par<-foreach(cf=1:length(compressed), .verbose=FALSE, .packages = c("data.table","seqinr"), .options.snow = opts) %dopar%{

  target.file<-gsub(".gz","",compressed[cf])
  if(file.exists(target.file)) file.remove(target.file)
  system(paste("gunzip -c ", compressed[cf], " > ",target.file,sep = ""))
  
  df<-fread(target.file, sep = "\t")
  df<-as.data.frame(df)
  
  lineages<-list()
  lineages.ohe<-list()
  if(nrow(df)>1){
  for (i in 1:nrow(df)) {
    ins<-which(unlist(strsplit(df$ref_msa[i], ""))=="-")
    if(length(ins)>0 ){
      dum.out<-unlist(strsplit(df$`#query_msa`[i], ""))[-ins][poi+1]
    }else{
      dum.out<-unlist(strsplit(df$`#query_msa`[i], ""))[poi+1]
    }
    out.ohe<-rep(0,(length(poi)*5))
    out.ohe[match(dum.out, c("A","T","C","G","-"))+seq(from=0, to=(length(poi)-1)*5, by=5)]<-1
    
    lineages[[i]]<-dum.out
    lineages.ohe[[i]]<-out.ohe
  }
  
  file.remove(target.file)
  
  lineages<-do.call( rbind, lineages)
  lineages.ohe<-do.call( rbind, lineages.ohe)
  lineages<-as.data.frame(lineages)
  lineages.ohe<-as.data.frame(lineages.ohe)
  
  colnames(lineages)<-paste("P",poi,sep = "")
  colnames(lineages.ohe)<-paste(rep(poi+1250,5)[order(rep(poi+1250,5))],c("A","T","C","G","D"),sep = "")
  

  # Pango Lineage probability calculation -----------------------------------------
  pmatrix<-read.csv(pmatrixfile) 
  rownames(pmatrix)<-pmatrix$X
  pmatrix$X<-NULL
  pmatrix<-pmatrix[which(row.names(pmatrix) %in% colnames(lineages.ohe)),]
  
  if(length(which(duplicated(lineages.ohe)))>0){
    lin.ohe.short<-lineages.ohe[-which(duplicated(lineages.ohe)),]
  }else{
    lin.ohe.short<-lineages.ohe
}
lin.mut.aa<-vector()
for (i in 1:nrow(lin.ohe.short)) {
  mlohe<- colnames(lin.ohe.short)[which(lin.ohe.short[i,]==1)]
  mut.ref.ohe<-reference
  mut.ref.ohe[as.numeric(gsub(".$","",mlohe))]<-tolower(gsub("[0-9]+","",mlohe))
  ref.trans<-translate(reference)
  mut.trans<-translate(mut.ref.ohe)
  
  if(length(which(ref.trans!=mut.trans))>0){
    mut.dum<-paste("S:", ref.trans[which(ref.trans!=mut.trans)],which(ref.trans!=mut.trans), mut.trans[which(ref.trans!=mut.trans)],sep = ""  )
    lin.mut.aa<-c(lin.mut.aa,paste(mut.dum,collapse = "/"))  
  }else{
    lin.mut.aa<-c(lin.mut.aa, NA)
  }
}

mut.aa<-vector()
for (i in 1:nrow(pmatrix)) {
  mut.ref.ohe<-reference
  mut.ref.ohe[as.numeric(gsub(".$","",row.names(pmatrix)[i]))]<-tolower(gsub("[0-9]+","",row.names(pmatrix)[i]))
  ref.trans<-translate(reference)
  mut.trans<-translate(mut.ref.ohe)
  
  if(length(which(ref.trans!=mut.trans))>0){
    mut.dum<-paste("S:", ref.trans[which(ref.trans!=mut.trans)],which(ref.trans!=mut.trans), mut.trans[which(ref.trans!=mut.trans)],sep = ""  )
    mut.aa<-c(mut.aa,mut.dum)  
  }else{
    mut.aa<-c(mut.aa, NA)
  }
  
}

lin.ohe.short$Bin<-apply(lin.ohe.short, 1,function(x) paste(x, collapse = ""))
lin.ohe.short$Mut<-lin.mut.aa
lin.ohe.short<-lin.ohe.short[,which(colnames(lin.ohe.short) %in% c("Bin", "Mut"))]

total.lin<-apply(lineages.ohe, 1,function(x) paste(x, collapse = ""))
total.lin.aa<-lin.ohe.short$Mut[match(total.lin, lin.ohe.short$Bin)]
ToMakeNA<-which(is.na(pmatrix),arr.ind = TRUE)
pmatrix[which(is.na(pmatrix), arr.ind = TRUE)]<-0
plinmatrix<-as.matrix(lineages.ohe) %*% as.matrix(pmatrix)

# plinvec<-as.numeric(apply(pmatrix, 2,sum))
# #plinmatrix2<-plinmatrix
# for (i in 1:nrow(plinmatrix)) {
#   plinmatrix[i,]<-as.numeric(plinmatrix[i,])/plinvec
# }

anti.pmatrix<-pmatrix-1

#Deletions do not count as negative scoring factors
anti.pmatrix[ToMakeNA]<-0

anti.plinmatrix<-as.matrix(lineages.ohe) %*% as.matrix(anti.pmatrix)

plinmatrix<-plinmatrix+anti.plinmatrix
  
  lin.index<-apply(plinmatrix, 1, function(x)which(x==max(x,na.rm = TRUE)))
  lin.max<-apply(plinmatrix, 1, function(x)max(x,na.rm = TRUE))
  lin.tab<-unlist(lapply(lin.index, function(x) paste(colnames(plinmatrix)[x],collapse = "/")))
  lin.tab.raw<-lin.tab
  lin.tab<-as.data.frame(table(lin.tab))
  
  
  colnames(lin.tab)<-c("Lineage", "Count")
  lin.tab$Ratio<-lin.tab$Count/sum(lin.tab$Count)
  lin.tab$Sample<- gsub("\\.sorted.*","",gsub(".*/","",target.file))

  #To return
  output.pango<-lin.tab

  #Lineages compact 
  
  lineages.compact<- as.data.frame(apply(lineages,1, function(x)paste(x, collapse = "")))
  colnames(lineages.compact)<-"Seq"
  
  lineages.compact$PangoLineage<-lin.tab.raw
  lineages.compact$Count<-0
  lineages.compact$Score<-lin.max
  lineages.compact$Mut.aa<-total.lin.aa
  lineages.compact.agg<-aggregate(Count ~Seq, lineages.compact, length)
  
  lineages.compact$Count<-NULL
  if(length(which(duplicated(lineages.compact$Seq)))>0){
    lineages.compact.unique<-lineages.compact[-which(duplicated(lineages.compact$Seq)),]
  }else{
    lineages.compact.unique<-lineages.compact
    }
  lineages.compact.agg<-merge(lineages.compact.agg, lineages.compact.unique, by="Seq")

  lineages.compact.agg$Positions<-paste(gsub("P","",colnames(lineages)),collapse = "/")
  
  lineages.compact.agg<-lineages.compact.agg[which(lineages.compact.agg$Count>0.001*nrow(lineages)),]
  
  lin.clean<-aggregate(Count~PangoLineage, lineages.compact.agg, sum)
  
  colnames(lin.clean)<-c("Lineage", "Count")
  lin.clean$Ratio<-lin.clean$Count/sum(lin.clean$Count)
  lin.clean$Sample<- gsub("\\.sorted.*","",gsub(".*/","",target.file))
  
  #To return
  output.pango.clean<-lin.clean

  plin.dum<-as.data.frame(plinmatrix)
  
  plin.dum$Seq<-lineages.compact$Seq
  if(length(which(duplicated(plin.dum)))>0) plin.dum<-plin.dum[-which(duplicated(plin.dum)),]
  
  index<-vector()
  for (i in 1:nrow(lineages.compact.agg)) {
    index<-c(index,which(lineages.compact$Seq==lineages.compact.agg$Seq[i])[1])
  }
  
  lineages.compact.agg$Sample<-gsub("\\.sorted.*","",gsub(".*/","",target.file))
  lineages.compact.agg$Ratio<-lineages.compact.agg$Count/sum(lineages.compact.agg$Count)
  
  #To return
  lineages.compact.out<-lineages.compact.agg
  
  plin.dum<-as.data.frame(plinmatrix[index,,drop=FALSE])
  plin.dum$Seq<-lineages.compact.agg$Seq
  plin.dum$Sample<-lineages.compact.agg$Sample
  
  #To return
  plin.out<-plin.dum

  # Single Mutation Calculation ----------------------------------------------

  
  colnames(lineages.ohe)<-mut.aa
  if(length(which(is.na(mut.aa)))>0)  lineages.ohe<-lineages.ohe[,-which(is.na(mut.aa))]
  
  sin.mut<-as.data.frame(colnames(lineages.ohe)) 
  colnames(sin.mut)<-"Mutation"
  sin.mut$Count<-NA
  sin.mut$Ratio<-NA
  
  for (i in 1:nrow(sin.mut)) {
    sin.mut$Count[i]<-sum(lineages.ohe[,which(colnames(lineages.ohe)==sin.mut$Mutation[i])])
    sin.mut$Ratio[i]<-sum(lineages.ohe[,which(colnames(lineages.ohe)==sin.mut$Mutation[i])])/nrow(lineages.ohe)
  }
  sin.mut$Sample<-gsub("\\.sorted.*","",gsub(".*/","",target.file))
  
  #To return
  mut.out<-sin.mut
  rm(plin.dum,df,lin.tab,pmatrix,total.lin,sin.mut,  lineages.compact.agg, lineages.ohe, anti.pmatrix, anti.plinmatrix)
  rtnlist<-list(output.pango,output.pango.clean,lineages.compact.out,plin.out,mut.out)
  rm(output.pango,output.pango.clean,lineages.compact.out,plin.out,mut.out)
  gc()
  names(rtnlist)<-c("Pango","PangoClean","Lineages","Mtrx","Mutations")
  return(rtnlist)
  }else{
    file.remove(target.file)
  }
}

stopCluster(cluster.cores)
gc()


Dlist<-lapply(out.par, function(x)x[1][[1]])
miss<-lapply(Dlist, is.data.frame)
if(length( which(!unlist(miss)))>0) Dlist<-Dlist[-which(!unlist(miss))]
pango.df<-do.call(rbind, Dlist)

Dlist<-lapply(out.par, function(x)x[2][[1]])
miss<-lapply(Dlist, is.data.frame)
if(length( which(!unlist(miss)))>0) Dlist<-Dlist[-which(!unlist(miss))]
pangoClean.df<-do.call(rbind, Dlist)

Dlist<-lapply(out.par, function(x)x[3][[1]])
miss<-lapply(Dlist, is.data.frame)
if(length( which(!unlist(miss)))>0) Dlist<-Dlist[-which(!unlist(miss))]
lineages.df<-do.call(rbind, Dlist)

Dlist<-lapply(out.par, function(x)x[4][[1]])
if(length( which(!unlist(miss)))>0) Dlist<-Dlist[-which(!unlist(miss))]
probMtrx<-do.call(rbind, Dlist)

Dlist<-lapply(out.par, function(x)x[5][[1]])
if(length( which(!unlist(miss)))>0) Dlist<-Dlist[-which(!unlist(miss))]
Mutation.df<-do.call(rbind, Dlist)


rm(out.par)
#Reference asignment method2

refids<-read.csv(refmsaid)

system("gunzip -c MSAFastas/Reference.b2f.tsv.gz > MSAFastas/Refs.tsv")

df<-fread("MSAFastas/Refs.tsv", sep = "\t")
df<-as.data.frame(df)

lineagesref<-list()

  for (i in 1:nrow(df)) {
    ins<-which(unlist(strsplit(df$ref_msa[i], ""))=="-")
    if(length(ins)>0 ){
      dum.out<-unlist(strsplit(df$`#query_msa`[i], ""))[-ins][poi+1]
    }else{
      dum.out<-unlist(strsplit(df$`#query_msa`[i], ""))[poi+1]
    }

    lineagesref[[i]]<-dum.out

  }
  
lineagesref<-lapply(lineagesref,function(x)paste(x, collapse = ""))
lineagesref<-as.data.frame(unlist(lineagesref))
colnames(lineagesref)<-"Sequence"
lineagesref$ID<-refids$x

file.remove("MSAFastas/Refs.tsv")
file.remove("MSAFastas/Reference.b2f.tsv.gz")



lineages.df.ref<-lineages.df[which(lineages.df$Sample=="Reference.b2f.tsv"),]
probMtrx.ref<-probMtrx[which(lineages.df$Sample=="Reference.b2f.tsv"),]

probMtrx<-probMtrx[-which(lineages.df$Sample=="Reference.b2f.tsv"),]
lineages.df<-lineages.df[-which(lineages.df$Sample=="Reference.b2f.tsv"),]

lineages.df.ref$ID<-lineagesref$ID[match(lineages.df.ref$Seq, lineagesref$Sequence)]

linx<-list()
for (i in 1:nrow(lineages.df.ref)) {
linx[[i]]<-lineages.df.ref$ID[which(lineages.df.ref$Seq==lineages.df.ref$Seq[i])]  
}

lineages.df.ref$ID2<-unlist(lapply(linx, function(x)paste(gsub("_SQ.*","",x),collapse = "/")))

lineages.df$PangoLineagesMatched<-lineages.df.ref$ID2[match(lineages.df$PangoLineage, lineages.df.ref$PangoLineage)]

if(length(which(is.na(lineages.df$PangoLineagesMatched)))>0){
  lineages.df$PangoLineagesMatched[which(is.na(lineages.df$PangoLineagesMatched))]<-"Unclassified"  
}


#Matching 


# UMAP --------------------------------------------------------------------


print("Dimension reduction...")
#umapping<-umap(as.matrix(probMtrx[which(lineages.df$Count>10),-c(which(colnames(probMtrx) %in% c("Seq","Sample")))]))
#lineages.clean<-lineages.df[which(lineages.df$Count>10),]

umapping<-uwot::umap(as.matrix(probMtrx[-c(which(colnames(probMtrx) %in% c("Seq","Sample")))]))
  
lineages.clean<-lineages.df

lineages.clean$X<-umapping[,1]
lineages.clean$Y<-umapping[,2]

#lineages.clean$X<-umapping$layout[,1]
#lineages.clean$Y<-umapping$layout[,2]

variant.hash<-read.csv(hasshes)

if(length(grep("NA\\..*\\.X", lineages.clean$PangoLineagesMatched))>0){
  to.match<-unique(lineages.clean$PangoLineagesMatched[grep("NA\\..*\\.X", lineages.clean$PangoLineagesMatched)])  
}

for(tm in c(1:length(to.match))){
  variant.temp<-variant.hash[which(variant.hash$Lineage==to.match[tm]),]
  variant.temp<-variant.temp[order(variant.temp$Freq, decreasing = TRUE),]
  variant.temp<-variant.temp[c(1:min(3,nrow(variant.temp))),]
  variant.temp$NewID<-paste(variant.temp$pangos,":" ,round((variant.temp$Freq)*100,0),"%" ,sep = "")
  newid<-paste(variant.temp$NewID, collapse = "/")
  lineages.clean$PangoLineagesMatched[which(lineages.clean$PangoLineagesMatched==to.match[tm])]<-newid
}


samples.vec<-unique(lineages.clean$Sample)
experiments<-unique(gsub("\\..*","",samples.vec))

intersetors<-intersect(c(1:nrow(lineages.clean))[-grep("/",lineages.clean$PangoLineage)], which(lineages.clean$PangoLineagesMatched=="Unclassified"))

if(length(intersetors)>0) lineages.clean$PangoLineagesMatched[intersetors]<-lineages.clean$PangoLineage[intersetors]


for (i in 1:length(experiments)) {
ggpl<-  ggplot(lineages.clean[grep(experiments[i], lineages.clean$Sample),])+
    geom_point(aes(X,Y,col=PangoLineagesMatched,size=Count,label=Mut.aa),alpha=0.3)+
    scale_colour_manual(values = rainbow(length(unique(lineages.clean$PangoLineagesMatched))))+
    theme_minimal()+
    xlab("Lineage Representation Dim1")+
    ylab("Lineage Representation Dim2")+
    facet_wrap(~Sample)

plty<-ggplotly(ggpl,tooltip = c("PangoLineagesMatched", "Count", "Mut.aa"))

saveWidget(partial_bundle(plty), paste("WWPlot_",experiments[i],"_lineages.html",sep = ""))
}

lineages.clean$Location<-gsub(".*_","",lineages.clean$Sample)
lineages.clean$Week<-gsub(".*\\.","",gsub("_.*","",lineages.clean$Sample))

uniqueaa<-lineages.clean[-which(duplicated(lineages.clean$Mut.aa)),]

ggplot(lineages.clean[-grep("Pos|Neg", lineages.clean$Sample),])+
  geom_point(aes(X,Y,col=Location,size=Ratio),alpha=0.3)+
  scale_colour_manual(values = rainbow(length(unique(lineages.clean$Location))))+
  theme_minimal()+
  facet_wrap(~as.numeric(Week))  


ptl2<- ggplotly(ggplot(lineages.clean[-grep("Pos|Neg", lineages.clean$Sample),])+
           geom_jitter(aes(X,Y,col=PangoLineagesMatched,size=Ratio, label=Mut.aa),alpha=0.3)+
           scale_colour_manual(values = rainbow(length(unique(lineages.clean$PangoLineagesMatched))))+
           theme_minimal()+
           facet_wrap(~Week+Location), tooltip = c("PangoLineagesMatched", "Count","Mut.aa") )

saveWidget(partial_bundle(ptl2), "TotalLineageMap.html")

ptl2<- ggplotly(ggplot(lineages.clean[-grep("Pos|Neg", lineages.clean$Sample),])+
                  geom_jitter(aes(X,Y,col=PangoLineagesMatched,size=Ratio, label=Mut.aa),alpha=0.3)+
                  scale_colour_manual(values = rainbow(length(unique(lineages.clean$PangoLineagesMatched))))+
                  theme_minimal()+
                  facet_wrap(~as.numeric(Week)), tooltip = c("PangoLineagesMatched", "Count","Mut.aa") )

saveWidget(partial_bundle(ptl2), "TotalLineageMap2.html")


weeks<-unique(lineages.clean$Week[-grep("Pos|Neg", lineages.clean$Sample)])
weeks<-weeks[order(as.numeric(weeks),decreasing = TRUE)]

# ptl3<- ggplotly(ggplot(lineages.clean[which(lineages.clean$Week %in% as.character(weeks[c(1:4)])),])+
#                   geom_jitter(aes(X,Y,col=PangoLineagesMatched,size=Count, label=Mut.aa),alpha=0.3)+
#                   scale_colour_manual(values = rainbow(length(unique(lineages.clean$PangoLineagesMatched))))+
#                   theme_minimal()+
#                   facet_wrap(~Week+Location), tooltip = c("PangoLineagesMatched", "Count", "Mut.aa"))

#saveWidget(partial_bundle(ptl3), "LastMonthLineageMap.html")



write_xlsx(lineages.clean[,c(6,1,8,2,9,4,5,7,12,13,10,11)], "LineagesClean.xlsx")

lineages.agg<-aggregate(Count~PangoLineagesMatched+ Sample, lineages.clean,sum)
total.agg<-aggregate(Count~Sample, lineages.clean,sum)
colnames(total.agg)[2]<-"TotalCount"
lineages.agg<-merge(lineages.agg, total.agg, by="Sample")
lineages.agg$Ratio<-lineages.agg$Count/lineages.agg$TotalCount

ggplot(lineages.agg)+
  geom_bar(aes( Sample, Ratio, fill=PangoLineagesMatched), stat = "identity")+
  theme_minimal()+
   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggsave("PangoLineageBarPlot.pdf" ,width = 20, height = 16)


lplt<-ggplotly(ggplot(lineages.agg)+
                 geom_bar(aes( Sample, Ratio, fill=PangoLineagesMatched), stat = "identity")+
                 theme_minimal()+
                 theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)))

saveWidget(partial_bundle(lplt), "PangoLineageBarPlot.html")
write_xlsx(lineages.agg, "PangoLineages.xlsx")

#Mutation Plot
Mutation.df<-Mutation.df[-which(Mutation.df$Sample=="Reference.b2f.tsv"),]
Mutation.df$Mutation<-gsub("\\.1$","",Mutation.df$Mutation)
if(length(grep("\\*$", Mutation.df$Mutation))>0) Mutation.df<-Mutation.df[-grep("\\*$", Mutation.df$Mutation),]
Mutation.df.agg<-aggregate(Count~Mutation+Sample, Mutation.df, sum)
Mutation.df.total<-aggregate(Count~Sample, pango.df, sum)
Mutation.df.total<-Mutation.df.total[-which(Mutation.df.total$Sample=="Reference.b2f.tsv"),]
colnames(Mutation.df.total)[2]<-"TotalCount"
Mutation.df.agg<-merge(Mutation.df.agg, Mutation.df.total, by="Sample")
Mutation.df.agg$Ratio<-Mutation.df.agg$Count/Mutation.df.agg$TotalCount
Mutation.df.agg$Position<-as.numeric(gsub("^.","",gsub(".$","",gsub("S:","",Mutation.df.agg$Mutation))))
Mutation.df.agg$AApos<-gsub(".$","",Mutation.df.agg$Mutation)
Mutation.df.agg$Mutation<-gsub(".*[0-9]+","",Mutation.df.agg$Mutation)

#Plotting single mutations
if(length(unique(Mutation.df.agg$Sample))<=10){
  ggplot(Mutation.df.agg)+
    geom_bar(aes(AApos,Ratio, fill=Mutation), stat = "identity")+
    theme_minimal()+
    facet_wrap(~Sample,ncol = 1)+
    scale_fill_manual(values = rainbow(length(unique(Mutation.df.agg$Mutation))))+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
    xlab("Position Spike")+
    ylab("Ratio mutant reads")+
    ylim(0,1)
  ggsave(paste(path,"SingleMutationsBarplots.pdf",sep = "") ,width =  8.27, height =  11.69)
}

if(length(unique(Mutation.df.agg$Sample))>10){
  files.to.plot<-unique(Mutation.df.agg$Sample)
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
    ggplot(Mutation.df.agg[which(Mutation.df.agg$Sample %in% files.to.plot[start:stop]),])+
      geom_bar(aes(AApos,Ratio, fill=Mutation), stat = "identity")+
      theme_minimal()+
      facet_wrap(~Sample,ncol = 1)+
      scale_fill_manual(values = rainbow(length(unique(Mutation.df.agg$Mutation))))+
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
      xlab("Position Spike")+
      ylab("Ratio mutant reads")+
      ylim(0,1)
    
    ggsave(paste("SingleMutationsBarplots_part",counter, ".pdf",sep = "") ,width =  8.27, height =  11.69)
    start<-stop+1
    stop<-start+9
    if(start>length(files.to.plot)) printing<-FALSE
  }
}

write_xlsx(Mutation.df.agg, "SingleMutationResults.xlsx")



# Sankeys -----------------------------------------------------------------


print("Printing sankeys")
samples<-unique(lineages.clean$Sample)
sk.list<-list()

for (sk in 1:length(samples)) {
  
  df<-lineages.clean[which(lineages.clean$Sample==samples[sk]),]
  df<-aggregate(Count~Mut.aa, df, sum)
  df<-df[order(df$Count, decreasing = TRUE),]
  
  if(exists("out"))try(rm(out))
  for (i in 1:min(nrow(df),20)) {
    if(exists("mutpos.df"))rm(mutpos.df)
    if(df$Mut.aa[i]!=""){
      mutations<-unlist(base::strsplit(df$Mut.aa[i], "/"))
      mutpos<- gsub("^.","",gsub(".$","",gsub("S:","",mutations)))
      
      mutpos.df<-as.data.frame(t(mutations))
      colnames(mutpos.df)<-mutpos
    }

    if(df$Mut.aa[i]=="" ){
      mutpos.df<-as.data.frame(matrix("Ref", nrow = 1, ncol = 1))
      colnames(mutpos.df)<-"Dummy"
    }
    
    if(exists("mutpos.df")){
      
      if(length(grep("-",mutpos.df[1,]))>0) mutpos.df<-mutpos.df[,-grep("-", mutpos.df[1,]),drop=FALSE]
      
      mutpos.df<-mutpos.df[round(rep(1,df$Count[i])),,drop(FALSE)]

      if(!exists("out")){
        out<-mutpos.df
      }else{
        
        if(length(which(colnames(mutpos.df) %in% colnames(out)))>0){
          missing.out<- colnames(mutpos.df)[-which(colnames(mutpos.df) %in% colnames(out))]
        }else{
          missing.out<- colnames(mutpos.df)
        }
        
        if(length(which(colnames(out) %in% colnames(mutpos.df)))>0){
          missing.df<- colnames(out)[-which(colnames(out) %in% colnames(mutpos.df))]
        }else{
          missing.df<- colnames(out)
        }
        
        if(length(missing.out)>0){
          pad.out<-as.data.frame(matrix(data = "Ref", nrow = nrow(out), ncol = length(missing.out)))
          colnames(pad.out)<-missing.out
          out<-cbind(out, pad.out)
        }
        
        if(length(missing.df)>0){
          pad.df<-as.data.frame(matrix(data = "Ref", nrow = nrow(mutpos.df), ncol = length(missing.df))) 
          colnames(pad.df)<-missing.df
          mutpos.df<- cbind(mutpos.df, pad.df)
        }
        
        out<-rbind(out, mutpos.df)  
      }
    }
  }
  
  try(out$Dummy<-NULL)
  
  out<-out[,-which(apply(out,2, function(x) length(unique(x))) ==1), drop=FALSE]
  
  if(ncol(out)>0){
    df.sankey <- make_long(out, colnames(out))
    
    ggplot(df.sankey, aes(x = x, 
                          next_x = next_x, 
                          node = node, 
                          next_node = next_node,
                          fill = factor(node))) +
      geom_sankey(flow.alpha = 0.5, node.color = 1) +
      scale_fill_viridis_d(option = "A", alpha = 0.5) +
      theme_sankey(base_size = 16)+
      xlab("Position")+
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
      ggtitle(paste("Mutation combinations in", samples[sk]))+
      labs(fill = "Mutations")
    
    ggsave(paste("SK_",samples[sk],"_Sankeyplot.Mutations.pdf"), width =   40, height = 20)
    
  }
}

file.remove("bam2msa")
file.remove("SpikeRef.fa")
