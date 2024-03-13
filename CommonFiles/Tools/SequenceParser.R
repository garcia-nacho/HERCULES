library(seqinr)
library(digest)

fastas<-seqinr::read.fasta("references.aligned.fasta")
metadata<-read.csv("references.csv",sep = ";", stringsAsFactors = FALSE)

metadata<-metadata[which(metadata$totalMissing<150),]
metadata$seqName<-gsub(" .*","",metadata$seqName)

if(length(which(duplicated(metadata$seqName)))>0) metadata<-metadata[-which(duplicated(metadata$seqName)),]
if(length(which(duplicated(names(fastas))))>0) fastas<-fastas[-which(duplicated(names(fastas)))]
fastas<-fastas[which(names(fastas) %in% metadata$seqName)]

fastas.trimmed<-lapply(fastas, function(x) x[c(21563:25384)] )
metadata$MashID<-NA
for (i in 1:length(fastas.trimmed)) {
  if(length(which(fastas.trimmed[[i]]=="-"))>0) fastas.trimmed[[i]] <- fastas.trimmed[[i]][-which(fastas.trimmed[[i]]=="-")] 
  names(fastas.trimmed)[i]<-paste("SID",i,"_",as.character(metadata$Nextclade_pango[which(metadata$seqName==names(fastas.trimmed)[i])]),sep = "")
  metadata$MashID[i]<-paste("SID",i,"_",as.character(metadata$Nextclade_pango[i]),sep = "")
  
}

#CleanDups
hashes<-lapply(fastas.trimmed, function(x) digest(paste(x,collapse = ""),algo="md5"))
hashes<-unlist(hashes)
dups<-names(hashes)[which( duplicated(hashes))]
fastas.trimmed<-fastas.trimmed[-which(names(fastas.trimmed) %in% dups )]
metadata<-metadata[which(metadata$MashID %in% names(fastas.trimmed)),]

write.fasta(fastas.trimmed,"/media/nacho/Data/wastewater/AFMode/refs/cleanrefs10012023.fasta", names = names(fastas.trimmed))

for (i in 1:length(fastas.trimmed)) {
  write.fasta(fastas.trimmed[[i]], paste("/media/nacho/Data/wastewater/AFMode/seqs/",names(fastas.trimmed)[i],".fasta",sep = ""),
              names = names(fastas.trimmed)[i])
}


#UMAP mapping. 

fastas<-seqinr::read.fasta("references.aligned.fasta")
fastas.umap<-lapply(fastas, function(x) x[c((21563+1250):(21563+2250))] )
fastas.umap<-unlist(lapply(fastas.umap, function(x)paste(as.character(x),collapse = "")))
metadata<-read.csv("refs/references.csv",sep = ";", stringsAsFactors = FALSE)
metadata$seqName<-gsub(" .*","",metadata$seqName)
umap.list<-list()
pb<-txtProgressBar(min = 1, max = length(fastas.umap), initial = 1)
for (i in 1:length(fastas.umap)) {
  setTxtProgressBar(pb,i)
  templist<-strsplit(system(paste('printf ">Read\n',fastas.umap[i],'\n" | ./mash dist ./reference10012023.msh -',sep = "" ), intern = TRUE),"\t") 
  templist<-as.data.frame(do.call("rbind", templist))
  umap.list[[i]]<-templist$V3
}

out<-do.call("rbind", umap.list)
out<-as.data.frame(out)
write.csv(out, "MatrixForUMAP.csv", row.names = FALSE)
library(Rtsne)
library(umap)
library(ggplot2)


custom.config <- umap.defaults
custom.config$min_dist<-0.2
custom.config$n_neighbors<-20
out<-apply(out,2,as.numeric)

umap.res<-umap(as.matrix(out), config = custom.config)
umap.mapping<-as.data.frame(umap.res$layout)
umap.mapping$Sample<-metadata$Nextclade_pango
umap.mapping$Clade<-metadata$clade
umap.mapping$WhoClade<-metadata$clade_who

ggplot(umap.mapping)+
  geom_jitter(aes(V1, V2, col=Sample),alpha=0.01)+
  theme_minimal()+
  theme(legend.position = "none")

ggplot(umap.mapping)+
  geom_jitter(aes(V1, V2, col=Clade),alpha=0.5)+
  scale_color_manual(values=rainbow(length(unique(umap.mapping$Clade))))+
  theme_minimal()

ggplot(umap.mapping)+
  geom_jitter(aes(V1, V2, col=WhoClade),alpha=0.5)+
  scale_color_manual(values=rainbow(length(unique(umap.mapping$WhoClade))))+
  theme_minimal()

new<-predict(umap.res, out[c(1:10000),])

saveRDS(umap.res, "UmapModel.rds")
umap.res<-readRDS("UmapModel.rds")

new<-predict(umap.res, out[rep(c(1:10000),10),])
