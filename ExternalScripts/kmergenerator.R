library(seqinr)

nc<-read.csv("/media/nacho/Data/wastewater/ReferencesProb/NextcladeTotal.csv")
ref<-read.fasta("/media/nacho/Data/DockerImages/FHISC2/FHI_SC2_Pipeline_Illumina/CommonFiles/reference_nc.fasta")
spikestart<-21563

pcrfw<-1250
pcrrv<-2250

nc<-nc[-which(duplicated(nc$seqName)),]


pangotab<-as.data.frame(table(nc$Nextclade_pango))

longpango<-nc$Nextclade_pango
longpangolist<-strsplit(longpango,"\\.")
sortpango<-unlist(lapply(longpangolist, function(x)paste(x[1:min(length(x),3)],collapse = ".") ))
nc$shortpango<-sortpango

shortpangotab<-as.data.frame(table(nc$shortpango))

kmerlist<-list()
pb<-txtProgressBar(max = nrow(shortpangotab))
for (i in 1:nrow(shortpangotab)) {
  setTxtProgressBar(pb,i)
  nc_dum<-nc[which(nc$shortpango==shortpangotab$Var1[i]),]
  mutations<-unlist(strsplit(nc_dum$substitutions,",") )
  mutations<-as.data.frame(table(mutations))
  mutations$ratio<-mutations$Freq/nrow(nc_dum)
  mutations<-mutations[which(mutations$ratio>0.9),]
  mutations$site<-as.numeric(gsub(".$","",gsub("^.","",mutations$mutations)))
  mutations<-mutations[which(mutations$site>spikestart+pcrfw & mutations$site<spikestart+pcrrv),]
  if(nrow(mutations)>2){

  combinations <- expand.grid(mutations$mutations, mutations$mutations)
  dups<-apply(combinations,1, function(x) max(table(x)) )
  combinations <- combinations[which(dups==1),]
  combinations$Var1<-as.character(combinations$Var1)
  combinations$Var2<-as.character(combinations$Var2)
  kmers<-vector()
  for (j in 1:nrow(combinations)) {
    muts<-as.character(combinations[j,])
    sites<-as.numeric(gsub(".$","",gsub("^.","",muts)))
    muts<-paste(muts[order(sites)], collapse = "-")
    kmers<-c(kmers, muts)
  }
  kmers<-kmers[-which(duplicated(kmers))] 
  kmerlist<-c(kmerlist, list(kmers))
  names(kmerlist)[length(kmerlist)]<-as.character(shortpangotab$Var1[i])
  }else{
    kmerlist<-c(kmerlist, list(NA))
    names(kmerlist)[length(kmerlist)]<-as.character(shortpangotab$Var1[i])
  }
  
}

kmertab<-as.data.frame(table(unlist(kmerlist)))
kmertab<-kmertab[which(kmertab$Freq==1),]
kmertab$Var1<-as.character(kmertab$Var1)
kmertab$Lineage<-NA
kmertab$Seq<-NA

pb<-txtProgressBar(max = nrow(kmertab))
for (i in 1:nrow(kmertab)) {
setTxtProgressBar(pb,i)
kmertab$Lineage[i]<- names(which(unlist(lapply(kmerlist, function(x) grep(paste("^",kmertab$Var1[i],"$",sep = ""), x)))>0))
}

shortpangotab$Covered2<-"NO"
shortpangotab$Covered2[which(shortpangotab$Var1 %in% unique(kmertab$Lineage))]<-"YES"



# Compression -------------------------------------------------------------
if(length(which(is.na(kmerlist)))>0) kmerlist<-kmerlist[-which(is.na(kmerlist))]

distance<-as.data.frame(matrix(NA, nrow = length(kmerlist), ncol = length(kmerlist)))
colnames(distance)<-names(kmerlist)
rownames(distance)<-names(kmerlist)

for (i in 1:length(kmerlist)) {
  for (j in 1:length(kmerlist)) {
    distance[i,j]<-length(which(kmerlist[[i]] %in% kmerlist[[j]]))/length(kmerlist[[i]])    
  }
}
  

uniquelin<-list()


for (i in 1:nrow(distance)) {
  
  lin2<-colnames(distance)[which(distance[i,]==1)]
  lin3<-strsplit(lin2,"\\.")
  headers<- unlist(lapply(lin3, function(x)x[1]))
  uniqueheads<- unique(headers)
  total.lin<-vector()
  for (j in 1:length(uniqueheads)) {
    lineages2<-lin2[grep(paste("^",uniqueheads[j],sep = ""), lin2)]
    lineages2<-strsplit(lineages2,"\\.")
    
    if(length(lineages2)>1){
      noe<-max(unlist(lapply(lineages2, length)))
      running<-TRUE
      while(running){
        new.lineages<- unlist(lapply(lineages2, function(x) paste(x[(1:min(length(x),noe))], collapse = ".") ))  
        if(length(unique(new.lineages))==1){
          new.lin<-unique(new.lineages)
          #unique.index<-c(unique.index,which(df$hash==uniqueregions[i])[1])
          running<-FALSE
          
        }else{
          new.lin<-NA
          if(noe==1) running<-FALSE
        }
        noe<-noe-1
      }
      
     new.lin<-paste(new.lin,".X",sep = "")
    }else{
     new.lin<-lin2[grep(paste("^",uniqueheads[j],sep = ""), lin2)]
    }
    total.lin<-c(total.lin,new.lin)
    
  }
  
  uniquelin<-c(uniquelin,list(total.lin))
  names(uniquelin)[length(uniquelin)]<-rownames(distance)[i]  
}


kmerlist2<-kmerlist
for (i in 1:length(kmerlist2)){
  #lins<-colnames(distance)[which(distance[i,]==1)]
  names(kmerlist2)[i]<-
    paste(unlist(uniquelin[which(names(uniquelin)==names(kmerlist2)[i])]),collapse = "|")
}
matchtable<-as.data.frame(names(kmerlist))
colnames(matchtable)<-"Old"
matchtable$New<-names(kmerlist2)

kmerlist2<-kmerlist2[-which(duplicated(names(kmerlist2)))]

kmertab2<-as.data.frame(table(unlist(kmerlist2)))
kmertab2<-kmertab2[which(kmertab2$Freq==1),]
kmertab2$Var1<-as.character(kmertab2$Var1)
kmertab2$Lineage<-NA
kmertab2$Seq<-NA

pb<-txtProgressBar(max = nrow(kmertab2))
for (i in 1:nrow(kmertab2)) {
  setTxtProgressBar(pb,i)
  kmertab2$Lineage[i]<- names(which(unlist(lapply(kmerlist2, function(x) grep(paste("^",kmertab2$Var1[i],"$",sep = ""), x)))>0))
}



for (i in 1:nrow(kmertab2)) {
  muts<-unlist(strsplit(kmertab2$Var1[i],"-"))
  sites<-as.numeric(gsub("^.","",gsub(".$","",muts)))
  subs<-as.character(gsub(".*[0-9]","",muts))
  subs<-subs[order(sites)]
  sites<-sites[order(sites)]
  mutref<-ref[[1]]
  
  kmerstoins<-unique(unlist(strsplit(unlist(kmerlist2[which(names(kmerlist2)==kmertab2$Lineage[i])]), "-") ))
  kmersites<-as.numeric(gsub("^.","",gsub(".$","",kmerstoins)))
  kmersubs<-as.character(gsub(".*[0-9]","",kmerstoins))
  
  mutref[kmersites]<-kmersubs
  
  if(sites[2]-sites[1]>3){
    mutref<-toupper(mutref[(sites[1]-1):(sites[2]+1)])
    mutref[4:(length(mutref)-3)]<-"."
  }else{
    mutref<-toupper(mutref[(sites[1]-2):(sites[2]+2)])
  }
  
  dels<-unique(unlist(strsplit(nc$deletions[which(nc$Nextclade_pango==kmertab2$Lineage[i])], ",")))
  
  if(length(dels)>0){
    del.sites<-vector()
    for (d in 1:length(dels)) {
      del.sites<-c(del.sites,c(as.numeric(gsub("-.*","",dels[d])):as.numeric(gsub(".*-","",dels[d]))))
    }
    if(length( which(c(sites[1]-1,sites[1]+1,sites[2]+1,sites[2]-1)%in%del.sites) )==0){
      kmertab2$Seq[i]<-paste(mutref,collapse = "")
    }else{
      kmertab2$Seq[i]<-"Deletion overlap"
    }
  }else{
    kmertab2$Seq[i]<-paste(mutref,collapse = "")
  }
  
}

write.csv(kmertab2, "/media/nacho/Data/DockerImages/HERCULES/CommonFiles/reference/kmerlist.csv")

