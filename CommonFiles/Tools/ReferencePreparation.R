library(seqinr)

#Reference preparation after nextclade

#total<-read.fasta("/media/nacho/Data/wastewater/ReferencesProb/update28082023/update28082023.fasta")
alig<-read.fasta("/media/nacho/Data/wastewater/ReferencesProb/update28082023/aligned_update29082023.fasta")
metadata<-read.csv("/media/nacho/Data/wastewater/ReferencesProb/update28082023/update28082023.csv",sep = ";")

metadata$seqName<-gsub(" .*", "", metadata$seqName)

metadata<-metadata[-which(duplicated(metadata$seqName)),]
alig<-alig[-which(duplicated(names(alig)))]

start.number<-5282

metadata$SID<-paste("SID", c((1:nrow(metadata))+5282),"_",metadata$Nextclade_pango,sep = "" )

pb<-txtProgressBar(max = length(alig))

for (i in 1:length(alig)) {
  setTxtProgressBar(pb,i)
  names(alig)[i]<-metadata$SID[which(metadata$seqName == names(alig)[i])]
  alig[[i]]<-alig[[i]][21563:25384]
}

refs<-read.fasta("/media/nacho/Data/wastewater/ReferencesProb/Spike17022022.fasta")

refs<-c(refs,alig)
write.fasta(refs, file.out = "/media/nacho/Data/wastewater/ReferencesProb/Spike28082023.fasta", names = names(refs))

