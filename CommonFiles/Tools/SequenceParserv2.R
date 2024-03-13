library(seqinr)
library(digest)

#newfiles

df.new<-read.csv("/media/nacho/Data/wastewater/ReferencesProb/update28082023/update28082023.csv", sep = ";")
df<-read.csv("/media/nacho/Data/wastewater/AFMode/refs/references.csv",sep = ";")

df<-df[,-which(colnames(df) %in% c("clade_legacy", "clade_display"))]
df.new<-df.new[,-which(colnames(df.new) %in% c("clade_legacy", "clade_display"))]

df<-rbind(df, df.new)

sqs<-read.fasta("/media/nacho/Data/wastewater/AFMode/refs/references.aligned.fasta")
sqs.new<-read.fasta("/media/nacho/Data/wastewater/ReferencesProb/update28082023/aligned_update29082023.fasta")

sqs<-c(sqs,sqs.new)

df<-df[-which(duplicated(df$seqName)),]
sqs<-sqs[-which(duplicated(names(sqs)))]

df<-df[which(df$totalMissing<300),]
sqs<-sqs[which(names(sqs) %in% gsub(" .*","",df$seqName))]

sqs<-lapply(sqs, function(x) x[c((21563+1250):(21563+2250))] )

hashes<-lapply(sqs, function(x) digest(paste(x,collapse = ""),algo="md5"))

df$hash<-unlist(hashes)

uniqueregions<-unique(hashes)
df$QC<-NA
df$Corrected<-NA
unique.index<-vector()


for (i in 1:length(uniqueregions)) {
  pangos<- df$Nextclade_pango[ which(df$hash== uniqueregions[i])]
  pangos.df<-as.data.frame(table(pangos))
  pangos.df$Freq<-pangos.df$Freq/sum(pangos.df$Freq)
  pangos.df$Hash<-paste(strsplit(uniqueregions[[i]],"")[[1]][1:10], collapse = "")
  
  
  if(length(unique(pangos))==1){
    df$QC[which(df$hash==uniqueregions[i])]<-"OK"
    df$Corrected[which(df$hash==uniqueregions[i])]<-df$Nextclade_pango[which(df$hash==uniqueregions[i])]
    #unique.index<-c(unique.index,which(df$hash==uniqueregions[i])[1])
    
    
  }else{
    df$QC[which(df$hash==uniqueregions[i])]<-"Warning"
    
    lineages<- df$Nextclade_pango[which(df$hash==uniqueregions[i])]
    lineages<-strsplit(lineages,"\\.")
    noe<-max(unlist(lapply(lineages, length)))
    
    running<-TRUE
    while(running){
      new.lineages<- unlist(lapply(lineages, function(x) paste(x[(1:min(length(x),noe))], collapse = ".") ))  
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
    if(is.na(new.lin)){
      df$Corrected[which(df$hash==uniqueregions[i])]<-paste(new.lin,".",paste(strsplit(uniqueregions[[i]],"")[[1]][1:10], collapse = "") ,".X",sep = "")
    }else{
      df$Corrected[which(df$hash==uniqueregions[i])]<-paste(new.lin,".X",sep = "")
    }
    
    
  }
  
  pangos.df$Lineage<-df$Corrected[which(df$hash==uniqueregions[i])][1]
  if(!(exists("pangos.out"))){
    pangos.out<-pangos.df
  }else{
    pangos.out<-rbind(pangos.out, pangos.df)
  }
  
}

df<-df[-which(duplicated(df$hash)),]

sqs<-sqs[which(names(sqs) %in% gsub(" .*","",df$seqName))]

lineages.to.parse<-unique(df$Corrected)

df$seqName<-gsub(" .*","",df$seqName)
for (i in 1:length(lineages.to.parse)) {
  names(sqs)[which(names(sqs) %in% df$seqName[which(df$Corrected==lineages.to.parse[i])])]<-
    paste(lineages.to.parse[i],"_SQ",c(1: length(df$seqName[which(df$Corrected==lineages.to.parse[i])])),sep = "")
}

write.fasta(sqs, "/media/nacho/Data/wastewater/ReferencesProb/ReferencesAligned.28082023.fasta",names = names(sqs))
write.csv(pangos.out, "/media/nacho/Data/wastewater/ReferencesProb/variant_hash.csv", row.names = FALSE)

spikes<-lapply(sqs,function(x) paste(toupper(x),collapse = ""))

spikes_tsv<-as.data.frame(unlist(spikes))

colnames(spikes_tsv)<-"#query_msa"
spikes_tsv$`#query_msa`<-paste("A",spikes_tsv$`#query_msa`,sep="")
spikes_tsv$ref_msa<-"agattgctgattataattataaattaccagatgattttacaggctgcgttatagcttggaattctaacaatcttgattctaaggttggtggtaattataattacctgtatagattgtttaggaagtctaatctcaaaccttttgagagagatatttcaactgaaatctatcaggccggtagcacaccttgtaatggtgttgaaggttttaattgttactttcctttacaatcatatggtttccaacccactaatggtgttggttaccaaccatacagagtagtagtactttcttttgaacttctacatgcaccagcaactgtttgtggacctaaaaagtctactaatttggttaaaaacaaatgtgtcaatttcaacttcaatggtttaacaggcacaggtgttcttactgagtctaacaaaaagtttctgcctttccaacaatttggcagagacattgctgacactactgatgctgtccgtgatccacagacacttgagattcttgacattacaccatgttcttttggtggtgtcagtgttataacaccaggaacaaatacttctaaccaggttgctgttctttatcaggatgttaactgcacagaagtccctgttgctattcatgcagatcaacttactcctacttggcgtgtttattctacaggttctaatgtttttcaaacacgtgcaggctgtttaataggggctgaacatgtcaacaactcatatgagtgtgacatacccattggtgcaggtatatgcgctagttatcagactcagactaattctcctcggcgggcacgtagtgtagctagtcaatccatcattgcctacactatgtcacttggtgcagaaaattcagttgcttactctaataactctattgccatacccacaaattttactattagtgttaccacagaaattctaccagtgtctatgaccaagacatcagtagattgtacaatgtacatttgtggtgattcaactgaatgcagc"

write.table(spikes_tsv, "/media/nacho/Data/wastewater/ReferencesProb/MSA_Refs.tsv", row.names = FALSE, quote = FALSE, sep = "\t")
write.csv(rownames(spikes_tsv),"/media/nacho/Data/wastewater/ReferencesProb/MSA_RefsID.csv",row.names = FALSE)

