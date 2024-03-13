library(seqinr)
sums<-list.files("/media/nacho/Data/wastewater/refs/", pattern = "report_v13", recursive = TRUE, full.names = TRUE)
fastas<-list.files("/media/nacho/Data/wastewater/refs/", pattern = "_ivar.consensus.masked.fa", recursive = TRUE, full.names = TRUE)

if(exists("out")) rm(out)
for (i in 1:length(sums)) {
  dummy<-read.csv(sums[i], sep = "\t")
  if(!exists("out")){
    out<-dummy
  }else{
    out<-rbind(out,dummy)
  }
}

out.clean<-out[which(out$WGS_pct20x>99),]

table(out$pangolin_ivar_lineage)

out<-out[which(out$pangolin_ivar_lineage %in% c("BA.2", "BA.4", "BA.5", "BE.1", "BQ.1.1","BR.1")),]
out<-out[-grep("TMP",out$Name),]
out<-out[-grep("hS",out$Name),]

lineages<-unique(out$pangolin_ivar_lineage)

fasta.to.get<-list()
for (i in 1:length(lineages)) {
    id<- out$Name[which(out$pangolin_ivar_lineage==lineages[i] & out$WGS_mean == max(  out$WGS_mean[which(out$pangolin_ivar_lineage==lineages[i])]) )][1]
    fasta.to.get[[i]]<- read.fasta(fastas[ grep(id,fastas) ])[[1]]
    names(fasta.to.get)[i]<-lineages[i]
}

write.fasta(fasta.to.get,paste("/media/nacho/Data/wastewater/refs/variantRefs", gsub("-","",Sys.Date()), ".fasta", sep = ""), names = names(fasta.to.get))