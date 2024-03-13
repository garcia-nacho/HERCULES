library(writexl)
library(ggplot2)

#Plots for AF mash

files<-list.files( pattern = "_AF.lineages.csv")
if(exists("out.af")) rm(out.af)
for (i in 1:length(files)) {
  dummy<-read.csv(files[i],stringsAsFactors = FALSE)
  #dummy<-dummy[order(dummy$Freq, decreasing = TRUE)[1:10],]
  #dummy<-dummy[-which(dummy$out.vec=="Unclassified"),]
  dummy$Sample<-gsub("_AF.lineages.csv","",files[i])
  
  if(!exists("out.af")){
    out.af<-dummy
  }else{
    out.af<-rbind(out.af,dummy)
  }
}

pango.vars<-read.csv("https://raw.githubusercontent.com/cov-lineages/pango-designation/master/pango_designation/alias_key.json", sep = ":", header = FALSE)

if(nrow(pango.vars)>0){
  
  pango.vars<-pango.vars[-c(1,nrow(pango.vars)),]
  colnames(pango.vars)<-c("Short","Long")
  pango.vars$Long<-gsub(",","",pango.vars$Long)
  pango.vars$Long<-gsub(" ","",pango.vars$Long)
  pango.vars$Short<-gsub(" ","",pango.vars$Short)
  
  if(length(which(pango.vars$Long==""))>0){
    pango.vars$Long[which(pango.vars$Long=="")]<-pango.vars$Short[which(pango.vars$Long=="")]
  }
  
  short.lineage<-out.af$out.vec
  short.lineage<-gsub("\\..*","",short.lineage)
  coil<-gsub("[[:upper:]]","",out.af$out.vec)
  pango.replacement<-pango.vars$Short
  names(pango.replacement)<-pango.vars$Long
  start<-names(pango.replacement)[match(short.lineage, pango.replacement)]
  out.af$Pangolin_full<-paste(start, coil, sep = "")
  if(length(grep("NA", out.af$Pangolin_full))>0) out.af$Pangolin_full[grep("NA", out.af$Pangolin_full)]<-NA
  out.af$Pangolin_full<-gsub("\\.$",".X",out.af$Pangolin_full)
  
}

out.af$Lineage<-out.af$out.vec

write_xlsx(out.af[,-which(colnames(out.af)=="out.vec")], "AF_lineages.xlsx")

try(rm(out.af))
for (i in 1:length(files)) {
  dummy<-read.csv(files[i],stringsAsFactors = FALSE)
  d.total<-sum(dummy$Count)
  dummy<-dummy[order(dummy$Freq, decreasing = TRUE)[1:10],]
  #dummy<-dummy[-which(dummy$out.vec=="Unclassified"),]
  pad<-as.data.frame(t(c("Other", as.character(1-sum(dummy$Freq)), as.character(d.total-sum(dummy$Count)))))
  colnames(pad)<-colnames(dummy)
  dummy<-rbind(dummy,pad)
  dummy$Sample<-gsub("_AF.lineages.csv","",files[i])
  
  if(!exists("out.af")){
    out.af<-dummy
  }else{
    out.af<-rbind(out.af,dummy)
  }
}

out.af$Freq<-as.numeric(out.af$Freq)
ggplot(out.af)+
  geom_bar(aes(Sample, Freq,fill=out.vec),stat = "identity" )+
  theme_minimal()+
  scale_fill_manual(values=rainbow(length(unique(out.af$out.vec))))+
  xlab("Sample")+
  ylab("Frequency")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave("RatioVariantMash_Total_Stacked.pdf",width = 8.27,height = 11.69)

out.af$Count<-as.numeric(out.af$Count)
ggplot(out.af)+
  geom_bar(aes(Sample, Count,fill=out.vec),stat = "identity" )+
  theme_minimal()+
  scale_fill_manual(values=rainbow(length(unique(out.af$out.vec))))+
  xlab("Sample")+
  ylab("Read Count")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggsave("CountVariantMash_Total_Stacked.pdf",width = 8.27,height = 11.69)




try(rm(out.af))
for (i in 1:length(files)) {
  dummy<-read.csv(files[i],stringsAsFactors = FALSE)
  
  dummy<-dummy[-which(dummy$out.vec=="Unclassified"),]
  d.total<-sum(dummy$Count)
  
  dummy$Freq<-dummy$Freq/sum(dummy$Freq)
  dummy<-dummy[order(dummy$Freq, decreasing = TRUE)[1:10],]
  #dummy<-dummy[-which(dummy$out.vec=="Unclassified"),]
  pad<-as.data.frame(t(c("Other", as.character(1-sum(dummy$Freq)), as.character(d.total-sum(dummy$Count)))))
  colnames(pad)<-colnames(dummy)
  dummy<-rbind(dummy,pad)
  dummy$Sample<-gsub("_AF.lineages.csv","",files[i])
  
  if(!exists("out.af")){
    out.af<-dummy
  }else{
    out.af<-rbind(out.af,dummy)
  }
}

out.af$Freq<-as.numeric(out.af$Freq)
ggplot(out.af)+
  geom_bar(aes(Sample, Freq,fill=out.vec),stat = "identity" )+
  theme_minimal()+
  xlab("Sample")+
  ylab("Frequency")+
  scale_fill_manual(values=rainbow(length(unique(out.af$out.vec))))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave("RatioVariantMash_Clean_Stacked.pdf",width = 8.27,height = 11.69)


out.af$Count<-as.numeric(out.af$Count)
ggplot(out.af)+
  geom_bar(aes(Sample, Count,fill=out.vec),stat = "identity" )+
  theme_minimal()+
  xlab("Sample")+
  ylab("Read Count")+
  scale_fill_manual(values=rainbow(length(unique(out.af$out.vec))))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave("CountVariantMash_Clean_Stacked.pdf",width = 8.27,height = 11.69)

ggplot(out.af)+
  geom_bar(aes(out.vec, Freq),stat = "identity", fill="red" )+
  theme_minimal()+
  xlab("Variant")+
  ylab("Frequency")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  facet_wrap(~Sample)

 pdfh<- 11.69
 if(length(unique(out.af$out.vec))>10) pdfh<- round(pdfh+(length(unique(out.af$out.vec))/4),2) 
 pdfw<-round(0.70744*pdfh,2)
ggsave("RatioVariantMash_Clean_bySample.pdf",width = pdfh,height = pdfw)


ggplot(out.af)+
  geom_bar(aes(out.vec, Count),stat = "identity", fill="red" )+
  theme_minimal()+
  xlab("Variant")+
  ylab("Read Count")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  facet_wrap(~Sample)

pdfh<- 11.69
if(length(unique(out.af$out.vec))>10) pdfh<- round(pdfh+(length(unique(out.af$out.vec))/4),2) 
pdfw<-round(0.70744*pdfh,2)
ggsave("CountVariantMash_Clean_bySample.pdf",width = pdfh,height = pdfw)


ggplot(out.af)+
  geom_bar(aes(Sample, Freq),stat = "identity", fill="red" )+
  theme_minimal()+
  xlab("Sample")+
  ylab("Frequency")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  facet_wrap(~out.vec,ncol = 4)

ggsave("RatioVariantMash_Clean_byVariant.pdf",width = 8.27,height = 11.69)


ggplot(out.af)+
  geom_bar(aes(Sample, Count),stat = "identity", fill="red" )+
  theme_minimal()+
  xlab("Sample")+
  ylab("Read Count")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  facet_wrap(~out.vec,ncol = 4)

ggsave("CountVariantMash_Clean_byVariant.pdf",width = 8.27,height = 11.69)
