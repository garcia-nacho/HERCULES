library(ggplot2)
library(seqinr)
library(writexl)

depths<-list.files(pattern = "*._voc_depth.csv")
ref<-read.fasta("/Data/VariantSpike.fasta")
ref<-names(ref)
if(exists("out")) rm(out)
for (i in 1:length(depths)) {
  dummy<-read.csv(depths[i])
  dummy$Sample<-gsub("_voc_depth.csv","", depths[i])
  if(!exists("out")){
    out<-dummy

  }else{
    out<-rbind(out,dummy)
  }
}

if(length(which(out$Var1=="*"))>0){
  out.clean<-out[-which(out$Var1=="*"),]
}else{
  out.clean<-out
}

out.clean$Var1<-gsub("_Spike","",out.clean$Var1)

ggplot(out.clean)+
  geom_bar(aes(Sample, Freq), fill="red", stat="identity")+
  theme_minimal()+
  facet_wrap(~Var1)+
  ylab("Mapped Read Count")+
  xlab("Sample")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave("CountVariant_groupped_byVariant.pdf",width = 8.27,height = 11.69)

ggplot(out.clean)+
  geom_bar(aes(Var1, Freq), fill="red", stat="identity")+
  theme_minimal()+
  facet_wrap(~Sample)+
  ylab("Mapped Read Count")+
  xlab("Variant")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave("CountVariant_groupSample.pdf",width = 8.27,height = 11.69)

agg.df<-aggregate(Freq~Sample   , out.clean, sum)
colnames(agg.df)[2]<-"TotalCount"

out.clean<-merge(out.clean,agg.df, by="Sample",all.x=TRUE)
out.clean$Proportio<-out.clean$Freq/out.clean$TotalCount

ggplot(out.clean)+
  geom_bar(aes(Sample, Proportio), fill="red", stat="identity")+
  theme_minimal()+
  facet_wrap(~Var1)+
  ylab("Mapped Read Ratio")+
  xlab("Sample")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave("RatioVariant_groupVariant.pdf",width = 8.27,height = 11.69)

ggplot(out.clean)+
  geom_bar(aes(Var1, Proportio), fill="red", stat="identity")+
  theme_minimal()+
  facet_wrap(~Sample)+
  ylab("Mapped Read Count")+
  xlab("Variant")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave("RatioVariant_groupSample.pdf",width = 8.27,height = 11.69)


ggplot(out.clean)+
  geom_bar(aes(Sample, Proportio, fill=Var1), stat="identity")+
  theme_minimal()+
  scale_fill_manual(values = rainbow(length(unique(out.clean$Var1))))+
  labs(fill = "Variant")+
  ylab("Mapped Read Ratio")+
  xlab("Sample")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggsave("RatioVariantMapped_Stacked.pdf",width = 8.27,height = 11.69)

colnames(out)[1]<-"Variant"
colnames(out)[2]<-"ReadCount"

if(length(which(out$Variant=="*"))>0) out$Variant[which(out$Variant=="*")]<-"Unmapped"

write_xlsx(out, "VariantMappedReads.xlsx")
