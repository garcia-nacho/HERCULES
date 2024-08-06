library(ggplot2)
path=commandArgs(TRUE)
path<-path[1]
co.n<-as.numeric(path[2])


noise.files<-list.files(path, full.names = TRUE, pattern = ".noise.tsv",recursive = TRUE)

for (i in 1:length(noise.files)) {
  dum<-read.csv(noise.files[i],sep = "\t", header = FALSE)
  dum$id<-gsub("noise.tsv","",gsub(".*/","",noise.files[i]))
  if(!exists("noise.table")){
    noise.table<-dum
  }else{
    noise.table<-rbind(noise.table, dum)
  }
}


ggplot(noise.table)+
  geom_line(aes(V1,V2))+
  xlab("Position")+
  ylab("Noise")+
  geom_hline(yintercept=co.n, linetype='dotted', col = 'red')+
  theme_minimal()+
  facet_wrap(~id,ncol = 4)+
  ggtitle("Noise")

ggsave(paste(path, "NoiseMosaic.pdf",sep = ""))

ggplot(noise.table)+
  geom_line(aes(V1,V3))+
  xlab("Position")+
  ylab("Depth")+
  facet_wrap(~id,ncol = 4)+
  theme_minimal()+
  ggtitle("Coverage")
ggsave(paste(path, "CoverageMosaic.pdf",sep = ""))