
path=commandArgs(TRUE)
co.n<-0.15
pathf<-path[1]
co.n<-as.numeric(path[2])
n.start<-as.numeric(path[3])
n.end<- as.numeric(path[4])

noise.path<-list.files(pathf, pattern = ".*noise.tsv",full.names = TRUE, recursive = TRUE)

poi<-vector()
for (i in 1:length(noise.path)) {
  noise.table.d<-read.csv(noise.path[i], sep = "\t", header = FALSE) 
  poi<-c(poi, noise.table.d$V1[which(noise.table.d$V2>=co.n)])
}
poi<-unique(poi)

write.csv(poi, "poi.temp",row.names = FALSE)