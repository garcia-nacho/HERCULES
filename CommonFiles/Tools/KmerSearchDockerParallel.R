library(data.table)
library(seqinr)
library("doParallel")
library("parallel")
library("foreach")
library(doSNOW)
library(progress)
library(ggplot2)
library(writexl)
library(plotly)
library(htmlwidgets)
arg=commandArgs(TRUE)

compressed<-list.files(full.names = TRUE, pattern = "b2f.tsv.gz",recursive = TRUE)
save.kmer<-FALSE
# ref<-read.fasta("/media/nacho/Data/DockerImages/Wastewater_SARS-CoV-2/CommonFiles/reference/SpikeRef.fa")
# kmers<-c("TAAGCATAGTG")
# kmers<-c("G23048A-C23202A", "C23453T-T23018G", "C23271T-C23423T-G23222A")


ref<-read.fasta("/home/docker/CommonFiles/reference/SpikeRef.fa")
ref<-unlist(ref)

kmers<-arg[1]
startn<-arg[2]
#startn<-1250
kmer<-unlist(strsplit(kmers,","))
kmername<-kmer
for (i in 1:length(kmer)) {
  
  if(length(grep("-",kmer[i]))>0){
    dummy<-unlist(strsplit(kmer[i],"-"))
    pos<-as.numeric(gsub("^.","",gsub(".$","",dummy)))
    if(min(pos)>3822) pos<-pos-21563+1
    subs<-gsub(".*[0-9]", "", dummy)
    temp.kmer<-ref
    temp.kmer[pos]<-tolower(subs)
    
#    if(length((min(pos)-1):(max(pos)+1))>7){
      
      to.unmask<-unique(c(pos,pos+1,pos-1))
      temp.kmer[-to.unmask]<-"."
      temp.kmer<-temp.kmer[c(startn:max(to.unmask))]
      temp.kmer<-c("^",temp.kmer)
      #temp.kmer<-temp.kmer[(min(pos)-1):(max(pos)+1)]
      #temp.kmer[4:(length(temp.kmer)-3)]<-"."
      #temp.kmer<-c("^",rep(".",  min(pos)-1 -startn ), temp.kmer)
    # }else{
    #   st<-max(pos)-min(pos)
    #   dif<-7-st
    #   temp.kmer<-temp.kmer[(min(pos)-dif):(max(pos)+dif)]
    #   temp.kmer<-c("^",rep(".",  min(pos)-1 -startn ), temp.kmer)
    # }
    kmer[i]<-toupper(paste(temp.kmer,collapse = ""))
    
  }else{
    kmerns<-strsplit(kmer[i],"")
    kmerns<-paste(unlist(lapply(kmerns, function(x) length(x)- length(which(toupper(x) %in% c("A","T","C","G") )) )),"N",sep = "")
    id<-unlist(strsplit(kmer[i],"\\."))
    
    if(length(which(id==""))>0) id<-id[-which(id=="")]
    if(length(id)==2){
      id<-gsub(" ","_",paste(id[1],kmerns,id[2]))
      if(kmerns[i]=="0N") id <-kmer[i]
    }
    
    if(!is.na(id)){
      kmername[i]<-id  
    }else{
      kmername[i]<-kmer[i]
    }
    
    
  }
  
}



dir.create("kmeruncompressed")

for (i in 1:length(compressed)) {
  target.file<-gsub(".*/","kmeruncompressed/",gsub(".gz","",compressed[i]))
  if(file.exists(target.file)) file.remove(target.file)
  system(paste("gunzip -c ", compressed[i], " > ",target.file,sep = ""))
}

files<-list.files("kmeruncompressed",full.names = TRUE )
kmer<-toupper(kmer)

for (k in 1:length(kmer)) {
  if(!dir.exists("ExtractedKmer")) dir.create("ExtractedKmer")
  df<-as.data.frame(files)
  print(paste("Searching for ", kmername[k],sep = ""))
 
  df$kmer<-kmer[k]
  df$kmername<-kmername[k]
  df$countQuery<-NA
  df$countRef<-NA
  df$RatioQuery<-NA
  df$RatioReference<-NA

  samples.to.analyze<-gsub(".*/","",gsub("\\.sorted.*","",files))
  
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
  
  cores<-as.numeric(detectCores())-2
  cores<-9
  cluster.cores<-makeCluster(cores)
  registerDoSNOW(cluster.cores)
  
  out.par<-foreach(i=1:nrow(df), .verbose=FALSE, .packages = c("data.table","seqinr"), .options.snow = opts) %dopar%{
     df2<-df
     gc()

     file.to.check<-fread(df$files[i])
     file.to.check<-as.data.frame(file.to.check)
     
     file.to.check$`#query_msa`<-toupper(file.to.check$`#query_msa`)
     file.to.check$ref_msa<-toupper(file.to.check$ref_msa)
     index.query<-grep(kmer[k],file.to.check$`#query_msa`)
     
     df2$countQuery[i]<-length(index.query)
     df2$countRef[i]<-length(grep(kmer[k],file.to.check$ref_msa))
     df2$RatioQuery[i]<-df2$countQuery[i]/nrow(file.to.check)
     df2$RatioReference[i]<-df2$countRef[i]/nrow(file.to.check)
     df2$TotalReadCount<-nrow(file.to.check)
     
     mutation.vector<-vector()
     if(length(index.query)>0){
       
     for (iq in index.query) {
       querydf<-unlist(strsplit(file.to.check$`#query_msa`[iq], ""))
       refdf<-unlist(strsplit(file.to.check$ref_msa[iq], "") )
       
       if(length(which(refdf=="-"))>0){
         querydf<-querydf[-which(refdf=="-")]
         refdf<-refdf[-which(refdf=="-")]
       }
       if(length(which(querydf=="-"))>0){
         delquery<-querydf
         delquery[which(querydf=="-")]<-refdf[which(querydf=="-")]
         
       }else{
         delquery<-querydf
       }
       
       mut.ref<-ref
       mut.ref[1250:2250]<-tolower(querydf)
       
       
       aa.mut<-translate(mut.ref)
       aa.ref<-translate(ref)
       mutation.vector<-c(mutation.vector,paste(aa.ref[which(aa.mut != aa.ref)],which(aa.mut != aa.ref), aa.mut[which(aa.mut != aa.ref)],sep = "") )
       
     }

     mutation.table<-as.data.frame(table(mutation.vector))
     colnames(mutation.table)<-c("Mutation","Count")
     mutation.table$Ratio<-mutation.table$Count/length(index.query)
     if(length(index.query) < 20){
       mutation.table<-mutation.table[which(mutation.table$Count>2),]
     }else{
       mutation.table<-mutation.table[which(mutation.table$Ratio>0.05),]  
     }
     if(nrow(mutation.table)>0){
     mutation.table$kmer<-kmer[k]
     mutation.table$kmer<-kmername[k]
     mutation.table$File<-gsub(".*/","",df$files[i])
     mutation.table$TotalCount<-nrow(file.to.check)
     }else{
       mutation.table<-as.data.frame(matrix(data = NA, nrow = 1, ncol = 5))
       colnames(mutation.table)<-c("Mutation","Count","Ratio", "kmer", "File")
       mutation.table$File<-gsub(".*/","",df$files[i])
       mutation.table$kmer<-kmer[k]
       mutation.table$kmer<-kmername[k]
       mutation.table$Count<-0
       mutation.table$Ratio<-0
       mutation.table$TotalCount<-nrow(file.to.check)
    }
     }else{
       mutation.table<-as.data.frame(matrix(data = NA, nrow = 1, ncol = 5))
       colnames(mutation.table)<-c("Mutation","Count","Ratio", "kmer", "File")
       mutation.table$File<-gsub(".*/","",df$files[i])
       mutation.table$kmer<-kmer[k]
       mutation.table$kmer<-kmername[k]
       mutation.table$Count<-0
       mutation.table$Ratio<-0
       mutation.table$TotalCount<-nrow(file.to.check)
     }
     
     if(length(index.query)>0 & save.kmer ){ 
       write.table(file.to.check[index.query,], 
                                            paste(gsub(".tsv",paste("_",kmername[k],".tsv",sep=""),gsub(".*/","ExtractedKmer/",df2$files[i])),sep = ""),sep = "\t", quote = FALSE, row.names = FALSE)
        system(paste("gzip ",gsub(".tsv",paste("_",kmername[k],".tsv",sep=""),gsub(".*/","ExtractedKmer/",df2$files[i])),sep = ""))
       }
     rm(file.to.check)
     return(list(df2[i,], mutation.table))
     
  }
  stopCluster(cluster.cores)
  outpar1<-lapply(out.par, function(x)x[[1]])
  outpar2<-lapply(out.par, function(x)x[[2]])
  dfunlist<-do.call(rbind, outpar1)
  mutlist<-do.call(rbind, outpar2)
  
  dfunlist$Sample<-gsub(".*/","",dfunlist$files)
  dfunlist$Sample<-gsub(".sorted.*","",dfunlist$Sample)
  
  mutlist$Sample<-gsub(".*/","",mutlist$File)
  mutlist$Sample<-gsub(".sorted.*","",mutlist$Sample)
  
  
  if(!exists("df.out")){
    df.out<-dfunlist
  }else{
    df.out<-rbind(df.out, dfunlist)
  }
  
  if(!exists("mut.out")){
    mut.out<-mutlist
  }else{
    mut.out<-rbind(mut.out, mutlist)
  }

  
  
}

file.remove(list.files("kmeruncompressed/",full.names = TRUE) )


write_xlsx(df.out, "KmerSearchResults.xlsx")
write_xlsx(mut.out, "KmerSearchMutations.xlsx")

df.out$Location<-gsub(".*_","",df.out$Sample)
df.out$Date<-as.Date(gsub(".*\\.","",gsub("_.*","",df.out$Sample)))

if(length(which(is.na(df.out$RatioQuery)))>0) df.out<-df.out[-which(is.na(df.out$RatioQuery)),]

df.agg.sum<-aggregate(countQuery ~ Date + Location + kmername, df.out,sum )
df.agg.sum$dateloc<-paste(df.agg.sum$Date, df.agg.sum$Location, df.agg.sum$kmername, sep = "_")
df.agg.total<-aggregate(TotalReadCount ~ Date + Location + kmername, df.out,sum )
df.agg.total$dateloc<-paste(df.agg.total$Date, df.agg.total$Location, df.agg.total$kmername, sep = "_")


df.agg.sum<-merge(df.agg.sum, df.agg.total[,c("TotalReadCount", "dateloc")], by="dateloc")
df.agg.sum$RatioQuery<-df.agg.sum$countQuery/df.agg.sum$TotalReadCount


ggplot(df.agg.sum)+
  geom_line(aes(Date,log(RatioQuery) ))+
  theme_minimal()+
  facet_wrap(~kmername+Location)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ylab("Ratio")

ggsave("KmerSearchPlotLog.pdf" ,width = 20, height = 16)

ggplot(df.agg.sum)+
  geom_line(aes(Date,RatioQuery ))+
  theme_minimal()+
  facet_wrap(~kmername+Location)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ylab("Ratio")

ggsave("KmerSearchPlot.pdf" ,width = 20, height = 16)

# ggplot(df.agg.sum)+
#   geom_line(aes(Date,-1/log10(RatioQuery) ))+
#   theme_minimal()+
#   facet_wrap(~kmer+Location)+
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
#   ylab("Ratio")

mut.out.bk<-mut.out
kmers<-unique(mut.out.bk$kmer)

for (i in 1:length(kmers)) {
  mut.out<-mut.out.bk[which(mut.out.bk$kmer==kmers[i]),]
  if(length(grep("X",mut.out$Mutation))>0) mut.out<-mut.out[-grep("X",mut.out$Mutation),]
  if(length(grep("\\*$",mut.out$Mutation))>0)mut.out<-mut.out[-grep("\\*$",mut.out$Mutation),]
  
  mut.out$Date<-as.Date(gsub(".*\\.","",gsub("_.*","",mut.out$Sample)))
  mut.out$Detected<-"YES"
  
  if(length(which(is.na(mut.out$Mutation)))>0) mut.out<-mut.out[-which(is.na(mut.out$Mutation)),]
  if(nrow(mut.out)>0){
  pad<-expand.grid(unique(mut.out$Date), unique(mut.out$Mutation))
  colnames(pad)<-c("Date", "Mutation")
  pad$Date<-as.character(pad$Date)
  pad$Mutation<-as.character(pad$Mutation)

  pad<-pad[-which(paste(pad$Date, pad$Mutation) %in% paste(mut.out$Date, mut.out$Mutation)),]
  pad$Ratio<-0

  toplot<-rbind(mut.out[,c("Date","Mutation","Ratio")], pad)
  toplot$Detected<-"NO"
  toplot$Detected[which(toplot$Ratio>0)]<-"YES"

  toplot$Date<-as.character(toplot$Date)
  toplot$Mutation<-paste("S:",toplot$Mutation,sep = "")
  toplot$Mutation<-factor(toplot$Mutation, levels = unique(toplot$Mutation)[order(as.numeric(gsub("S:.","",gsub(".$","",unique(toplot$Mutation)))))])
  ggplot(toplot)+
  geom_tile(aes(Mutation, Date, fill=Detected,alpha=Ratio))+
  #scale_fill_gradient(low = "blue", high = "red")+
  scale_fill_manual(values =c("white","red"))+
  theme_minimal()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ggtitle(paste("Kmer:",kmers[i]))

ggsave(paste("KmerMutations",kmers[i],".pdf",sep = "") ,width = 40, height = 20)

}
}

