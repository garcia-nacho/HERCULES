library(readxl)
library(ggsankey)
library(tidyverse)
library(writexl)

r.path<-"/Data/results/"

include.deletions<-FALSE

results<-list.files(r.path, full.names = TRUE)

excels<-results[grep("xlsx", results)]
excels<-excels[grep("results.xlsx", excels)]
dummync<-read.csv(paste(r.path,"Nextclade.results.csv",sep = ""), sep = ";")
#dummync<-read.csv(paste(r.path,"Variants.nextclade.csv",sep = ""), sep = ";")
colnames(dummync)[1]<-"variantID"
dummync$variantID<-gsub("ID-","", dummync$variantID)
if(length(which(duplicated(dummync)))>0 ) dummync<-dummync[-which(duplicated(dummync)),]

if(exists("results.out"))try(rm(results.out))

if(length(grep("dependent",excels))>0){
  excels<-excels[grep("dependent",excels)]
}

for (i in 1:length(excels)) {
  dummy<-read_xlsx(excels[i])
  if(nrow(dummy)>0){
  dummy<- merge(dummy, dummync[,c("variantID", "clade","Nextclade_pango", "aaSubstitutions","substitutions")], by ="variantID", all.x = TRUE)
  samplename<-gsub("\\.results.*","",gsub(".*/","",excels[i]))
  dummy$Sample<-samplename
  dummy$ImputedDeletions<-NA
  for (j in 1:nrow(dummy)) {
    deletions<-colnames(dummy)[which(dummy[j,]=="D")]
    if(length(deletions)>0){
      deletions<-as.numeric(gsub("P_","",deletions))
      deletions<-deletions/3
      if( length(which(deletions %% 3 == 1))>0) deletions[which(deletions %% 3 == 1)]<-deletions[which(deletions %% 3 == 1)]+1
      deletions<-round(deletions)
      dummy$ImputedDeletions[j]<-paste(paste("S:",deletions,"-",sep = ""), collapse = ",")
      rm(deletions)
    }
  }
  if(!exists("results.out")){
    results.out<-dummy[,c( "variantID" ,  "count","ratio",  "clade","Nextclade_pango","aaSubstitutions","Sample","ImputedDeletions","substitutions")]
  }else{
    results.out<-rbind(results.out, dummy[,c( "variantID" ,  "count","ratio",  "clade","Nextclade_pango","aaSubstitutions","Sample","ImputedDeletions","substitutions")])
  }
  }
}



#results.v<-list.files(paste(r.path,"variants",sep = ""), full.names = TRUE)
#variantfiles<-results.v[grep("variants.csv", results.v)]
results.out$Mutations<-results.out$aaSubstitutions
if(length(which(results.out$Mutations==""))>0) results.out$Mutations[which(results.out$Mutations=="")]<-"None"



results.to.plot<-results.out[which(results.out$ratio>0.015),c("Sample", "Mutations", "ratio")]
results.to.plot<-aggregate(ratio~Mutations + Sample,  results.to.plot, sum)

samples.toplot<-unique(results.to.plot$Sample)
if(exists("dummy.out"))try(rm(dummy.out))
for (i in 1:length(samples.toplot)) {
  dummy<-results.to.plot[1,,drop=FALSE]
  dummy$Sample<-samples.toplot[i]
  dummy$Mutations<-"Other Combination"
  dummy$ratio<-1-sum(results.to.plot$ratio[which(results.to.plot$Sample==samples.toplot[i])])
  if(!exists("dummy.out")){
    dummy.out<-dummy
  }else{
    dummy.out<-rbind(dummy.out,dummy)
  }
}

numberofreads<-aggregate(count~Sample, results.out, sum )
results.to.plot<-rbind(results.to.plot,dummy.out)

ggplot(results.to.plot)+
  geom_bar(aes(Sample, ratio, fill=Mutations),position = "stack", stat = "identity")+
  scale_fill_manual(values =  c("grey",rainbow(length(unique(results.to.plot$Mutations))-1)))+
  theme(legend.position="bottom")+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ylab("Ratio")
  
ggsave(paste(r.path, "MutationsBarplot.pdf",sep=""), width =   25, height = 12)
agg.mut<-results.to.plot

agg.mut<-merge(results.to.plot, numberofreads, by="Sample")
agg.mut$count<-agg.mut$count*agg.mut$ratio

 resultsagg<-aggregate(ratio~clade+Sample, results.out, sum)
 if(length(which(resultsagg$clade==""))>0) resultsagg$clade[which(resultsagg$clade=="")]<-"NA"
# ggplot(resultsagg)+
#   geom_bar(aes(Sample, ratio, fill=clade),position = "stack", stat = "identity")+
#   scale_fill_manual(values =  rainbow(length(unique(resultsagg$clade))))+
#   theme_minimal()+
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
# ggsave(paste(r.path, "Clade_Barplot.pdf",sep=""), width =   12, height = 12)
 agg.clades<-resultsagg
 agg.clades<-merge(agg.clades, numberofreads, by="Sample")
 agg.clades$count<-agg.clades$count*agg.clades$ratio

 resultsagg<-aggregate(ratio~Nextclade_pango+Sample, results.out, sum)
 if(length(which(resultsagg$Nextclade_pango==""))>0) resultsagg$Nextclade_pango[which(resultsagg$Nextclade_pango=="")]<-"NA"
# ggplot(resultsagg)+
#   geom_bar(aes(Sample, ratio, fill=Nextclade_pango),position = "stack", stat = "identity")+
#   scale_fill_manual(values =  rainbow(length(unique(resultsagg$Nextclade_pango))))+
#   theme_minimal()+
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
# ggsave(paste(r.path, "Pangolin_Barplot.pdf",sep=""), width =   12, height = 12)

resultsagg<-merge(resultsagg, numberofreads, by="Sample")
resultsagg$count<-resultsagg$count*resultsagg$ratio
write_xlsx(list("Mutations"=agg.mut, "Clades"=agg.clades, "NCPangolin"=resultsagg ),paste(r.path, "Results.Aggregated.xlsx",sep=""))

samples<-unique(results.out$Sample)
sk.list<-list()


for(s.method in c("Unique","Total")){
  for (sk in 1:length(samples)) {
    
    df<-results.out[which(results.out$Sample==samples[sk]),]
    
    df<-df[order(df$count, decreasing = TRUE),]
    try(rm(out))
    for (i in 1:min(nrow(df),20)) {
      try(rm(mutpos.df))
      if(df$aaSubstitutions[i]!=""){
        mutations<-unlist(base::strsplit(df$aaSubstitutions[i], ","))
        mutpos<- gsub("S:.","",gsub(".$","",mutations))
        
        mutpos.df<-as.data.frame(t(mutations))
        colnames(mutpos.df)<-mutpos
      }
      try(rm(delpos))
      if(!is.na(df$ImputedDeletions[i])){
        deletions<-unlist(base::strsplit(df$ImputedDeletions[i], ","))
        delpos<- gsub("S:","",gsub("-$","",deletions))
        delpos.df<-as.data.frame(t(deletions))
        colnames(delpos.df)<-delpos
        if(exists("mutpos.df")) mutpos.df<-cbind(mutpos.df, delpos.df)
        if(!exists("mutpos.df")) mutpos.df<-delpos.df
      }
      
      if(df$aaSubstitutions[i]=="" & is.na(df$ImputedDeletions[i])){
        mutpos.df<-as.data.frame(matrix("Ref", nrow = 1, ncol = 1))
        colnames(mutpos.df)<-"Dummy"
      }
      if(length(grep("-",mutpos.df[1,]))>0)mutpos.df<-mutpos.df[,-grep("-", mutpos.df[1,]),drop=FALSE]
      
      mutpos.df<-mutpos.df[round(rep(1,df$count[i])),,drop(FALSE)]
      
      
      
      
      
      if(!exists("out")){
        out<-mutpos.df
      }else{
        
        if(length(which(colnames(mutpos.df) %in% colnames(out)))>0){
          missing.out<- colnames(mutpos.df)[-which(colnames(mutpos.df) %in% colnames(out))]
        }else{
          missing.out<- colnames(mutpos.df)
        }
        
        if(length(which(colnames(out) %in% colnames(mutpos.df)))>0){
          missing.df<- colnames(out)[-which(colnames(out) %in% colnames(mutpos.df))]
        }else{
          missing.df<- colnames(out)
        }
        
        
        # 
        # if(length(missing.out)==0 & length(missing.df)==0){
        #   pad.out<-as.data.frame(matrix(data = "Ref", nrow = nrow(out), ncol = ncol(mutpos.df)))
        #   colnames(pad.out)<-colnames(mutpos.df)
        #   outnew<-cbind(out, pad.out)
        #   
        #   pad.df<-as.data.frame(matrix(data = "Ref", nrow = nrow(mutpos.df), ncol = ncol(out))) 
        #   colnames(pad.df)<-colnames(out)
        #   mutpos.df<- cbind(mutpos.df, pad.df)
        #   out<-outnew
        # }
        # 
        if(length(missing.out)>0){
          pad.out<-as.data.frame(matrix(data = "Ref", nrow = nrow(out), ncol = length(missing.out)))
          colnames(pad.out)<-missing.out
          out<-cbind(out, pad.out)
        }
        
        if(length(missing.df)>0){
          pad.df<-as.data.frame(matrix(data = "Ref", nrow = nrow(mutpos.df), ncol = length(missing.df))) 
          colnames(pad.df)<-missing.df
          mutpos.df<- cbind(mutpos.df, pad.df)
        }
        
        out<-rbind(out, mutpos.df)  
      }
      
    }
    
    try(out$Dummy<-NULL)
    if(s.method=="Unique"){
      c.to.remove<- apply(out,2,function(x)length(unique(x)))
      out<-out[,-which(c.to.remove==1), drop=FALSE]
    }
    
    if(ncol(out)>0){
      df.sankey <- make_long(out, colnames(out))
      
      ggplot(df.sankey, aes(x = x, 
                            next_x = next_x, 
                            node = node, 
                            next_node = next_node,
                            fill = factor(node))) +
        geom_sankey(flow.alpha = 0.5, node.color = 1) +
        scale_fill_viridis_d(option = "A", alpha = 0.5) +
        theme_sankey(base_size = 16)+
        xlab("Position")+
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
        ggtitle(paste("Mutation combinations in", samples[sk]))+
        labs(fill = "Mutations")
      ggsave(paste(r.path, samples[sk],"_Sankeyplot",s.method,".AA.pdf",sep = ""), width =   20, height = 10)
      
      out.clean<-out[,-which(apply(out, 2, function(x) length(grep("-",x)))>0),drop=FALSE] 
      
      if(ncol(out.clean)>0){
        df.sankey2 <- make_long(out.clean, colnames(out.clean))
        
        ggplot(df.sankey2, aes(x = x, 
                               next_x = next_x, 
                               node = node, 
                               next_node = next_node,
                               fill = factor(node))) +
          geom_sankey(flow.alpha = 0.5, node.color = 1) +
          scale_fill_viridis_d(option = "A", alpha = 0.5) +
          theme_sankey(base_size = 16)+
          xlab("Position")+
          theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
          ggtitle(paste("Mutation combinations in", samples[sk]))+
          labs(fill = "Mutations")
        ggsave(paste(r.path, samples[sk],"_SankeyplotNoDel",s.method,".AA.pdf",sep = ""), width =   20, height = 10)
        
      }
    }
  }
}




samples<-unique(results.out$Sample)
sk.list<-list()

for (sk in 1:length(samples)) {
  
  df<-results.out[which(results.out$Sample==samples[sk]),]
  
  df<-df[order(df$count, decreasing = TRUE),]
  if(exists("out"))try(rm(out))
  for (i in 1:min(nrow(df),20)) {
    try(rm(mutpos.df))
    if(df$substitutions[i]!=""){
      mutations<-unlist(base::strsplit(df$substitutions[i], ","))
      mutpos<- gsub("^.","",gsub(".$","",mutations))
      
      mutpos.df<-as.data.frame(t(mutations))
      colnames(mutpos.df)<-mutpos
    }
    try(rm(delpos))
    
    
    if(df$aaSubstitutions[i]=="" & is.na(df$ImputedDeletions[i])){
      mutpos.df<-as.data.frame(matrix("Ref", nrow = 1, ncol = 1))
      colnames(mutpos.df)<-"Dummy"
    }
    
    if(exists("mutpos.df")){
    
    if(length(grep("-",mutpos.df[1,]))>0) mutpos.df<-mutpos.df[,-grep("-", mutpos.df[1,]),drop=FALSE]
    
    mutpos.df<-mutpos.df[round(rep(1,df$count[i])),,drop(FALSE)]
    
    
    
    
    
    if(!exists("out")){
      out<-mutpos.df
    }else{
      
      if(length(which(colnames(mutpos.df) %in% colnames(out)))>0){
        missing.out<- colnames(mutpos.df)[-which(colnames(mutpos.df) %in% colnames(out))]
      }else{
        missing.out<- colnames(mutpos.df)
      }
      
      if(length(which(colnames(out) %in% colnames(mutpos.df)))>0){
        missing.df<- colnames(out)[-which(colnames(out) %in% colnames(mutpos.df))]
      }else{
        missing.df<- colnames(out)
      }
      
      
      # 
      # if(length(missing.out)==0 & length(missing.df)==0){
      #   pad.out<-as.data.frame(matrix(data = "Ref", nrow = nrow(out), ncol = ncol(mutpos.df)))
      #   colnames(pad.out)<-colnames(mutpos.df)
      #   outnew<-cbind(out, pad.out)
      #   
      #   pad.df<-as.data.frame(matrix(data = "Ref", nrow = nrow(mutpos.df), ncol = ncol(out))) 
      #   colnames(pad.df)<-colnames(out)
      #   mutpos.df<- cbind(mutpos.df, pad.df)
      #   out<-outnew
      # }
      # 
      if(length(missing.out)>0){
        pad.out<-as.data.frame(matrix(data = "Ref", nrow = nrow(out), ncol = length(missing.out)))
        colnames(pad.out)<-missing.out
        out<-cbind(out, pad.out)
      }
      
      if(length(missing.df)>0){
        pad.df<-as.data.frame(matrix(data = "Ref", nrow = nrow(mutpos.df), ncol = length(missing.df))) 
        colnames(pad.df)<-missing.df
        mutpos.df<- cbind(mutpos.df, pad.df)
      }
      
      out<-rbind(out, mutpos.df)  
    }
    }
  }
  
  try(out$Dummy<-NULL)
  
  out<-out[,-which(apply(out,2, function(x) length(unique(x))) ==1), drop=FALSE]
  
  if(ncol(out)>0){
  df.sankey <- make_long(out, colnames(out))
  
  ggplot(df.sankey, aes(x = x, 
                        next_x = next_x, 
                        node = node, 
                        next_node = next_node,
                        fill = factor(node))) +
    geom_sankey(flow.alpha = 0.5, node.color = 1) +
    scale_fill_viridis_d(option = "A", alpha = 0.5) +
    theme_sankey(base_size = 16)+
    xlab("Position")+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
    ggtitle(paste("Mutation combinations in", samples[sk]))+
    labs(fill = "Mutations")
  ggsave(paste(r.path, samples[sk],"_Sankeyplot.Nucleotide.pdf"), width =   40, height = 20)
  
  out.clean<-out[,-which(apply(out, 2, function(x) length(grep("-",x)))>0),drop=FALSE] 
  }
}


