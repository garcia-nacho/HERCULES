

if(length(grep(grep("N",lineages.df.ref$Seq)))>0)lineages.df.ref<-lineages.df.ref[-grep("N",lineages.df.ref$Seq),]
if(length(which(is.na(lineages.df.ref$Mut.aa)))>0) lineages.df.ref$Mut.aa[which(is.na(lineages.df.ref$Mut.aa))]<-"None"

variant.hash<-read.csv(hasshes)

if(length(grep("NA\\..*\\.X", lineages.df.ref$ID2))>0){
  to.match<-unique(lineages.df.ref$ID2[grep("NA\\..*\\.X", lineages.df.ref$ID2)])  
  
  for(tm in c(1:length(to.match))){
    variant.temp<-variant.hash[which(variant.hash$Lineage==to.match[tm]),]
    variant.temp<-variant.temp[order(variant.temp$Freq, decreasing = TRUE),]
    variant.temp<-variant.temp[c(1:min(3,nrow(variant.temp))),]
    
    lineages2<- variant.temp$pangos
    lineages2<-strsplit(lineages2,"\\.")
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
    if(is.na(new.lin)){
      newid<-paste(variant.temp$pangos,collapse =  "/")
    }else{
      newid<-paste(new.lin,".X",sep = "")
    }
  
    lineages.df.ref$ID2[which(lineages.df.ref$ID2==to.match[tm])]<-newid
  }
}

to.col.ref<-as.data.frame(unique(lineages.df.ref$Mut.aa))

colnames(to.col.ref)<-"Mut.aa"
to.col.ref$Lineage<-NA

for (i in 1:nrow(to.col.ref)) {
  
  if(length(which(lineages.df.ref$Mut.aa==to.col.ref$Mut.aa[i]))>1){
  
  lineages2<- unique(gsub("X$","",unlist(strsplit(lineages.df.ref$ID2[which(lineages.df.ref$Mut.aa==to.col.ref$Mut.aa[i])],"/") )))
  
  lineages2<-strsplit(lineages2,"\\.")
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
  if(is.na(new.lin)){
    newid<-paste(variant.temp$pangos,collapse =  "/")
  }else{
    newid<-paste(new.lin,"\\.X",sep = "")
  }
  to.col.ref$Lineage[i]<-paste(lineages.df.ref$ID2[which(lineages.df.ref$Mut.aa==to.col.ref$Mut.aa[i])],collapse = "/")
  if(length(unique(lineages.df.ref$ID2[which(lineages.df.ref$Mut.aa==to.col.ref$Mut.aa[i])]))==1){
    to.col.ref$Lineage[i]<-unique(lineages.df.ref$ID2[which(lineages.df.ref$Mut.aa==to.col.ref$Mut.aa[i])])
  }
  }else{
    to.col.ref$Lineage[i]<-lineages.df.ref$ID2[which(lineages.df.ref$Mut.aa==to.col.ref$Mut.aa[i])]
  }
  
}

ref.list<-strsplit(lineages.df.ref$Mut.aa, "/")
names(ref.list)<-lineages.df.ref$ID2


mut.c<-unlist(lapply(ref.list, length))
lineages.clean$PangoMismatches<-NA
lineages.clean$PangoDiscrepancies<-NA
lineages.clean$CleanMut.aa<-NA
PangoLineagesMatched<-NA

pb<-txtProgressBar(max = nrow(lineages.clean))
for(i in 1:nrow(lineages.clean)){
  setTxtProgressBar(pb,i)
  mut.aa<-unlist(strsplit(lineages.clean$Mut.aa[i], "/") )
  if(length(grep("X$",mut.aa))>0) mut.aa<-mut.aa[-grep("X$",mut.aa)]
  lineages.clean$CleanMut.aa[i]<-mut.aa
  match.vect<-vector()
  mismatch.vect<-vector()
  
  for (k in 1:length(ref.list)) {
    match.vect<- c(match.vect,length(which(mut.aa %in% ref.list[[k]] )))
    mismatch.vect<-c(mismatch.vect, length(which(table(c(mut.aa,ref.list[[k]]))==1)))  
  }
  #PerfectMatch
  if(length(which(mismatch.vect==0))==1){
    lineages.clean$PangoLineagesMatched[i]<-names(ref.list)[which(mismatch.vect==0)]
    lineages.clean$PangoMismatches[i]<-0
  #DiscrepancyMatch
  }else{
    pangonames<-names(ref.list)[which(mismatch.vect==min(mismatch.vect))]
    pangonames.max<-match.vect[which(mismatch.vect==min(mismatch.vect))]
    index.min<-which(mismatch.vect==min(mismatch.vect))
    
    if(length(which(pangonames.max==max(pangonames.max)))==1){
      index.min<-index.min[which(pangonames.max==max(pangonames.max))]
      mut.ref<-ref.list[[index.min]]
      lineages.clean$PangoLineagesMatched[i]<-pangonames[which(pangonames.max==max(pangonames.max))]
      lineages.clean$PangoMismatches[i]<-min(mismatch.vect)
    }else{
      mut.ref<-unique(unlist(ref.list[which(mismatch.vect==min(mismatch.vect))]))
      lineages.clean$PangoLineagesMatched[i]<-paste(unique(names(ref.list)[which(mismatch.vect==min(mismatch.vect))]),collapse = "/")
      lineages.clean$PangoMismatches[i]<-min(mismatch.vect)
    }
    
    plus<- mut.aa[-which(mut.aa %in% mut.ref)]
    minus<-mut.ref[-which(mut.ref %in% mut.aa)]
    if(length(plus)>0){
      plus<-paste("+",plus,sep = "")
      plus<-paste(plus,collapse = "")
      lineages.clean$PangoDiscrepancies[i]<-plus
    } 
    if(length(minus)>0){
      minus<-paste("-",minus,sep = "")
      minus<-paste(minus,collapse = "")
      if(!is.na(lineages.clean$PangoDiscrepancies[i])){
        lineages.clean$PangoDiscrepancies[i]<-paste(lineages.clean$PangoDiscrepancies[i], minus,collapse = "")
      }else{
        lineages.clean$PangoDiscrepancies[i]<-minus
      }
    } 
    
  }
}

if(length(which(is.na(lineages.clean$PangoDiscrepancies)))>0){
  lineages.clean$PangoDiscrepancies[which(is.na(lineages.clean$PangoDiscrepancies))]<-"NA"  
}

if(length(which(lineages.clean$PangoMismatches==0))>0){
  lineages.clean$PangoDiscrepancies[which(lineages.clean$PangoMismatches==0)]<-"None"  
}

if(length(which(lineages.clean$PangoMismatches>3))>0){
  lineages.clean$PangoLineagesMatched[which(lineages.clean$PangoMismatches>3)]<-"Unclassified"
  lineages.clean$PangoDiscrepancies[which(lineages.clean$PangoMismatches>3)]<-"NA"
  }


#Collapse lineages
lineages<-unique(lineages.clean$PangoLineagesMatched)
lineages<-lineages[grep("/", lineages)]

if(length(lineages)>0){
  for (i in 1:length(lineages)) {
    lineages2<- unique(gsub("\\.X$","",unlist(strsplit(lineages[i],"/") )))
    lineages2<-strsplit(lineages2,"\\.")
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
    if(is.na(new.lin)){
      newid<-lineages[i]
    }else{
      newid<-paste(new.lin,".X",sep = "")
    }
    lineages.clean$PangoLineagesMatched[which(lineages.clean$PangoLineagesMatched==lineages[i])]<-newid
  }  
  

}




