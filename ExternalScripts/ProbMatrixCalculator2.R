library(seqinr)

refs<-read.fasta("/media/nacho/Data/wastewater/ReferencesProb/Spike28082023.fasta")

poi<-c(1250:2250)

#Cleaning identical spikes.

lineages<-unique(gsub(".*_","", names(refs)))

bases<-c("a","t","c","g","-")
#Generation of the probability matrix
#Parallel here
pb<-txtProgressBar(min = 1, max = length(poi),initial = 1)
for (i in 1:length(poi)) {
  setTxtProgressBar(pb,i)
  dummym<-matrix(0, ncol = length(lineages), nrow = 5 )
  
  base.m<-lapply(refs, function(x)x[poi[i]])  
  names(base.m)<-gsub(".*_","", names(base.m))
  base.m<-unlist(base.m)
  base.m2<-match(base.m, c("a","t","c","g","-"))
  names(base.m2)<-names(base.m)
  
  for (j in 1:length(lineages)) {
    for (k in c(1:5)) {
      P.L_M<-length(which(names(base.m2)==lineages[j] & base.m2==k))
      AP.L_M<-length(which(names(base.m2)!=lineages[j] & base.m2==k))
      
      #Priors based on data
      #P.L<-length(which(names(base.m2)==lineages[j]))/length(base.m2)
      #AP.L<- length(which(names(base.m2)!=lineages[j]))/length(base.m2)
      
      #No prior
      P.L<-1/length(lineages)
      AP.L<-(length(lineages)-1)/length(lineages)
      
      #No prior
      #P.L<-1/10
      #AP.L<-(10-1)/10
      
      #Include penalty
      
      #dummym[k,j]<-   (P.L_M * P.L) /((P.L_M * P.L)+(AP.L_M * AP.L))
      dummym[k,j]<-   P.L_M / (P.L_M + AP.L_M)
    }
  }
  dummym<-as.data.frame(dummym)
  colnames(dummym)<-lineages
  rownames(dummym)<-paste(poi[i],c("A","T","C","G","D"),sep = "")
  if(!exists("dummy.out")){ 
    dummy.out<-dummym
  }else{
    dummy.out<-rbind(dummy.out, dummym)
  }
}

#Normalization
dummy.out<-as.data.frame(apply(dummy.out, 2, FUN = function(x) x/sum(x,na.rm=TRUE)))

write.csv(dummy.out,"/media/nacho/Data/DockerImages/Wastewater_SARS-CoV-2/CommonFiles/reference/ProbMatrixSep2023.csv" )


refs.ohe<-list()
poi<-c(1250:2250)
refs.to.poi<-do.call(rbind,refs)
refs.to.poi<-refs.to.poi[,poi]

discriminant<-function(x){
  most.common<-names(table(x)[which(table(x)==max(table(x)))])[1]
  ratio<-length(which(x!=most.common))/length(x)
  return(ratio)
}
vari.refs<-apply(refs.to.poi, 2,FUN=discriminant )


refs.ohe<-do.call(rbind,refs.ohe)

vari.refs<-apply(refs.ohe,2,FUN=function(x)sum(x)/nrow(refs.ohe))
which(vari.refs>0.035)

write.csv(which(vari.refs>0.035),"/media/nacho/Data/DockerImages/Wastewater_SARS-CoV-2/CommonFiles/reference/NoiseRefs.csv",row.names = FALSE )
