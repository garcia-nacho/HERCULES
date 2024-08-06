#Corona Swiss-Army-Knife Docker Image
#Nacho Garcia 2021 / iggl@fhi.no

library(msa)
library(seqinr)

#Docker

spike <-"/home/docker/CommonFiles/reference/SpikeRef.fa"
refs<-"/Data/variantRefs.fasta"
sequences<-readDNAStringSet(c(refs,spike))
  
  sequences.aln<-list()
  pb<-txtProgressBar(min = 1, max = length(sequences), initial = 1)
  for (i in 1:(length(sequences)-1)) {
    
    setTxtProgressBar(pb, i)    
    seqs.to.aln<-sequences[c(i, grep("Spike", names(sequences)))]
    alignment<-msa(seqs.to.aln, "Muscle")
    x<-DNAMultipleAlignment(alignment)
    DNAStr = as(x, "DNAStringSet")
    Aligned.samples <-  unlist(base::strsplit(as.character(DNAStr[grep("Spike", names(seqs.to.aln))]),"") )
    Query <-  unlist(base::strsplit(as.character(DNAStr[-grep("Spike", names(seqs.to.aln))]),"") )
    start.read<-which(Aligned.samples!="-")
    start.read<-start.read[which(start.read==min(start.read))]
    
    end.read<-which(Aligned.samples!="-")
    end.read<-end.read[which(end.read==max(end.read))]
    
    Aligned.samples<-as.character(Query[start.read:end.read])
    if(length(which(Aligned.samples=="-"))>0) Aligned.samples<-Aligned.samples[-which(Aligned.samples=="-")]
    sequences.aln[[i]]<-Aligned.samples
    names(sequences.aln)[i]<-names(DNAStr)[-grep("Spike", names(DNAStr))]
    }
  
  write.fasta(sequences = sequences.aln, names = paste(names(sequences.aln),"_Spike",sep = ""),
              "/Data/VariantSpike.fasta")


