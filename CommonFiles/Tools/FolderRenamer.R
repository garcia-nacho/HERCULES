library(readxl)


files<-list.files(pattern = "xlsx")
if(length(files)==1){
  df<-read_xlsx(files)
  df$Date<-gsub("/","-",df$Date)

  all.files<-list.files(recursive = TRUE)
  barcodes<-gsub("/.*", "", all.files)
  
  for (i in 1:nrow(df)) {
    if(length(which(barcodes == df$Barcode[i]))==1){
      dir.create(paste(df$Experiment[i],".", df$Date[i],"_",df$Location[i],sep = ""))
      oldname<-all.files[which(barcodes == df$Barcode[i])]
      newname<-paste(gsub("/.*","",oldname),"/","read",c(1:length(oldname)),".fastq.gz",sep = "")
      newname<-gsub(".*/",paste(df$Experiment[i],".", df$Date[i],"_",df$Location[i],"/",sep = ""),newname)
      file.rename(oldname,newname)
      
      file.remove(df$Barcode[i])
    }
  }
  
}



