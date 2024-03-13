library(readxl)

inputdir<-"/MountingPoint/"
#inputdir<-"/media/nacho/Data6G/wastewater_raw_fastq/"
outputdir<-"/Data/"
#outputdir<-"/media/nacho/Data6G/test/"

if(dir.exists(inputdir)){

input<-list.files(inputdir, pattern ="*.xlsx$",full.names = TRUE )
if(length(input)==1){

df<-read_xlsx(input)

if(length(which(df$Status=="FAIL"))){
  df<-df[-which(df$Status=="FAIL"),]
}

df$`NGS Oppsett`<-tolower(df$`NGS Oppsett`)

ww<-unique(df$`NGS Oppsett`)

for (i in 1:length(ww)) {
  temp.df<-df[which(df$`NGS Oppsett`==ww[i]),]
  temp.df$NewName<-paste(temp.df$`NGS Oppsett`, temp.df$`Geo-Navn`, sep = ".")
  dirs<-list.dirs(paste(inputdir,ww[i],sep = ""), full.names = TRUE)
   for (d in 1:nrow(temp.df)) {
     dir.create(paste(outputdir,temp.df$NewName[d],sep = ""))
     system(paste("ln -s ", inputdir, temp.df$`NGS Oppsett`[d] , "/", temp.df$Barcode[d], "/*.fastq.gz ", 
                  outputdir,temp.df$NewName[d],"/" ,sep = ""))
   }
}

}
}