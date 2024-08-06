library(keras)
library(data.table)
library(ggplot2)

#Initialization of layers

df<-fread("/media/nacho/Data/wastewater/TestW14/results/QC/45_Oslo_MashMatrix.csv")
df<-as.data.frame(df)
activation<-"sigmoid"

input.size<-ncol(df)-1

#neg<- 1-as.matrix(df[,-ncol(df)])

# Input <- layer_input(shape = input.size)
# 
# L1<- layer_dense( units = input.size, activation = activation)
# L2<- layer_dense( units = 300, activation = activation)
# L3<- layer_dense( units = 50, activation = activation)
# L4<- layer_dense( units = 10, activation = activation)
# 
# LS<- layer_dense( units = 2, activation = activation, name = "LatentSpace")
# 
# L5<- layer_dense( units = 10, activation = activation)
# L6<- layer_dense( units = 50, activation = activation)
# L7<- layer_dense( units = 300, activation = activation)
# L8<- layer_dense( units = input.size, activation = activation)
# 
# 
# #Concatenation of layers 
# IO1<-L1(Input)
# IO4<-L2(IO1)
# IO5<-L3(IO4)
# IO6<-L4(IO5)
# LSO<-LS(IO6)
# IO7<-L5(LSO)
# IO8<-L6(IO7)
# IO9<-L7(IO8)
# Output<-L8(IO9)
  
act.f<-"sigmoid"
input <- layer_input(shape = 3425)
aencoder = input %>% 
  layer_dense(units=3425, activation = act.f) %>% 
  layer_dense(units=600, activation = act.f) %>%
  layer_dense(units=60, activation = act.f) %>% 
  layer_dense(units=10, activation = act.f) %>% 
  layer_dense(units=2, activation = act.f, name = "LatentSpace") %>% 
  layer_dense(units=10, activation = act.f) %>%
  layer_dense(units=60, activation = act.f) %>%
  layer_dense(units=600, activation = act.f) %>%
  layer_dense(units=3425, activation = act.f)

autoencoder = keras_model(input, aencoder)

summary(autoencoder)

autoencoder %>% compile(
  loss = 'mean_squared_error',
  optimizer =optimizer_adagrad(),
  metrics = c('mean_squared_error')
)

ES<-list(callback_early_stopping(
  monitor = "mean_squared_error",
  min_delta = 0,
  patience = 8,
  verbose = 0,
  mode =  "min",
  restore_best_weights = TRUE
))


history <- autoencoder %>% fit(
  as.matrix(df[,-length(df)]),as.matrix(df[,-length(df)]) , 
  epochs = 100, batch_size = 128, 
  callbacks=ES,
  shuffle=TRUE,
  validation_split = 0.1
)

layer_name <- 'LatentSpace'
projection <- keras_model(inputs = autoencoder$input,
                          outputs = get_layer(autoencoder, layer_name)$output)


#predictions<-projection %>% predict(as.matrix(neg))
predictions<-projection %>% predict(as.matrix(df[,-length(df)]))
predictions<-as.data.frame(predictions)
predictions$Lineage<-df$Lineage

ggplot(predictions)+
  geom_point(aes(V1,V2, col=Lineage),alpha=0.1)+
  theme_minimal()+
  facet_wrap(~Lineage)

ggplot(predictions)+
  geom_jitter(aes(V1,V2, col=Lineage),alpha=0.1)+
  theme_minimal()+
  ylim(0.599,0.6)+
  xlim(0.6820,0.6822)+
  facet_wrap(~Lineage)


files<-list.files("/media/nacho/Data/wastewater/TestW14/results/QC/",pattern = "MashMatrix",full.names = TRUE)
pb<-txtProgressBar(min = 1, max = length(files),initial = 1)
for (i in 1:length(files)) {
  gc()
  setTxtProgressBar(pb,i)
  df<-fread(files[i])
  df<-as.data.frame(df)
  dummy<-projection %>% predict(as.matrix(df[,-length(df)]))
  dummy<-as.data.frame(dummy)
  dummy$Lineage<-df$Lineage
  dummy$Sample<-gsub(".*/","",gsub("_MashMatrix.*","",files[i]))
  
  if(!exists("df.out")){
    df.out<-dummy
  }else{
    df.out<-rbind(df.out,dummy)
  }
  rm(dummy)
  rm(df)
}

ggplot(df.out)+
  geom_jitter(aes(V1,V2,col=Lineage),alpha=0.01)+
  theme_minimal()+
  facet_wrap(~Sample)

table(df.out$Lineage[ which(df.out$V1>0.3212)])

#
df<-fread("/media/nacho/Data/wastewater/TestW14/results/QC/45_Oslo_MashMatrix.csv")
df<-as.data.frame(df)

deviation<-apply(df[,-ncol(df)], 2, mean )
names(deviation)[which(deviation==max(deviation))]
names(deviation)[which(deviation==min(deviation))]

names(deviation)[order(deviation ,decreasing = TRUE)][c(1:2)]

ggplot(df)+
  geom_jitter(aes(V631,V3037,col=Lineage),alpha=0.01, pch=".")+
  theme_minimal()

pca<-prcomp(df[,-ncol(df)])
