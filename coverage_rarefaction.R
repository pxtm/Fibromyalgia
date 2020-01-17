# coverage 
## sum taxa per sample !=0 of both non_rarefied and rarefied datasets
rar_otus<-as.data.frame(t(otus))
rar_otus[]<-sapply(rar_otus, function(z){as.numeric(as.character(z))})
taxa_rar<-sapply(rar_otus, function(x){sum(x!=0)})

otus_norar<-read.csv2('non_rar_otus.csv', header=T, row.names = 1)
otus_norar<-otus_norar[,names(otus_norar)%in%row.names(metadata)]
taxa_norar<-sapply(otus_norar, function(x){sum(x!=0)})

## divide the rarefied taxa by the non-rarefied and compute basical stats data
rar_coverage<-taxa_rar/taxa_norar
sd(rar_coverage)
summary(rar_coverage)