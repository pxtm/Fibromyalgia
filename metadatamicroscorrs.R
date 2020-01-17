## pain correlations ####
pain <- read.csv2('pain_indices.csv', header = T, stringsAsFactors = F)
pain <- pain[pain$X%in%micros.fibro$Sample,]
pain <- pain[order(pain$X),]
rownames(pain) <- pain$X
pain$X <- NULL
micros.fibro <- micros.fibro[order(micros.fibro$Sample),]
row.names(micros.fibro) <- micros.fibro$Sample
micros.fibro$Sample <- NULL
micros.fibro$Group <- NULL
row.names(pain) == row.names(micros.fibro)

library(corrplot)
pain.corr <- data.frame(pain, micros.fibro)
names(pain.corr) <- sapply(names(pain.corr), function(c)gsub('\\.', '-', c))
names(pain.corr)[2] <- gsub('fatiga', 'fatigue', names(pain.corr[2]))
names(pain.corr)[3] <- gsub('sueño', 'sleepness', names(pain.corr[3]))
names(pain.corr)[4] <- gsub('cognitivos', 'cognitive', names(pain.corr[4]))



for (i in seq_along(pain.corr[,1:8])){
  pdf(paste0('corr_', names(pain.corr)[i], '.pdf'))
  for (a in seq_along(pain.corr[,9:57]))
    print(ggplot(pain.corr, aes_string(x=names(pain.corr)[i], y=names(pain.corr)[8+a])) + geom_smooth(method='lm') + geom_point()+stat_cor(method='pearson'))
  dev.off()
}

pain.corr2 <- pain.corr
pain.corr2 <- sapply(pain.corr2, NAfunc)
pain.corr2 <- cor(pain.corr2)
res.corr <- cor.mtest(pain.corr2)

par(mfrow=c(1,1))
tiff('micros_corr.tiff', width = 3500, height = 5500, res=600, units = 'px')
corrplot(pain.corr2[9:57,1:8], p.mat = res.corr$p[9:57,1:8], insig = 'blank', cl.pos='n', tl.col='black')
dev.off()

tiff('corrplot_all.tiff', width = 5000, height = 5000, res=600, units = 'px')
corrplot(pain.corr2, p.mat = res.corr$p, insig = 'blank', order='hclust')#, cl.pos='b')
dev.off()

## metadata correlations ####
metadata.forms <- metadata[,1:48]
metadata.forms <- metadata.forms[row.names(metadata.forms)%in%row.names(abcam.data.log), ]
micros.metadata <- abcam.data.log[row.names(abcam.data.log)%in%row.names(metadata.forms),]
row.names(metadata.forms) == row.names(micros.metadata) 
metadata.forms <- metadata.forms[,!(names(metadata.forms)%in%c('Status2', 'Status1', 'Subject', 'Age'))]
micros.metadata$Group <- NULL

metadata.forms2<- sapply(metadata.forms,function(a)as.numeric(a))
row.names(metadata.forms2) <- row.names(metadata.forms)
metadata.corr <- data.frame(metadata.forms2, micros.fibro)
metadata.corr2 <- cor(metadata.corr)
res.corr <- cor.mtest(metadata.corr)

tiff('corrplot_metadata.tiff', width = 6000, height = 6000, res=600, units='px')
corrplot(metadata.corr2[1:45,46:94], p.mat = res.corr$p[1:45,46:94], insig = 'blank')#, order='hclust')
dev.off()

## metadata associations with miRNAs ####
metadata.bps <- data.frame(metadata.forms, micros.metadata)
for (x in seq_along(metadata.bps[,1:44])){
  pdf(paste0('metadata/', names(metadata.bps)[x], '_corr.pdf'))
  for (y in seq_along(metadata.bps[,45:93]))
    print(ggplot(metadata.bps, aes_string(x=names(metadata.bps)[x], y=names(metadata.bps)[44+y])) + geom_boxplot() + geom_beeswarm()+
            stat_compare_means())
  dev.off()
}
summary(metadata.bps[,1:44])


## correlations micro-cytokines ####
cytos.summ <- readRDS('cytosdata.rds')
cytos.summ <- aggregate(cytos.summ[,4:55], list(cytos.summ$Replicate), mean)
rownames(cytos.summ) <- cytos.summ$Group.1
cytos.summ$Group.1 <- NULL
micros.cytos <- abcam.data.log[order(abcam.data.log$Group),]
row.names(cytos.summ) == row.names(micros.cytos)
cytos.summ <-  cytos.summ[match(rownames(micros.cytos), rownames(cytos.summ)),]
row.names(cytos.summ) == row.names(micros.cytos)

micros.cytos <- data.frame(cytos.summ, micros.cytos)
row.names(micros.cytos) <- row.names(cytos.summ)
micros.cytos$Class <- NULL
micros.cytos$Group <- NULL
micros.cytos <- sapply(micros.cytos, NAfunc)
rename <- colnames(micros.cytos)
rename[1:51] <- sapply(rename[1:51], function(b)(substr(b, 1, nchar(b)-8)))
rename <- sapply(rename, function(c)gsub('\\.', '-', c))
colnames(micros.cytos) <- rename
micros.cytos.corr <- cor(micros.cytos)
micros.cytos.res <- cor.mtest(micros.cytos)


tiff('microsvscytos.tiff', height = 10000, width = 10000, res=600, units='px')
corrplot(micros.cytos.corr, p.mat = micros.cytos.res$p, insig = 'blank', tl.col='black', cl.pos = 'b')#, order='hclust')
dev.off()
tiff('microsvscytos_sepLABMEETING.tiff', height = 5500, width = 5500, res=600, units='px')
corrplot(micros.cytos.corr[53:101,1:52], p.mat = micros.cytos.res$p[53:101,1:52], insig = 'blank')
dev.off()
