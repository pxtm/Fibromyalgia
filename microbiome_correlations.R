## load microbiome data ####
otus <- read.csv2('otus_tax.csv', header = T, row.names = 1)
abcam.data.log
taxon <- data.frame(otus$Taxon)
row.names(taxon) <- row.names(otus)
otus <- otus[,names(otus)%in%abcam.data.log$Sample]
otus <- as.data.frame(t(otus))
otus1 <- otus[,colSums(otus[,3:833]==0)<24]
row.names(otus1) == abcam.data.log$Sample


otus.associations <- data.frame(abcam.data.log, otus1)
otus.associations$Sample <- NULL
otus.associations$Group <- NULL

library(corrplot)
otus.cors <- cor(otus.associations)
otus.pvals <- cor.mtest(otus.associations)
corrplot(otus.cors[1:49,50:65], p.mat = otus.pvals$p[1:49,50:65], insig = 'blank')

############################
library(corrplot)
otus.genus <- read.csv2('level6.csv', header = T, sep = '\t', row.names=1)
micros.log <- abcam.data.log
row.names(micros.log) <- micros.log$Sample
micros.log$Sample <- NULL ; micros.log$Group <- NULL
otus.genus <- otus.genus[,names(otus.genus)%in%row.names(micros.log)]
otus.genus <- as.data.frame(t(otus.genus))
otus.genus[] <- sapply(otus.genus, function(g)as.numeric(as.character(g)))
otus.genus <- otus.genus[, colSums(otus.genus!=0)>36]
row.names(otus.genus) == row.names(micros.log)
otus.genus <- otus.genus[match(row.names(micros.log), row.names(otus.genus)),]

otus.associations <- data.frame(micros.log, otus.genus)
otus.cors <- cor(otus.associations)
otus.pvals <- cor.mtest(otus.associations)
pdf('corrplot_microbiome_micros.pdf', height = 15, width=15)
corrplot(otus.cors[1:49,50:100], p.mat = otus.pvals$p[1:49, 50:100], insig = 'blank')
dev.off()

otus.associations$hsa.mir.335.5p
otus.associations$c__Bacilli.o__Bacillales.__.__
otus.associations$f__Enterobacteriaceae.g__Enterobacter
ggplot(otus.associations, aes(x=hsa.mir.335.5p, y=c__Bacilli.o__Bacillales.__.__)) + geom_point() + geom_smooth()
ggplot(otus.associations, aes(x=hsa.mir.335.5p, y=f__Enterobacteriaceae.g__Enterobacter)) + geom_point() + geom_smooth()

pdf('corrplot_microbiome_microsCLUSTERED.pdf', height = 25, width=25)
corrplot(otus.cors, p.mat=otus.pvals$p, insig = 'blank', order='hclust')
dev.off()
