## MICROS FIBRO DATA ANALYSIS ##
# Libraries required ####
library(ggplot2); library(ggpubr); library(ggthemes);  library(ggbeeswarm); library(xlsx); library(UnvaR)

## Data upload and group generation ####
abcam.data <- read.csv2('Data/NORMALIZED_DATA_working.csv', header=T, dec='.')
abcam.data <- abcam.data[-73,]
sapply(abcam.data[,3:70], function(x)sum(x < 0.56)) ## count how many values per column don't reach the LOD, which is 0.56 for this test
lod.all <- sapply(abcam.data[,3:70], function(x)(sum(x < 0.56)/72)*100)
abcam.data.ctl <- subset(abcam.data, abcam.data$Group=='Control')
abcam.data.fib <- subset(abcam.data, abcam.data$Group!='Control')
lod.ctl <- sapply(abcam.data.ctl[,3:70], function(x)(sum(x < 0.56)/36)*100)
lod.fib <- sapply(abcam.data.fib[,3:70], function(x)(sum(x < 0.56)/36)*100)
lod.samples <- data.frame(lod.all, lod.fib, lod.ctl)
#write.xlsx(lod.samples, 'lod_samples.xlsx', sheetName = 'LOD recount')
remove <- names(lod.all[lod.all>30]) ## remove all those miRNAs with >30% of values below LOD
abcam.data <- abcam.data[,!(names(abcam.data)%in%remove)] 

## normality assesment ####
shapiro.micros <- sapply(abcam.data[,3:51], function(x)shapiro.test(x)$p.val)
pdf('normality_assesment.pdf')
for (i in seq_along(abcam.data[,3:51])){
  print(ggqqplot(abcam.data[,2+i], title = paste(names(abcam.data[2+i]), ' shapiro ', shapiro.micros[i])))
  print(hist(abcam.data[,2+i]))
}
dev.off()

## log-renormalization ####
abcam.data.log <- abcam.data
abcam.data.log[,3:51] <- sapply(abcam.data.log[,3:51], function(x)(log(x)))
shapiro.micros.log <- sapply(abcam.data.log[,3:51], function(x)shapiro.test(x)$p.val)
pdf('qqplots_micros_log.pdf')
for (i in seq_along(abcam.data.log[,3:51])){
  print(ggqqplot(abcam.data.log[,2+i], title = paste(names(abcam.data.log[2+i]), ' shapiro ', shapiro.micros.log[i])))
  print(hist(abcam.data.log[,2+i]))
}
dev.off()

## Boxplots ####
pdf('boxplots.pdf')
for (i in seq_along(abcam.data.log[,3:51])){
  tiff(paste0('sig_boxplots/',names(abcam.data.log[2+i]), '.tiff'), res=96, height = 480, width = 320)
  print(ggplot(abcam.data.log, aes_string(x='Group', y=names(abcam.data.log[2+i])))+geom_boxplot() + geom_beeswarm()+ 
    labs(title=names(abcam.data[2+i]), y='Value') + stat_compare_means(method='t.test')+theme_bw())
  dev.off()
}
dev.off()


## Volcano plot ####
abcam.data.log <- do.call(data.frame,lapply(abcam.data.log, function(x) replace(x, is.infinite(x),0)))

micros.control <- subset(abcam.data.log, abcam.data.log$Group=='Control')
micros.fibro <- subset(abcam.data.log, abcam.data.log$Group!='Control')
control.mean <- sapply(micros.control[,3:51], function(x)(mean(x, na.rm=T)))
fibro.mean <- sapply(micros.fibro[,3:51], function(x)(mean(x, na.rm=T)))

fold.change <- fibro.mean/control.mean
log2.fc <- log2(fold.change)
pvals <- sapply(abcam.data.log[,3:51], function(y)t.test(y~Group, data=abcam.data.log)$p.val)
log10.pvals <- -log10(pvals)

volcano <- data.frame('name'= names(log2.fc), log2.fc, log10.pvals, 'colour'=cut(log2.fc, c(-Inf, -.51, .37, Inf)))
tiff('volcano_micros1.tiff', res=600, height = 4500, width=4500)
ggplot(volcano, aes(log2.fc, log10.pvals)) + geom_point(aes(colour=ifelse(log10.pvals>1.30103, colour, '')))+
  xlab(bquote(~log[2]~ ' Fold Change'))+ ylab(bquote(~-log[10]~ ' p-values'))+
#    geom_text(aes(label=ifelse(log10.pvals>1.30103, as.character(name), '')), hjust=1.05, yjust=0)+
  geom_hline(yintercept=-log10(.05),linetype="dashed",col="red")+
  annotate("text",label="p-val=.05",y=1.35,x=-1.5,fontface=4)+
  geom_hline(yintercept=2,linetype="dashed",col="orange")+annotate("text",label="p-val=.01",fontface=4,y=2.05,x=-1.5)+
  geom_hline(yintercept=3,linetype="dashed",col="green")+annotate("text",label="p-val=.001",y=3.05,x=-1.5,fontface=4)+
  geom_vline(xintercept=0,linetype="dashed",col="black")+theme_base(base_size=17)+
  scale_color_manual(values=c('black', 'darkolivegreen3', 'black', 'darkolivegreen3')) + guides(colour=F) + xlim(-1.6,1.6)
dev.off()
