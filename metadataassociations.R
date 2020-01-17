## pain correlations 
pain <- read.csv2('pain_indices.csv', header=T, stringsAsFactors = F)
cytos.fibro.replicatemean$Group.1 <- as.character(cytos.fibro.replicatemean$Group.1) 
pain <- pain[pain$X%in%cytos.fibro.replicatemean$Group.1,]
pain <- pain[order(pain$X), ]
rownames(pain) <- pain$X
pain$X <- NULL
cytos.fibro.replicatemean <- cytos.fibro.replicatemean[order(cytos.fibro.replicatemean$Group.1)]
rownames(cytos.fibro.replicatemean) <- cytos.fibro.replicatemean$Group.1
cytos.fibro.replicatemean$Group.1 <- NULL
row.names(pain) == row.names(cytos.fibro.replicatemean)

library(corrplot)
pain.corr <- data.frame(pain, cytos.fibro.replicatemean)

for (i in seq_along(pain.corr[,1:8])){
  pdf(paste0('corr_', names(pain.corr)[i], '.pdf'))
  for (a in seq_along(pain.corr[,9:43]))
    print(ggplot(pain.corr, aes_string(x=names(pain.corr)[i], y=names(pain.corr)[8+a])) + geom_smooth(method='lm') + geom_point()+stat_cor(method='pearson'))
  dev.off()
}

pain.corr2 <- pain.corr
pain.corr2 <- sapply(pain.corr2, NAfunc)
rename <- colnames(pain.corr2)
rename[9:43] <- sapply(rename[9:43], function(b)(substr(b, 1, nchar(b)-8)))
rename <- sapply(rename, function(c)gsub('\\.', '-', c))
rename[2] <- gsub('fatiga', 'fatigue', rename[2])
rename[3] <- gsub('sueño', 'sleepness', rename[3])
rename[4] <- gsub('cognitivos', 'cognitive', rename[4])
colnames(pain.corr2) <- rename
pain.corr2 <- cor(pain.corr2)
res.corr <- cor.mtest(pain.corr2)

tiff('corrplot_pain.tiff', width = 3500, height = 5500, res=600, units = 'px')
corrplot(pain.corr2[9:43,1:8], p.mat = res.corr$p[9:43,1:8], insig = 'blank', tl.col='black', cl.pos='n')
dev.off()

par(mfrow=c(1,1))
tiff('corrplot_all.tiff', width = 5000, height = 5000, res=600, units = 'px')
corrplot(pain.corr2, p.mat = res.corr$p, insig = 'blank', order='hclust')
dev.off()


pain.corr
for (i in seq_along(metadata)){
  levels(metadata[,i])[levels(metadata[,i])==''] <- NA
}
metadata$Age <- NULL
metadata$Status1 <- NULL
metadata$Status2 <- NULL
metadata$Subject <- NULL
metadata.assoc <- metadata[row.names(metadata)%in%row.names(pain.corr),]
pain.assoc <- pain.corr[,c('CXCL5..pg.mL.', 'G.CSF..pg.mL.')] 
pain.assoc <- pain.assoc[row.names(pain.assoc)%in%row.names(metadata.assoc),]
row.names(pain.assoc) == row.names(metadata.assoc)

pain.assoc <- data.frame(pain.assoc, metadata.assoc)

sink('CXCL5_associations.txt')
for (x in seq_along(pain.assoc[,3:68])){
  if (length(levels(pain.assoc[,x]))==2){
    tryCatch({print(paste0(names(pain.assoc[x]), ' t-test: ', 
                           t.test(pain.assoc$CXCL5..pg.mL.~pain.assoc[,x], na.action = na.omit)$p.val))
    }, error=function(a){cat(names(pain.assoc[x]), ' ERROR: ', conditionMessage(a), '\n')})
  }
  else if (length(levels(pain.assoc[,x]))==3){
    tryCatch({
      aovmodel <- aov(CXCL5..pg.mL.~pain.assoc[,x], data=pain.assoc);
      print(paste0(names(pain.assoc)[x], ' ANOVA: ', summary(aovmodel)[[1]][['Pr(>F)']][1]));
      tryCatch({
        print(TukeyHSD(aovmodel)[[1]][,4])
      }, error=function(b){cat('ERROR: ', conditionMessage(b), '\n')})
    }, error=function(a){cat('ERROR: ', conditionMessage(a), '\n')})
  }
  else {}
}
sink()

sink('GCSF_associations.txt')
for (x in seq_along(pain.assoc[,3:68])){
  if (length(levels(pain.assoc[,x]))==2){
    tryCatch({print(paste0(names(pain.assoc[x]), ' t-test: ', 
                           t.test(pain.assoc$G.CSF..pg.mL.~pain.assoc[,x], na.action = na.omit)$p.val))
    }, error=function(a){cat(names(pain.assoc[x]), ' ERROR: ', conditionMessage(a), '\n')})
  }
  else if (length(levels(pain.assoc[,x]))==3){
    tryCatch({
      aovmodel <- aov(G.CSF..pg.mL.~pain.assoc[,x], data=pain.assoc);
      print(paste0(names(pain.assoc)[x], ' ANOVA: ', summary(aovmodel)[[1]][['Pr(>F)']][1]));
      tryCatch({
        print(TukeyHSD(aovmodel)[[1]][,4])
      }, error=function(b){cat('ERROR: ', conditionMessage(b), '\n')})
    }, error=function(a){cat('ERROR: ', conditionMessage(a), '\n')})
  }
  else {}
}
sink()

tiff('cytosassociations.tiff', height = 4500, width = 1500, res=600, units = 'px')
multiplot(
  ggplot(data=subset(pain.assoc, !is.na(Vegetables)), aes(x=Vegetables, y=CXCL5..pg.mL.)) + geom_boxplot() + geom_beeswarm() + 
    stat_compare_means(method='t.test', label.x = 1, label.y = 0.5, aes(label=paste0('p= ', round(..p.adj.., 3))),fontface='bold.italic' ) + 
    theme(axis.line = element_line(size = 0.6, linetype = "solid"), axis.ticks = element_line(colour = "gray16"),
          axis.title = element_text(size = 13, face = "bold"), axis.text = element_text(size = 11, face = "bold", colour = "black"),
          panel.background = element_rect(fill = NA)) +labs(y = "CXCL5 (pg/mL)"),
  ggplot(data=subset(pain.assoc, !is.na(Urate)), aes(x=Urate, y=CXCL5..pg.mL.)) + geom_boxplot() + geom_beeswarm() + 
    stat_compare_means(method='anova', label.y = 4, aes(label=paste0('p= ', round(..p.adj.., 3))),fontface='bold.italic') + 
    stat_compare_means(comparisons = list(c('Below Range', 'In Range'), c('Below Range', 'Over Range'), c('In Range', 'Over Range')), aes(label=paste0('p= ', round(..p.adj.., 3))),fontface='bold.italic') + 
    theme(axis.line = element_line(size = 0.6, linetype = "solid"), axis.ticks = element_line(colour = "gray16"), 
          axis.title = element_text(size = 13, face = "bold"), axis.text = element_text(size = 11, face = "bold", colour = "black"), 
          panel.background = element_rect(fill = NA)) +labs(y = "CXCL5 (pg/mL"), 
  ggplot(data=subset(pain.assoc, !is.na(Antiinflamatorio)), aes(x=Antiinflamatorio, y=G.CSF..pg.mL.)) + geom_boxplot() + geom_beeswarm() + 
    stat_compare_means(method='t.test', label.x = 1, label.y = 0.5, aes(label=paste0('p= ', round(..p.adj.., 3))),fontface='bold.italic' ) + 
    theme(axis.line = element_line(size = 0.6, linetype = "solid"), axis.ticks = element_line(colour = "gray16"),
          axis.title = element_text(size = 13, face = "bold"), axis.text = element_text(size = 11, face = "bold", colour = "black"),
          panel.background = element_rect(fill = NA)) +labs(y = "G-CSF (pg/mL)"),
  cols=1)
dev.off()
