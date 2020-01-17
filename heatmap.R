## cytos vs drug regime ##

## load metadata file ####
metadata <- read.csv2('metadata.csv', header=T, row.names=1)
metadata.drugs <- metadata[,50:70]

metadata.drugs <- metadata.drugs[row.names(metadata.drugs)%in%unique(cytos.fibro$Replicate),]

cytos.fibro.replicatemean <- aggregate(cytos.fibro[, 4:55], list(cytos.fibro$Replicate), mean)
row.names(cytos.fibro.replicatemean) <- cytos.fibro.replicatemean$Group.1
cytos.fibro.replicatemean$Group.1 <- NULL
metadata.drugs <- metadata.drugs[row.names(metadata.drugs)%in%row.names(cytos.fibro.replicatemean),]
cytos.fibro.replicatemean <- cytos.fibro.replicatemean[row.names(cytos.fibro.replicatemean)%in%row.names(metadata.drugs),]

row.names(metadata.drugs)==row.names(cytos.fibro.replicatemean)


cytos.drugs <- data.frame(cytos.fibro.replicatemean, metadata.drugs)

for (x in seq_along(cytos.drugs[,53:73])){
  pdf(paste0(names(cytos.drugs)[52+x], '.pdf'))
    for (y in seq_along(cytos.drugs[,1:52])){
    print(ggplot(cytos.drugs, aes_string(x=names(cytos.drugs)[x+52], y=names(cytos.drugs)[y])) + geom_boxplot() + 
            geom_beeswarm() + stat_compare_means(method='t.test'))
    
    }
  dev.off()  
}


cytos.summ <- aggregate(abcam.cytos.log[,4:55],list(abcam.cytos.log$Replicate), mean)
row.names(cytos.summ) <- cytos.summ$Group.1
library(stringi); library(plyr)
classification <- data.frame('names'=cytos.summ$Group.1, 'class'=revalue(as.factor(sapply(cytos.summ$Group.1, function(x)(stri_sub(x, 8,8)))), c('C'='Control', 'F' = 'Fibromyalgia')))
cytos.summ$Group.1 <- NULL
row.names(cytos.summ) == classification$names
cytos.summ <- data.frame('Class' = classification$class, cytos.summ)
cytos.summ <-  cytos.summ[order(cytos.summ$Class),]

classification <- classification[order(classification$class),]
row.names(cytos.summ) == classification$names
row.names(classification) <- classification$names
annotations <- classification
annotations$names <- NULL

library(gplots)
heatmap.2(as.matrix(log(cytos.summ[,2:53])), distfun = dist,  Rowv=F, hclustfun = hclust,
           scale = 'column', trace = 'none', col=redblue(100), breaks = seq(-3,3, length.out = 101), key = T)

library(pheatmap)
pheatmap(as.matrix(log(cytos.summ[,2:53])), cluster_rows = F, cluster_cols = T, annotation_row = annotations, scale = 'column', gaps_row = 36, col=redblue(100))

write.xlsx(cytos.summ, 'cytos_summ.xlsx', sheetName = 'Sheet 1')
cytos.summ.fib <- subset(cytos.summ, cytos.summ$Class=='Fibromyalgia')
cytos.summ.ctl <- subset(cytos.summ, cytos.summ$Class!='Fibromyalgia')

pheatmap(as.matrix(log(cytos.summ.fib[,2:53])), cluster_rows = T, cluster_cols = F, scale = 'column', col=redblue(100), 
         main = 'Fibro cytokines',  fontsize_col = 7.5)#, filename = 'fibrocytosheatmap.pdf')
pheatmap(as.matrix(log(cytos.summ.ctl[,2:53])), cluster_rows = T, cluster_cols = F, scale = 'column', col=redblue(100), 
         main = 'Control cytokines', fontsize_col = 7.5, filename = 'controlcytosheatmap.pdf')

heatmap.2(as.matrix(log(cytos.summ.fib[,2:53])), distfun = dist,  Rowv=T, hclustfun = hclust,
          scale = 'column', trace = 'none', col=redblue(100), breaks = seq(-3,3, length.out = 101), key = T)
heatmap.2(as.matrix(log(cytos.summ.ctl[,2:53])), distfun = dist,  Rowv=T, hclustfun = hclust,
          scale = 'column', trace = 'none', col=redblue(100), breaks = seq(-3,3, length.out = 101), key = T)

cytos.info <- read.xlsx2('cytos_summ.xlsx', sheetIndex = 2, header = T, stringsAsFactors=F)
cytos.info$Function <- as.factor(cytos.info$Function)
cytos.info$Function <- ordered(cytos.info$Function, levels=c('Antiinflammatory', 'Proinflammatory', 'Others'))
cytos.info <- cytos.info[order(cytos.info$Function),]
col.order <- cytos.info$Function
names(col.order) <- cytos.info$Label
col.order <- as.data.frame(col.order)
names(col.order) <- 'CytoFunction'

  
cytos.summ1 <- cytos.summ[c(cytos.info$Names)]
cytos.summ1 <- data.frame('Class'=cytos.summ$Class, cytos.summ1) 
names(cytos.summ1)[2:53] == cytos.info$Names
names(cytos.summ1) <- c('Class', c(cytos.info$Label))
cytos.summ.fib <- subset(cytos.summ1, cytos.summ1$Class=='Fibromyalgia')
cytos.summ.ctl <- subset(cytos.summ1, cytos.summ1$Class!='Fibromyalgia')
pheatmap(as.matrix(log(cytos.summ.fib[,2:53])), cluster_rows = T, cluster_cols = F, scale = 'column', col=redblue(100), 
         main = 'Fibro cytokines',  fontsize_col = 7.5, annotation_col = col.order)#, filename = 'fibrocytosheatmap.pdf')
pheatmap(as.matrix(log(cytos.summ.ctl[,2:53])), cluster_rows = T, cluster_cols = F, scale = 'column', col=redblue(100), 
         main = 'Control cytokines', fontsize_col = 7.5, annotation_col = col.order) #, filename = 'controlcytosheatmap.pdf')

row.leg <- cytos.drugs[,53:73]
row.leg <- row.leg[,c('Benzodiazepina', 'Antidepresivo_SerNora', 'Morfico', 'Antiinflamatorio', 'Proton.Pump.Inhbitor')]
names(row.leg) <- c('Benzodiazepina', 'Antidepresivo', 'Morfico', 'Antiinflamatorio', 'Omeoprazol')
annot.col <-  list(
  CytoFunction = c(Antiinflammatory='#c50f01', Proinflammatory='#049afd', Others='#00ac07'),
  Benzodiazepina = c(NO='#c50f01', YES = '#049afd'),
  Antidepresivo = c(NO='#c50f01', YES = '#049afd'),
  Morfico = c(NO='#c50f01', YES = '#049afd'),
  Antiinflamatorio = c(NO='#c50f01', YES = '#049afd'),
  Omeoprazol = c(NO='#c50f01', YES = '#049afd'))


pheatmap(as.matrix(log(cytos.summ.fib[,2:53])), cluster_rows = T, cluster_cols = F, scale = 'column', col=redblue(100), 
         main = 'Fibro cytokines',  fontsize_col = 7.5, annotation_col = col.order, annotation_row = row.leg, 
         gaps_col = c(6,49), cutree_rows = 3, annotation_colors = annot.col, annotation_legend = F)#, filename = 'fibrocytosheatmap.pdf')

pheatmap(as.matrix(log(cytos.summ.ctl[,2:53])), cluster_rows = T, cluster_cols = F, scale = 'column', col=redblue(100), 
         main = 'Control cytokines', fontsize_col = 7.5, annotation_col = col.order, 
         gaps_col = c(6,49), cutree_rows = 4, annotation_colors = annot.col) #, filename = 'controlcytosheatmap.pdf')

pdf('heatmap_complete.pdf', width = 12, height = 12)
pheatmap(as.matrix(log(cytos.summ1[,2:53])), cluster_rows = T, cluster_cols = F, scale = 'column', col=redblue(100), 
         main = 'Cytokines expression',  fontsize_col = 9, annotation_col = col.order, 
         gaps_col = c(6,49),  annotation_legend = T, cutree_rows = 3)
dev.off()
