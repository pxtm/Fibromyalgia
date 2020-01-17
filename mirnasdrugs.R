## miRNAs distribution by drug regime ##

## metadata loading ####
metadata <- read.csv2('metadata.csv', header=T, row.names = 1)
metadata.drugs <- metadata[,50:70]
metadata.drugs <- metadata.drugs[row.names(metadata.drugs)%in%micros.fibro$Sample,]

micros.fibro <- micros.fibro[micros.fibro$Sample%in%row.names(metadata.drugs),]
row.names(metadata.drugs) == micros.fibro$Sample

micros.drugs <- data.frame(micros.fibro, metadata.drugs)

## boxplots of drugs vs micros ####
for (x in seq_along(micros.drugs[,52:72])){
  pdf(paste0(names(micros.drugs)[51+x], '.pdf'))
  for (y in seq_along(micros.drugs[,3:51])){
    print(ggplot(micros.drugs, aes_string(x=names(micros.drugs)[x+51], y=names(micros.drugs)[2+y])) + geom_boxplot() + 
            geom_beeswarm() + stat_compare_means(method='t.test'))
    
  }
  dev.off()  
}

## heatmap by drugs ####
library(pheatmap); library(gplots)
row.names(abcam.data.log) <- abcam.data.log$Sample

row.leg <- micros.drugs[,52:72]
row.leg <- row.leg[,c('Benzodiazepina', 'Antidepresivo_SerNora', 'Morfico', 'Antiinflamatorio', 'Proton.Pump.Inhbitor')]
names(row.leg) <- c('Benzodiazepina', 'Antidepresivo', 'Morfico', 'Antiinflamatorio', 'Omeoprazol')
annot.col <-  list(
  Benzodiazepina = c(NO='#c50f01', YES = '#049afd'),
  Antidepresivo = c(NO='#c50f01', YES = '#049afd'),
  Morfico = c(NO='#c50f01', YES = '#049afd'),
  Antiinflamatorio = c(NO='#c50f01', YES = '#049afd'),
  Omeoprazol = c(NO='#c50f01', YES = '#049afd'))

row.names(micros.drugs) <- micros.drugs$Sample

pdf('heatmapall.pdf', width = 12, height = 9)
pheatmap(as.matrix(abcam.data.log[,2:50]), scale = 'column', cluster_cols = T, cutree_rows = 3)
dev.off()

pdf('heatmapfibros.pdf', width = 12, height = 9)
pheatmap(as.matrix(micros.drugs[,3:51]), scale = 'column', cluster_cols = T, col=redblue(100), main='Fibro micros', fontsize_col = 7.5, 
         annotation_row = row.leg, annotation_colors = annot.col, cutree_rows = 3, annotation_legend = F)
dev.off()


## multivariate analysis ####
library(ropls); library(ggrepel)
abcam.data.log$Sample <- NULL
micros.pca <- opls(abcam.data.log[,2:50])
plot(micros.pca, typeVc='x-score', parAsColFcVn=abcam.data.log$Group)
micros.pls <- opls(abcam.data.log[,2:50], abcam.data.log$Group, predI=3)

loadings.pls <- as.data.frame(cbind('label'=rownames(getLoadingMN(micros.pls)), getLoadingMN(micros.pls)))
loadings.pls[,2:4] <- sapply(loadings.pls[,2:4], function(x)as.numeric(as.character(x)))
ggplot(loadings.pls, aes(p1, p2)) + geom_point() + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) +
  geom_text_repel(aes(label=label), hjust=1.05, yjust=0) + ylim(0.3, -0.3) + xlim(-0.4,0.4) + 
  theme(axis.line = element_line(size = 0.8, linetype = "solid"), axis.title = element_text(face = "bold"),
        axis.text = element_text(face = "bold"), panel.background = element_rect(fill = NA))

getSummaryDF(micros.pls)
VIPs <- as.data.frame(sort(getVipVn(micros.pls)))
names(VIPs) <- 'VIPscore'
VIPs$order <- factor(row.names(VIPs), levels=row.names(VIPs))

ggplot(VIPs, aes(x=VIPscore, y=order)) + geom_point(size=3) + 
  theme(axis.line = element_line(size = 0.3, linetype = "solid"), axis.text = element_text(size = 11, face = "bold", colour = "black"), 
        panel.background = element_rect(fill = NA)) +labs(x = "VIP score", y = "Variable") +
  theme(axis.title = element_text(size = 14,  face = "bold"), axis.ticks.y = element_blank()) + geom_segment(aes(yend=order), xend=0) +labs(y = NULL)
