setwd("C:/Users/mcgma/Desktop")
library(xlsx)
qpcr.data <- read.xlsx2('qPCRfull.xlsx', header=T, sheetIndex = 1)
qpcr.data[,3:8] <- sapply(qpcr.data[,3:8], function(x)as.numeric(as.character(x)))
library(ggplot2); library(ggpubr)

## multiplot function ####
multiplot<-function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

## plotting ####

plots <- list()
for (i in (names(qpcr.data)[3:8])){
  plots[[i]] <- multiplot(
  ggplot(qpcr.data, aes_string(x=i, group="Status", fill=as.factor(qpcr.data$Status))) + geom_density(position='identity', alpha=.5) + 
    scale_fill_manual(name='Status', values=c('#00BFC4', '#F8766D')) + theme(plot.subtitle = element_text(vjust = 1), 
    plot.caption = element_text(vjust = 1), legend.position = "none", axis.line = element_line(linetype = "solid"), 
    panel.grid.major = element_line(colour = NA), panel.grid.minor = element_line(colour = NA), 
    axis.title = element_text(size = 13, face = "bold"), 
    axis.text = element_text(size = 12, face = "bold", colour = "black"), plot.title = element_text(face = "bold"), 
    panel.background = element_rect(fill = NA))+labs(x = paste0(i, " Cts")),
  
  ggplot(qpcr.data, aes_string(y=i, x="Status",  fill=as.factor(qpcr.data$Status))) + geom_boxplot(alpha=.5) + 
    stat_compare_means(method='t.test', aes(label=paste0('p= ', round(..p.adj.., digits = 3))), fontface=4) + 
    scale_fill_manual(name='Status', values=c('#00BFC4', '#F8766D')) +
    theme(plot.subtitle = element_text(vjust = 1), plot.caption = element_text(vjust = 1), 
          panel.grid.major = element_line(linetype = "blank"), panel.grid.minor = element_line(linetype = "blank"), 
          axis.title = element_text(size = 15, face = "bold"), 
          axis.text = element_text(size = 12, face = "bold", colour = "black"), 
          legend.text = element_text(size = 12, face = "bold"), legend.title = element_text(size = 12, face = "bold"), 
          panel.background = element_rect(fill = NA), legend.key = element_rect(fill = NA), 
          legend.background = element_rect(fill = NA), axis.line = element_line(linetype = "solid")) +
    labs(fill = "Status", y = paste0(i, " Cts")),
  cols = 2)
}



## relative abundance plot ####
healthy <- subset(qpcr.data, qpcr.data$Status=='Control')
fibro <- subset(qpcr.data, qpcr.data$Status=='Fibromyalgia')

means.genesH <- sapply(healthy[,3:8], function(x)median(x, na.rm=T))
means.genesF <- sapply(fibro[,3:8], function(x)median(x, na.rm=T))

#relab <- 2^-(ct-means.genes[i])
relab = function(ct, meangene){
  return(2^-(ct-meangene))
}

relab2 = function(ctF, ctH){
  return(2^-(ctF-ctH))
}

# relab(healthy[1,4], means.genes[1])
# 
# a=data.frame(qpcr.data$Status, sapply(qpcr.data[,3], function(x)relab(x, means.genes[1])))
# ggplot(a, aes(y=a[,2], x=a[,1],  fill=as.factor(a[,1]))) + geom_boxplot(alpha=.5) + geom_jitter()+
#   stat_compare_means(method='t.test', aes(label=paste0('p= ', round(..p.adj.., digits = 3))), fontface=4) +
#  ylim(0,1)


## RE-DID WITH SEPARATION BY PLATES ####
setwd('qpcrs_384_good')
qpcr <- read.xlsx('qpcr_fullGOOD.xlsx', sheetIndex = 1, header=T)
sapply(qpcr, class)
qpcr <- subset(qpcr, qpcr$Status != 'Negative') ## remove the negative controls
qpcr$Status <- droplevels(qpcr$Status)
split(qpcr, qpcr$Plate)
ggplot(qpcr, aes(x=Status, y=gadC, fill=Status)) + geom_boxplot(alpha=.5) + facet_grid(rows = vars(Plate))

library(reshape2)
qpcr.m <- melt(qpcr, id.vars = c('Sample', 'Plate', 'Status'), measure.vars = c('gadC', 'glnA', 'glsA1', 'glsA2', 'gadB1', 'gadB2'))
qpcr.m$Status <- ordered(qpcr.m$Status, levels = c('Positive', 'C', 'F'))
qpcr.m$Plate <- as.factor(qpcr.m$Plate)
levels(qpcr.m$Plate) <- c('Plate 1', 'Plate 2', 'Plate 3', 'Plate 4')
windows()
ggplot(qpcr.m, aes(x=variable, y=value, fill=Status)) + geom_boxplot(alpha=.5) + facet_grid(rows = vars(Plate)) + 
  stat_compare_means(method='t.test', comparisons = c('C', 'F')) + labs(x = "Genes", y = "Cts") +
  scale_fill_manual(name='Status', values=c('#FFFFFF', '#00BFC4', '#F8766D')) + 
  geom_segment(aes(x=1, xend=1.3, y=40, yend=40)) +

  
  theme(plot.subtitle = element_text(vjust = 1), plot.caption = element_text(vjust = 1), axis.line = element_line(linetype = "solid"), 
         axis.ticks = element_line(colour = "gray0"), axis.title = element_text(size = 16, face = "bold"), 
         axis.text = element_text(size = 13, face = "bold", colour = "black"), legend.text = element_text(size = 12, face = "bold"), 
         legend.title = element_text(size = 12, face = "bold"), panel.background = element_rect(fill = NA), 
         legend.key = element_rect(fill = NA), legend.background = element_rect(fill = NA), strip.background = element_rect(fill='white'), 
        strip.text.y = element_text(size=20))

qpcr.pvals <- subset(qpcr, qpcr$Status!='Positive')    
qpcr.pvals$Status <- droplevels(qpcr.pvals$Status)
qpcr.pvals <- split(qpcr.pvals, qpcr.pvals$Plate)

qpcr.pvals.labels <- data.frame()
for (i in length(qpcr.pvals)){
  sapply(qpcr.pvals[[i]][,5:10], function(x){
    print(print(p.adjust(t.test(x~Status, data=qpcr.pvals[[i]])$p.val, method = 'bonferroni')))
  })
  # for (y in seq(from=5, to=10)){
  #   print(p.adjust(t.test(qpcr.pvals[[i]][,y]~Status, data=qpcr.pvals[[i]])$p.val, method = 'bonferroni'))
  # }
}
      

# ggplot(qpcr.m, aes(x=Status, y=value, fill=Status)) + geom_boxplot(alpha=.5) + facet_grid(Plate~variable) +  
#   stat_compare_means(method='t.test', comparisons = list(c('C', 'F')), aes(label=paste0('p= ', round(..p.adj.., digits = 3))), fontface=4) +
#   theme(plot.subtitle = element_text(vjust = 1), plot.caption = element_text(vjust = 1), axis.line = element_line(linetype = "solid"), 
#         axis.ticks = element_line(colour = "gray0"), axis.title = element_text(size = 16, face = "bold"), 
#         axis.text = element_text(size = 13, face = "bold", colour = "black"), legend.text = element_text(size = 12, face = "bold"), 
#         legend.title = element_text(size = 12, face = "bold"), panel.background = element_rect(fill = NA), 
#         legend.key = element_rect(fill = NA), legend.background = element_rect(fill = NA)) +
#   labs(x = "Genes", y = "Cts") +
#   scale_fill_manual(name='Status', values=c('#FFFFFF', '#00BFC4', '#F8766D'))

results <- list()
for(i in 1:seq_along(length(qpcr.pvals))){
  results[[i]] <- sapply(qpcr.pvals[[i]][,5:10], function(x){t.test(x~Status, data=qpcr.pvals[[i]])$p.val})
}

pvals.p1 <- sapply(qpcr.pvals[[1]][,5:10], function(x){p.adjust(t.test(x~Status, data=qpcr.pvals[[1]])$p.val, method = 'bonferroni')})
pvals.p2 <- sapply(qpcr.pvals[[2]][,5:10], function(x){p.adjust(t.test(x~Status, data=qpcr.pvals[[2]])$p.val, method = 'bonferroni')})
pvals.p3 <- sapply(qpcr.pvals[[3]][,5:10], function(x){p.adjust(t.test(x~Status, data=qpcr.pvals[[3]])$p.val, method = 'bonferroni')})
pvals.p4 <- sapply(qpcr.pvals[[4]][,5:10], function(x){p.adjust(t.test(x~Status, data=qpcr.pvals[[4]])$p.val, method = 'bonferroni')})
write.xlsx(data.frame(rbind(pvals.p1,pvals.p2,pvals.p3, pvals.p4)), 'pvals.xlsx', sheetName = 'pvals') 
p.labels <- read.xlsx('pvals.xlsx', header=T, sheetIndex = 1)

p <- ggplot(qpcr.m, aes(x=variable, y=value, fill=Status)) + geom_boxplot(alpha=.5) + facet_grid(rows = vars(Plate)) + 
  stat_compare_means(method='t.test', comparisons = c('C', 'F')) + labs(x = "Genes", y = "Cts") +
  scale_fill_manual(name='Status', values=c('#FFFFFF', '#00BFC4', '#F8766D')) + 

  theme(plot.subtitle = element_text(vjust = 1), plot.caption = element_text(vjust = 1), axis.line = element_line(linetype = "solid"), 
        axis.ticks = element_line(colour = "gray0"), axis.title = element_text(size = 16, face = "bold"), 
        axis.text = element_text(size = 13, face = "bold", colour = "black"), legend.text = element_text(size = 12, face = "bold"), 
        legend.title = element_text(size = 12, face = "bold"), panel.background = element_rect(fill = NA), 
        legend.key = element_rect(fill = NA), legend.background = element_rect(fill = NA), strip.background = element_rect(fill='white'), 
        strip.text.y = element_text(size=20))

p + geom_text(data = p.labels, label=p.labels$label)


## geom_point() to indentify what to remove ####
data.clean <- qpcr.m
data.clean <- subset(data.clean, data.clean$Status!='Positive')
data.clean$Status <- droplevels(data.clean$Status)
library(ggrepel)
windows()
ggplot(data.clean, aes(x=variable, y=value, col=Status)) + geom_jitter() + facet_grid(rows = vars(Plate)) + labs(x = "Genes", y = "Cts") +
  scale_color_manual(name='Status', values=c('#00BFC4', '#F8766D')) + geom_vline(xintercept = 1.5) + geom_vline(xintercept = 2.5) +
  geom_vline(xintercept = 3.5) + geom_vline(xintercept = 4.5) + geom_vline(xintercept = 5.5) + 
  geom_text_repel(aes(label=ifelse(value>35, as.character(Sample), ''), hjust=0, vjust=0)) +
            
  
  theme(plot.subtitle = element_text(vjust = 1), plot.caption = element_text(vjust = 1), axis.line = element_line(linetype = "solid"), 
        axis.ticks = element_line(colour = "gray0"), axis.title = element_text(size = 16, face = "bold"), 
        axis.text = element_text(size = 13, face = "bold", colour = "black"), legend.text = element_text(size = 12, face = "bold"), 
        legend.title = element_text(size = 12, face = "bold"), panel.background = element_rect(fill = NA), 
        legend.key = element_rect(fill = NA), legend.background = element_rect(fill = NA), strip.background = element_rect(fill='white'), 
        strip.text.y = element_text(size=20))
  
## plate-normalized by control ####
setwd('qpcrs_384_good')
qpcr.normal <- read.xlsx2('qpcr_fullGOOD.xlsx', header=T, sheetIndex = 2)
sapply(qpcr.normal, class)
qpcr.normal[,5:10] <- sapply(qpcr.normal[,5:10], function(a)(as.numeric(as.character(a))))

library(reshape2)
qpcr.normal.m <- melt(qpcr.normal, id.vars = c('Sample', 'Plate', 'Status'), measure.vars = c('gadC', 'glnA', 'glsA1', 'glsA2', 'gadB1', 'gadB2'))
qpcr.normal.m$Plate <- as.factor(qpcr.normal.m$Plate)
levels(qpcr.normal.m$Plate) <- c('Plate 1', 'Plate 2', 'Plate 3', 'Plate 4')
windows()
ggplot(qpcr.normal.m, aes(x=variable, y=value, fill=Status)) +geom_boxplot(alpha=.5) + 
  facet_grid(rows = vars(Plate)) + 
 # stat_compare_means(method='t.test', comparisons = c('C', 'F')) + labs(x = "Genes", y = "Cts") +
  scale_fill_manual(name='Status', values=c('#00BFC4', '#F8766D')) + 
  # geom_segment(aes(x=1, xend=1.3, y=40, yend=40)) +
  
  
  theme(plot.subtitle = element_text(vjust = 1), plot.caption = element_text(vjust = 1), axis.line = element_line(linetype = "solid"), 
        axis.ticks = element_line(colour = "gray0"), axis.title = element_text(size = 16, face = "bold"), 
        axis.text = element_text(size = 13, face = "bold", colour = "black"), legend.text = element_text(size = 12, face = "bold"), 
        legend.title = element_text(size = 12, face = "bold"), panel.background = element_rect(fill = NA), 
        legend.key = element_rect(fill = NA), legend.background = element_rect(fill = NA), strip.background = element_rect(fill='white'), 
        strip.text.y = element_text(size=20))

qpcr.pvals.nor <- split(qpcr.normal, qpcr.normal$Plate)
pvals.nor.p1 <- sapply(qpcr.pvals.nor[[1]][,5:10], function(x){p.adjust(t.test(x~Status, data=qpcr.pvals.nor[[1]])$p.val, method = 'bonferroni')})
pvals.nor.p2 <- sapply(qpcr.pvals.nor[[2]][,5:10], function(x){p.adjust(t.test(x~Status, data=qpcr.pvals.nor[[2]])$p.val, method = 'bonferroni')})
pvals.nor.p3 <- sapply(qpcr.pvals.nor[[3]][,5:10], function(x){p.adjust(t.test(x~Status, data=qpcr.pvals.nor[[3]])$p.val, method = 'bonferroni')})
pvals.nor.p4 <- sapply(qpcr.pvals.nor[[4]][,5:10], function(x){p.adjust(t.test(x~Status, data=qpcr.pvals.nor[[4]])$p.val, method = 'bonferroni')})
write.xlsx(data.frame(rbind(pvals.nor.p1,pvals.nor.p2,pvals.nor.p3, pvals.nor.p4)), 'pvals.nor.xlsx', sheetName = 'pvals') 

##plot the full dataset
library(ggbeeswarm)
tiff('qpcr_plots.tiff', res=600, units='px', compression = 'lzw', height = 5000, width = 3500)
ggplot(qpcr.normal.m, aes(Status, value, fill=Status)) + geom_boxplot() + facet_grid(rows=vars(variable)) + 
  stat_compare_means(method='t.test', aes(label=paste0('p = ', round(..p.adj.., digits=3))), fontface=4, label.x = 2.1) + 
  # geom_beeswarm(alpha=.3) +   
  scale_fill_manual(name='Status', values=c('#00BFC4', '#F8766D')) + 
  theme(plot.subtitle = element_text(vjust = 1), plot.caption = element_text(vjust = 1), axis.line = element_line(linetype = "solid"), 
        axis.ticks = element_line(colour = "black"), axis.title = element_text(size = 19, face = "bold"),
        axis.text = element_text(size = 16, face = "bold", colour = "black"), legend.text = element_text(size = 17, face = "bold"), 
        legend.title = element_text(size = 17, face = "bold"), panel.background = element_rect(fill = NA), 
        legend.key = element_rect(fill = NA, size = 0.8), legend.background = element_rect(fill = NA), 
        strip.background = element_rect(fill='white'), strip.text.y = element_text(size=20)) +
  labs(x = NULL, y = "Ct (normalized)")
dev.off()

## deltaCT using average control ####
qpcr.delta <- read.xlsx2('qpcr_fullGOOD.xlsx', header=T, sheetIndex = 3)
sapply(qpcr.delta, class)
qpcr.delta[,5:10] <- sapply(qpcr.delta[,5:10], function(a)(as.numeric(as.character(a))))

qpcr.delta.m <- melt(qpcr.delta, id.vars = c('Sample', 'Plate', 'Status'), measure.vars = c('gadC', 'glnA', 'glsA1', 'glsA2', 'gadB1', 'gadB2'))
qpcr.delta.m$Plate <- as.factor(qpcr.delta.m$Plate)
levels(qpcr.delta.m$Plate) <- c('Plate 1', 'Plate 2', 'Plate 3', 'Plate 4')
windows()
ggplot(qpcr.delta.m, aes(x=variable, y=value, fill=Status)) + geom_boxplot(alpha=.5) + #geom_jitter(aes(col=Status)) + #geom_boxplot(alpha=.5) + 
  facet_grid(rows = vars(Plate)) + 
  stat_compare_means(method='t.test', comparisons = c('C', 'F')) + labs(x = "Genes", y = "Cts") +
  scale_fill_manual(name='Status', values=c('#FFFFFF', '#00BFC4', '#F8766D')) + ylim(0,15) +
  # geom_segment(aes(x=1, xend=1.3, y=40, yend=40)) +
  
  
  theme(plot.subtitle = element_text(vjust = 1), plot.caption = element_text(vjust = 1), axis.line = element_line(linetype = "solid"), 
        axis.ticks = element_line(colour = "gray0"), axis.title = element_text(size = 16, face = "bold"), 
        axis.text = element_text(size = 13, face = "bold", colour = "black"), legend.text = element_text(size = 12, face = "bold"), 
        legend.title = element_text(size = 12, face = "bold"), panel.background = element_rect(fill = NA), 
        legend.key = element_rect(fill = NA), legend.background = element_rect(fill = NA), strip.background = element_rect(fill='white'), 
        strip.text.y = element_text(size=20))


## PCA analysis ####
pca.qpcr <- qpcr.normal
pca.qpcr[,5:10] <- sapply(pca.qpcr[,5:10], NAfunc)
pca.qpcr <- prcomp(pca.qpcr[,5:10], scale=T, center=T)
plot(pca.qpcr)
plot(pca.qpcr, type='l')
library(factoextra)
fviz_eig(pca.qpcr)
fviz_pca_ind(pca.qpcr, geom=c('point'), addEllipses = F, col.ind = qpcr.normal$Plate, ellipse.level=.95) + theme_minimal()
fviz_pca_ind(pca.qpcr, geom=c('point'), addEllipses = F, col.ind = qpcr.normal$Status, ellipse.level=.95) + theme_minimal()
fviz_pca_var(pca.qpcr, geom='arrow')
fviz_pca_biplot(pca.qpcr, geom=c('point'), label = 'var', fill.ind = qpcr.normal$Plate, 
                col.ind = 'white', addEllipses = T, pointshape=27, pointsize=3) + scale_color_brewer(palette = 'Set1')
