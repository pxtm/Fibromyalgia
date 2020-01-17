## metabolomics fibro final dataset --------------------------------------------------------------------------------------
library(xlsx); library(ropls); library(ggplot2); library(ggrepel); library(ggthemes); library(ggpubr) 
library(corrplot); library(reshape2); library(igraph); library(DESeq2); library(mixOmics); library(venn)
library(dplyr); library(amritr);library(made4)

#setwd('C:/Users/mclos/Desktop')
# upload functions
NAfunc<-function(x){
  set.seed(1)
  x[x==0.000000]<-NA
  n<-sum(is.na(x))
  m<-min(x,na.rm=T)
  r<-m*0.1
  x[is.na(x)]<-r
  x
}


## upload and prepare metabolomics data ----
metabolites <- read.xlsx('final_metabolites_fibro.xlsx', header=T, sheetIndex = 2, row.names=1)
metabolites <- metabolites[!row.names(metabolites)%in%'FIBR033C',] ##remove bad sample FIBR033C
metabolites$Class <- as.factor(metabolites$Class)
remove <- row.names(subset(metabolites, metabolites$Class=='A'))
metabolites <- metabolites[!row.names(metabolites)%in%remove,] ##remove 'Related' samples
metabolites$Class <- droplevels(metabolites$Class)

## compute FoldChange, p-value and generate VolcanoPlot ------------------------------------------------------------------
fib <- subset(metabolites, metabolites$Class=='F')
ctl <- subset(metabolites, metabolites$Class!='F')
fib$Class <- NULL
ctl$Class <- NULL

mean.fib <- sapply(fib, mean)
mean.ctl <- sapply(ctl, mean)
FC <- mean.fib/mean.ctl

pval <- sapply(metabolites[,1:14], function(a)t.test(a~Class, data=metabolites)$p.val)

write.xlsx(data.frame(log2(FC), pval), 'IPA_finalfibro.xlsx', sheetName = '1')
volcano <- data.frame('FC'=log2(FC), 'pval'=-log10(pval), 'colour'=cut(log2(FC), c(-Inf, -.51, .37, Inf)))
names <- read.xlsx('final_metabolites_fibro.xlsx', header=T, sheetIndex = 4)
row.names(volcano) == names$Code
volcano <- data.frame('label' = names$Name, volcano)
tiff('final_metabolomics_names.tiff', res=600, height = 5500, width=5000)
ggplot(volcano, aes(FC, pval)) +  geom_point(aes(size = abs(FC)))+#geom_point(aes(colour=ifelse(volcano$pval>1.30103, colour, ''), size = abs(FC)))+
  xlab(bquote(bold(~log[2]~ ' Fold Change')))+ ylab(bquote(bold(~-log[10]~ ' p-values')))+
  geom_text_repel(aes(label=ifelse(pval>1.30103, as.character(label), '')), hjust=1.25, yjust=5)+
  geom_hline(yintercept=-log10(.05),linetype="dashed",col="red")+annotate("text",label="p-val=.05",y=1.48,x=-0.72,fontface=4)+
  geom_hline(yintercept=2,linetype="dashed",col="orange")+annotate("text",label="p-val=.01",fontface=4,y=2.17,x=-0.72)+
  geom_hline(yintercept=3,linetype="dashed",col="green")+annotate("text",label="p-val=.001",y=3.17,x=-0.71,fontface=4)+
  geom_vline(xintercept=0,linetype="dashed",col="black")+theme_base(base_size=17) + guides(colour=F, size= F) + xlim(-.75,.75) +
  theme(axis.title = element_text(size=17), axis.text=element_text(color = 'black', face = 'bold', size = 12), 
        axis.line = element_line(color='black', size=1), panel.border = element_rect(color='black', size=1.25, fill=NA)) +
  labs(caption='Positive FC means higher in Fibromyalgia patients')

dev.off()

## test data distribution ------------------------------------------------------------------------------------------------
for(i in seq_along(metabolites[,1:14])){
  multiplot(
    ggplot(metabolites, aes_string(x=names(metabolites)[i], fill='Class')) + geom_density(col=NA, alpha=.35),
    ggplot(log(metabolites[,1:14]), aes_string(x=names(metabolites)[i], fill=metabolites$Class)) + geom_density(col=NA, alpha=.35),
    cols=2)
}

library(ggpubr)
for(i in seq_along(metabolites[,1:14])){
  print(
    multiplot(
      ggqqplot(metabolites, x=names(metabolites)[i], title = names(metabolites)[i]),
      ggqqplot(log(metabolites[,1:14]), x=names(metabolites)[i], title = paste0('log ', names(metabolites)[i])),
      cols = 2))
}


## multivariate analysis of metabolomics ---------------------------------------------------------------------------------
browseVignettes('ropls')
## non supervised PCA ====================================================================================================
metabolites.pca <- opls(as.matrix(metabolites[,1:14]))
disease <- metabolites[,15]
plot(metabolites.pca, typeVc='x-score', parAsColFcVn = disease)
## supervised PLS-DA =====================================================================================================
metabolites.plsda <- opls(as.matrix(metabolites[,1:14]), disease, predI=3)
metabolites.plsda@loadingMN
loadings.pls <- as.data.frame(cbind('label'=rownames(getLoadingMN(metabolites.plsda)), getLoadingMN(metabolites.plsda)))
loadings.pls[,2:4] <- sapply(loadings.pls[,2:4], function(x)as.numeric(as.character(x)))
ggplot(loadings.pls, aes(p1, p2)) + geom_point() + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) +
  geom_text_repel(aes(label=label), hjust=1.05, yjust=0) + ylim(0.3, -0.3) + xlim(-0.4,0.4) + 
  theme(axis.line = element_line(size = 0.8, linetype = "solid"), axis.title = element_text(face = "bold"),
        axis.text = element_text(face = "bold"), panel.background = element_rect(fill = NA))
getSummaryDF(metabolites.plsda)
VIPs <- as.data.frame(sort(getVipVn(metabolites.plsda)))
names(VIPs) <- 'VIPscore'
VIPs$order <- factor(row.names(VIPs), levels=row.names(VIPs))

ggplot(VIPs, aes(x=VIPscore, y=order)) + geom_point(size=3) + 
  theme(axis.line = element_line(size = 0.3, linetype = "solid"), axis.text = element_text(size = 11, face = "bold", colour = "black"), 
        panel.background = element_rect(fill = NA)) +labs(x = "VIP score", y = "Variable") +
  theme(axis.title = element_text(size = 14,  face = "bold"), axis.ticks.y = element_blank()) + geom_segment(aes(yend=order), xend=0) +labs(y = NULL)


sapply(names(metabolites[1:14]), function(x){
  owl[grep(x, owl$KEGG), 3]
})


## correlations with cytokines, miRNAs and microbiome --------------------------------------------------------------------
## microbiome correlations ===============================================================================================
otus <- read.csv2('level6.csv', header=T, sep='\t', row.names = 1)
otus[] <- sapply(otus, function(x)as.numeric(as.character(x)))
otus <- as.data.frame(t(otus))
otus <- otus[, colSums(otus!=0)>36]
otus <- otus[row.names(otus)%in%row.names(metabolites),]
metabolites1 <- metabolites[row.names(metabolites)%in%row.names(otus),]
row.names(metabolites1) == row.names(otus)
otus <- otus[match(row.names(metabolites1), row.names(otus)),]
otus_deseq <- t(otus)
coldata <- data.frame('Status'=metabolites1$Class)
row.names(coldata) <- row.names(metabolites1)
ddstable <- DESeqDataSetFromMatrix(countData = otus_deseq, colData = coldata, design = ~Status)
gm_mean = function(x, na.rm=TRUE){  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
geoMeans = apply(counts(ddstable), 1, gm_mean)
statdds = estimateSizeFactors(ddstable, geoMeans = geoMeans)
statdds2 = DESeq(statdds, fitType="local")
otus_normal <- counts(statdds2, normalized = T)
otus_normal <- as.data.frame(t(otus_normal))

row.names(metabolites1) == row.names(otus_normal)
otus.metabolites <- data.frame(metabolites1, otus_normal)
otus.metabolites$Class <- NULL
cor.otusmetabolites <- cor(otus.metabolites)
res.otusmetabolites <- cor.mtest(otus.metabolites)
corrplot(cor.otusmetabolites[15:89,1:14], p.mat = res.otusmetabolites$p[15:89,1:14], insig='blank')

## cytokines correlations ================================================================================================
cytos <- readRDS('cytos.rds')
cytos$Class <- NULL
cytos <- cytos[row.names(cytos)%in%row.names(metabolites), ]
metabolites2 <- metabolites[row.names(metabolites)%in%row.names(cytos),]
row.names(metabolites2) == row.names(cytos)
cytos <- cytos[match(row.names(metabolites2), row.names(cytos)),]
cytos.metabolites <- data.frame(metabolites2, cytos)
cytos.metabolites$Class <- NULL
cytos.metabolites[] <- sapply(cytos.metabolites, NAfunc)
cor.cytosmetabs <- cor(cytos.metabolites)
res.cytosmetabs <- cor.mtest(cytos.metabolites)
corrplot(cor.cytosmetabs[15:66, 1:14], p.mat=res.cytosmetabs$p[15:66,1:14], insig = 'blank')

## miRNAs correlations ===================================================================================================
mirnas <- readRDS('micros.rds')
row.names(mirnas) <- mirnas$Sample
mirnas$Group <- NULL; mirnas$Sample <- NULL
mirnas <- mirnas[row.names(mirnas)%in%row.names(metabolites),]
metabolites3 <- metabolites[row.names(metabolites)%in%row.names(mirnas), ]
row.names(metabolites3) == row.names(mirnas)
mirnas <- mirnas[match(row.names(metabolites3), row.names(mirnas)),]
mirnas.metabolites <- data.frame(metabolites3, mirnas)
mirnas.metabolites$Class <- NULL
cor.microsmetabs <- cor(mirnas.metabolites)
res.microsmetabs <- cor.mtest(mirnas.metabolites)
corrplot(cor.microsmetabs[15:63, 1:14], p.mat=res.microsmetabs$p[15:63,1:14], insig = 'blank')

## metadata correlations (pain, drugs, lifestyle) ========================================================================
metadata <- read.xlsx('metadata_full.xlsx', header=T, row.names=1, sheetIndex = 1)
metadata <- metadata[row.names(metadata)%in%row.names(metabolites),]
pain <- read.csv2('pain_indices.csv', header=T, row.names=1)
pain <- pain[row.names(pain)%in%row.names(metadata),]
pain <- pain[match(row.names(metadata), row.names(pain)),]
metabolites4 <- metabolites[row.names(metabolites)%in%row.names(metadata),]
metabolites4 <- metabolites4[match(row.names(pain), row.names(metabolites4)),]
#metabolomics.vs.metadata <- data.frame(metabolites4, metadata, pain)
metabolomics.pain <- data.frame(metabolites4, pain)
keep.pain <- row.names(subset(metadata, metadata$Diagnostico=='Fibromyalgia'))
metabolomics.pain <- metabolomics.pain[row.names(metabolomics.pain)%in%keep.pain,]
metabolomics.pain$Class <- NULL
names.metabolites <- read.xlsx('final_metabolites_fibro.xlsx', header=T, sheetIndex = 4)
names(metabolomics.pain)[1:14] <- as.character(names.metabolites$Name)
metabolomics.pain$SS1[is.na(metabolomics.pain$SS1)] <- median(metabolomics.pain$SS1, na.rm = T)
library(corrplot)
pain.cor <- cor(metabolomics.pain)
pain.res <- cor.mtest(metabolomics.pain)
corrplot(pain.cor, sig.level = pain.res$p, insig = 'blank')
tiff('metabolitesvspain.tiff', width = 3500, height = 5500, res=600, units = 'px')
corrplot(pain.cor[1:14, 15:22], sig.level=pain.res$p[1:14, 15:22], insig='blank', tl.col='black')
dev.off()

## mixOmics analysis - multivariate rCCA correlation analysis between datasets -------------------------------------------
## mixOmics package instructions followed ##
keep <- row.names(cytos)
## combine 4 datasets with identified metabolites ========================================================================
  ## prepare data in matrix format #######################################################################################
otus1 <- otus_normal[row.names(otus_normal)%in%keep,]
otus1m <- as.matrix(otus1); colnames(otus1m) <- names(otus1); rownames(otus1m) <- row.names(otus1)
metabs <- metabolites[row.names(metabolites)%in%keep,]
metabsm <- as.matrix(metabs); colnames(metabsm) <- names(metabs); rownames(metabsm) <- row.names(metabs) ##before doing this, remove Class column
cytosm <- as.matrix(cytos); colnames(cytosm) <- names(cytos); rownames(cytosm) <- row.names(cytos)
mirnasm <- as.matrix(mirnas); colnames(mirnasm) <- names(mirnas); rownames(mirnasm) <- row.names(mirnas)
  ## outcome
Y <- metabs$Class 
metabs$Class <- NULL
 ## prepare data list
data <- list(Microbiome = otus1m,
             Metabolites = metabsm,
             Cytokines = cytosm,
             miRNAs = mirnasm)
lapply(data, dim)
lapply(data, class)

  ## parameter choice ####################################################################################################
design <- matrix(0.1, ncol=length(data), nrow = length(data), dimnames= list(names(data), names(data)))
diag(design)=0
design

sgccda.res <- block.splsda(data, Y, ncomp=10, design = design, near.zero.var = T)

set.seed(2810) # reproducibility
set.seed(28101)
t1 <- proc.time()
perf.diablo <- perf(sgccda.res, validation = 'Mfold', folds=10, nrepeat = 10, auc = T, cpus = 4)
t2 <- proc.time()
running.time <- t2-t1; running.time

perf.diablo  # lists the different outputs
plot(perf.diablo) 
perf.diablo$choice.ncomp$WeightedVote # number of components selection
ncomp = perf.diablo$choice.ncomp$WeightedVote["Overall.ER", "centroids.dist"]

set.seed(0110) # reproducibility
test.keepX <- list(Microbiome = sample(1:75, 10),
                   Metabolites = sample(1:14, 5),
                   Cytokines = sample(1:52, 5),
                   miRNAs = sample(1:49, 5))
t1 <- proc.time() ## takes looooooooooooooooong time! selects 17 variables per component. Try reducing!
tune.data = tune.block.splsda(X = data, Y = Y, ncomp = ncomp, 
                              test.keepX = test.keepX, design = design,
                              validation = 'Mfold', folds = 10, nrepeat = 10, dist = "centroids.dist", cpus=4)#, near.zero.var = T)
t2 = proc.time()
running_time = t2 - t1; running_time

list.keepX = tune.data$choice.keepX # which variables to keep?
list.keepX
save(tune.data,list.keepX, running_time, file = '.RData')
saveRDS(tune.data, 'fulldataidmets.rds')
  ## final model generation ##############################################################################################
sgccda.res = block.splsda(X = data, Y = Y, ncomp = ncomp, #generate the new model
                          keepX = list.keepX, design = design)
  ## loadings of each dataset for the first component
sgccda.res
sgccda.res$design
selectVar(sgccda.res, block = 'Microbiome', comp = 1) 
selectVar(sgccda.res, block = 'Metabolites', comp = 1)
selectVar(sgccda.res, block = 'Cytokines', comp = 1)
selectVar(sgccda.res, block = 'miRNAs', comp = 1)

diablopanel <- mapply(function(x,y){
  c(x,y)},
  x = lapply(selectVar(sgccda.res, comp=1)[-5], function(i)unlist(i[[1]])),
  y = lapply(selectVar(sgccda.res, comp=2)[-5], function(i)unlist(i[[1]])))
sapply(diablopanel, length)
diablopanel <- lapply(diablopanel, unique)

 ## SPLSDA plots of the different pairwise combination of datasets + correlation of one whole dataset vs another one 
tiff('diablo_scoreplots.tiff', res=600, height = 4500, width=4500)
plotDiablo(sgccda.res, legend = T, cex=20) + scale_shape_manual(values = c(1,2))
dev.off()
plotDiablo(sgccda.res, legend = T, cex=20, ncomp=2)
plotDiablo(sgccda.res, legend = T, cex=20, ncomp=3)

vardat <- do.call(rbind, sgccda.res$variates[1:length(data)]) %>%
  as.data.frame %>% mutate(subj=row.names(.))
vardat$Dataset <- rep(names(data), e=length(Y))
vardat$Class <- Y
vardat <- vardat %>% group_by(Class, subj) %>% 
  summarise_all(funs(mean)) %>% 
  mutate(Dataset = "Consensus") %>% 
  dplyr::select(`comp 1`, `comp 2`, `comp 3`,  subj, Dataset, Class) %>% 
  as.data.frame() %>% 
  rbind(., vardat) %>% 
  mutate(Dataset = factor(Dataset, c("Consensus", "Microbiome", "Metabolites", "Cytokines", "miRNAs"))) 

consensus.plot <- filter(vardat, Dataset == "Consensus") %>% 
  ggplot(aes(x = `comp 1`, y = `comp 2`, color = Class, shape=Class)) + 
  geom_point(size = 5) +
  stat_ellipse(data = filter(vardat, Dataset == "Consensus"), size = 1) +
  xlab("Component 1") + ylab("Component 2") + 
  theme(strip.text.x = element_text(size=26, face = "bold")) + 
  scale_color_manual(values=color.mixo(1:4)) + scale_shape_manual(values=c(1,2))

tiff('diablo_consensusplot.tiff', res=600, height = 3000, width = 3000)
consensus.plot + 
  theme(axis.ticks = element_line(colour = "black", size=.6),axis.title = elemen?plott_text(size = 13, face = "bold"), 
        axis.text = element_text(face = "bold", colour = "black", size =12), 
        plot.title = element_text(face = "bold", size = 15, hjust=.5), legend.text = element_text(size = 12, face = "bold"), 
        legend.title = element_text(size = 12, face = "bold"), legend.key = element_rect(fill = NA),
        legend.background = element_rect(fill = NA), axis.line = element_line(size = 0.5, linetype = "solid"),
        panel.background = element_rect(fill = NA)) +
  labs(title = "Consensus plot") + guides(color=guide_legend(override.aes=list(linetype=0)))
dev.off()
 ## Individual PLS-DA plots for each omics, to identify the contributions of each one to the whole separation
tiff('diablo_individualcontributions.tiff', res=600, height = 4500, width=4500)
plotIndiv(sgccda.res, ind.names = F, legend = TRUE,  style="ggplot2", ellipse = T)
dev.off()
 ## ARROW plot: start of the arrow indicates the centroid between all datasets for a given sample and the tip of the arrow
 ## the location of the same sample in each block. This plot highlights the agreement between all datasets at the sample level
 ## when modeled with DIABLO
plotArrow(sgccda.res, ind.names = FALSE, legend = TRUE, title = 'DIABLO')
 ## Variable plot: correlation between variables and components. Most important variables should be closer to the large circle.
plotVar(sgccda.res, var.names = FALSE, style = 'graphics', legend = TRUE, 
        pch = c(15, 16, 17, 12), cex = c(2,2,2,2), col = c('darkgreen', 'blue4', 'chocolate1', 'darkorchid2'))
plotVar(sgccda.res, var.names = FALSE, style = 'graphics', legend = TRUE, 
        pch = c(15, 16, 17, 12), cex = c(2,2,2,2), col = c('darkgreen', 'blue4', 'chocolate1', 'darkorchid2'),
        comp.select = 1)
plotVar(sgccda.res, var.names = FALSE, style = 'graphics', legend = TRUE, 
        pch = c(15, 16, 17, 12), cex = c(2,2,2,2), col = c('darkgreen', 'blue4', 'chocolate1', 'darkorchid2'),
        comp.select = 2)
plotVar(sgccda.res, var.names = FALSE, style = 'ggplot2', legend = TRUE, 
        pch = c(15, 16, 17, 12), cex = c(3.5,3.5,3.5,3.5), col = c('darkgreen', 'blue4', 'chocolate1', 'darkorchid2'))
plotVar(sgccda.res, var.names = F, style = '3d', legend = T, blocks = c('Microbiome', 'Metabolites', 'Cytokines', 'miRNAs'))
 ## Circos plot also represents the correlations BETWEEN differents variables.
#tiff('circosfirstcomp.tiff', res=600, height = 4500, width = 4500)
pdf('diablo_circosfirstcomp.pdf', height = 12, width = 12)
circosPlot(sgccda.res, cutoff=.5, line = T, 
           color.blocks= c('darkgreen', 'blue4', 'chocolate1', 'darkorchid2'),
           color.cor = c("cadetblue4","coral4"), showIntraLinks = F, size.variables = 0.6, comp = 1)
dev.off()
 ## Relevant network: unites the correlated variables, indicating which kind of correlation is established between each
 ## one of them. Can be exported to Cytoscape.
my.network <- network(sgccda.res, blocks = c(1,2,3,4), color.node = c('darkgreen', 'blue4', 'chocolate1', 'darkorchid2'), cutoff = .4)
write.graph(my.network$gR, file='diablo_networknormal.gml', format='gml')
 ## Plot loadings of each block per each component.
plotLoadings(sgccda.res, contrib = 'max', method = 'median', size.name = .9)#, xlim=c(-.6,.6))
plotLoadings(sgccda.res, contrib = 'max', method = 'median', size.name = .9, comp=2)
plotLoadings(sgccda.res, contrib = 'max', method = 'median', size.name = .9, comp=3)
 ## HEatmap of the combination of the distinct blocks per component.
cimDiablo(sgccda.res, color.blocks = c('darkgreen', 'blue4', 'chocolate1', 'darkorchid2'),
          margins=c(8,20), legend.position = "right", transpose = F, col.names=F) #+ 

## performance of the model computation ##################################################################################
set.seed(129)
t1 <- proc.time()
perf.diablo <- perf(sgccda.res, validation = 'Mfold', folds = 10, nrepeat = 10, dist = 'centroids.dist')
t2 <- proc.time()
running_time <- t2-t1; running_time
perf.diablo
perf.diablo$MajorityVote.error.rate
perf.diablo$WeightedVote.error.rate
 ## Extract the ROC curve for each OMICS block and plot the 4 of them together
micro.plot <- (auroc(sgccda.res)$graph.Microbiome)$comp1$data
metab.plot <- (auroc(sgccda.res, roc.block=2)$graph.Metabolites)$comp1$data
cytos.plot <- (auroc(sgccda.res, roc.block=3)$graph.Cytokines)$comp1$data
mirna.plot <- (auroc(sgccda.res, roc.block=4)$graph.miRNAs)$comp1$data
aucs <- melt(list(Microbiome=micro.plot, Metabolites=metab.plot, Cytokines=cytos.plot, miRNAs=mirna.plot), 
             id=c('Sensitivity', 'Specificity'))
aucs$L1 <- as.factor(aucs$L1)
aucs$L1 <- ordered(aucs$L1, levels=c('Microbiome', 'Metabolites', 'Cytokines', 'miRNAs'))
tiff('diablo_aucs.tiff', res=600, height = 4500, width = 4100)
ggplot(aucs, aes(Specificity, Sensitivity, color=L1)) + geom_line(size=1) +#), linetype='dashed') + 
  scale_color_manual('Omics type', values=c('Cytokines' = 'chocolate1', 'Metabolites'='blue', 
                                            'Microbiome' = 'chartreuse4', 'miRNAs'='darkorchid2' )) +
  geom_abline(intercept=0, linetype="dashed") + 
  theme(axis.line = element_line(size = 0.5, linetype = "solid"), axis.title = element_text(size = 13, face = "bold"), 
        axis.text = element_text(size = 12, face = "bold", colour = "black"), legend.title = element_text(size = 12,face = "bold"), 
        panel.background = element_rect(fill = NA), legend.key = element_rect(fill = "white"), 
        legend.background = element_rect(fill = "white"), legend.text = element_text(face = "bold", size = 11), legend.position = c(0.9,0.25)) +
  labs(x = "100 - Specificity (%)", y = "Sensitivity (%)") +
  annotate('text', label='AUC = 0.771', colour = 'chocolate1', x=95, y=3, fontface=2) + 
  annotate('text', label='AUC = 0.810', colour='chartreuse4', x=95, y=9, fontface=2) + 
  annotate('text', label='AUC = 0.747', colour='blue', x=95, y=6, fontface=2) + 
  annotate('text', label='AUC = 0.727', colour='black', x=95, y=0, fontface=2) + theme(legend.text = element_text(face = "bold"))
dev.off()

 ## ropls for whole data #################################################################################################
all.data <- data.frame(cbind(otus1m, metabsm, cytosm, mirnasm))
plsda.all <- opls(all.data, Y, predI=3, permI=20)
scores <- as.data.frame(plsda.all@scoreMN)
scores$class <- Y
 ## scores plot
tiff('diablo_scoresropls.tiff', height = 4300, width = 4300, res=600)
ggplot(scores, aes(p1, p2, color=class)) + geom_point() + stat_ellipse(geom = 'polygon', aes(fill=class), alpha=0.3) + 
  labs(x = "t1 (4%)", y = "t2 (3%)") + scale_fill_manual(values=c('cadetblue1', 'coral2')) + 
  scale_color_manual(values=c('blue2', 'red')) + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + 
  labs(title = "Scores (PLS-DA)") + xlim(-7.5,7.5) + ylim(-7.5,7.5) +
  theme(axis.line = element_line(size = 0.6, linetype = "solid"), axis.ticks = element_line(colour = "gray0"), 
       axis.title = element_text(size = 13, face = "bold"), axis.text = element_text(face = "bold", colour = "gray0"), 
       panel.background = element_rect(fill = NA), plot.background = element_rect(colour = NA, linetype = "solid"), legend.text = element_text(face = "bold"), 
       legend.title = element_blank(), plot.title = element_text(hjust = .5), 
       panel.grid.major = element_line(colour = "gray89"), panel.grid.minor = element_line(linetype = "blank"), 
       panel.border=element_rect(colour = 'black', fill=NA, size = 3), legend.position = c(0.05,0.94))
dev.off() 
## loadings plot -> we only display the names of those variables that are over ABS(#comp)>.15
loadings.plsda <- as.data.frame(cbind('label'=rownames(getLoadingMN(plsda.all)), getLoadingMN(plsda.all)))
loadings.plsda[,2:4] <- sapply(loadings.plsda[,2:4], function(x)as.numeric(as.character(x)))
tiff('diablo_loadings_plsdaALL.tiff', res=600, height = 6500, width = 6500)
ggplot(loadings.plsda, aes(p1, p2)) + geom_point(aes(colour=ifelse(abs(loadings.plsda$p1)>.15|abs(loadings.plsda$p2)>.15, 'red', 'black'))) + 
  geom_hline(yintercept = 0) + geom_vline(xintercept = 0) +
  geom_text_repel(data=subset(loadings.plsda, abs(p1)>0.15 | abs(p2)>0.15), aes(label=label), cex=3, hjust=1.05, yjust=0) + ylim(0.3, -0.3) + xlim(-0.4,0.4) + 
  theme(axis.line = element_line(size = 0.8, linetype = "solid"), axis.title = element_text(face = "bold"),
        axis.text = element_text(face = "bold"), panel.background = element_rect(fill = NA)) + scale_color_manual(values=c('black', 'red')) +
  guides(colour=F)
dev.off()
getSummaryDF(plsda.all)
 ## VIPs plot -> we generated two distinct plots, as there were 190 variables, plotting all of them in one and the
 ## 50 most contributing variables to allow and easier inspection.
VIPs <- as.data.frame(sort(getVipVn(plsda.all)))
names(VIPs) <- 'VIPscore'
VIPs$order <- factor(row.names(VIPs), levels=row.names(VIPs))
ggplot(VIPs, aes(x=VIPscore, y=order)) + geom_point(size=3) + 
  theme(axis.line = element_line(size = 0.3, linetype = "solid"), axis.text = element_text(size = 11, face = "bold", colour = "black"), 
        panel.background = element_rect(fill = NA)) +labs(x = "VIP score", y = "Variable") +
  theme(axis.title = element_text(size = 14,  face = "bold"), axis.ticks.y = element_blank()) + geom_segment(aes(yend=order), xend=0) +labs(y = NULL)
VIPs.subset <- tail(VIPs, 50) ## plot only the 50 most contributing variables
tiff('diablo_VIPs.tiff', height = 6000, width = 4500, res=600)
ggplot(VIPs.subset, aes(x=VIPscore, y=order)) + geom_point(size=3) + 
  theme(axis.line = element_line(size = 0.3, linetype = "solid"), axis.text = element_text(size = 11, face = "bold", colour = "black"), 
        panel.background = element_rect(fill = NA)) +labs(x = "VIP score", y = "Variable") +
  theme(axis.title = element_text(size = 14,  face = "bold"), axis.ticks.y = element_blank()) + geom_segment(aes(yend=order), xend=0) +labs(y = NULL)
dev.off()

#####

## combine 4 datasets with full metabolomics (no identified metabolites) =================================================
full.metabolomics <- readRDS('fullmetabolomics.rds')
full.metabolomics$Status <- NULL
full.metabolomics <- full.metabolomics[row.names(full.metabolomics)%in%keep,]
full.metabolomics <- full.metabolomics[match(rownames(otus1m), row.names(full.metabolomics)),]
full.metabolomicsm <- as.matrix(full.metabolomics); colnames(full.metabolomicsm) <- names(full.metabolomics); rownames(full.metabolomicsm) <- row.names(full.metabolomics)

data2 <- list(Microbiome = otus1m,
             Metabolites = full.metabolomicsm,
             Cytokines = cytosm,
             miRNAs = mirnasm)
lapply(data2, dim)

## parameter choice ####################################################################################################
design2 <- matrix(0.1, ncol=length(data2), nrow = length(data2), dimnames= list(names(data2), names(data2)))
diag(design2)=0
design2

sgccda.res <- block.splsda(data2, Y, ncomp=25, design = design2, near.zero.var = T)

set.seed(1992) # reproducibility
t1 <- proc.time()
perf.diablo <- perf(sgccda.res, validation = 'Mfold', folds=10, nrepeat = 10, auc = T, cpus = 4)
t2 <- proc.time()
running.time <- t2-t1; running.time

perf.diablo  # lists the different outputs
plot(perf.diablo) 
perf.diablo$choice.ncomp$WeightedVote # number of components selection
ncomp = perf.diablo$choice.ncomp$WeightedVote["Overall.ER", "max.dist"]

set.seed(1001) # reproducibility
test.keepX <- list(Microbiome = sample(1:75, 10),
                   Metabolites = sample(1:1070, 20),
                   Cytokines = sample(1:52, 5),
                   miRNAs = sample(1:49, 5))
t1 <- proc.time() ## takes looooooooooooooooong time! selects 17 variables per component. Try reducing!
tune.data = tune.block.splsda(X = data2, Y = Y, ncomp = ncomp, 
                              test.keepX = test.keepX, design = design2,
                              validation = 'Mfold', folds = 10, nrepeat = 10, dist = "max.dist", cpus=4)#, near.zero.var = T)
t2 = proc.time()
running_time = t2 - t1; running_time

list.keepX = tune.data$choice.keepX # which variables to keep?
list.keepX
save(tune.data,list.keepX, running_time, file = '.RData')
saveRDS(tune.data, 'fullmetabolomicstuned.RDS')
tune.data <- readRDS('fullmetabolomicstuned.RDS')
## final model generation ##############################################################################################
sgccda.res = block.splsda(X = data2, Y = Y, ncomp = ncomp, #generate the new model
                          keepX = list.keepX, design = design2)
## loadings of each dataset for the first component
sgccda.res
sgccda.res$design
selectVar(sgccda.res, block = 'Microbiome', comp = 1) 
metabolites.used <- c(selectVar(sgccda.res, block = 'Metabolites', comp = 1)$Metabolites$name, 
                      selectVar(sgccda.res, block = 'Metabolites', comp = 2)$Metabolites$name)

selectVar(sgccda.res, block = 'Cytokines', comp = 1)
selectVar(sgccda.res, block = 'miRNAs', comp = 1)

## number of variables, overlapping and variable analysis of the model
diablopanel <- mapply(function(x,y){
  c(x,y)},
  x = lapply(selectVar(sgccda.res, comp=1)[-5], function(i)unlist(i[[1]])),
  y = lapply(selectVar(sgccda.res, comp=2)[-5], function(i)unlist(i[[1]])))
sapply(diablopanel, length)
diablopanel <- lapply(diablopanel, unique)


## SPLSDA plots of the different pairwise combination of datasets + correlation of one whole dataset vs another one 
tiff('fullmetabolomics_scoreplots.tiff', res=600, height = 4500, width=4500)
plotDiablo(sgccda.res, legend = T, style='ggplot2')
dev.off()
plotDiablo(sgccda.res, legend = T, cex=20, ncomp=2)

vardat <- do.call(rbind, sgccda.res$variates[1:length(data2)]) %>%
  as.data.frame %>% mutate(subj=row.names(.))
vardat$Dataset <- rep(names(data2), e=length(Y))
vardat$Class <- Y
vardat <- vardat %>% group_by(Class, subj) %>% 
  summarise_all(funs(mean)) %>% 
  mutate(Dataset = "Consensus") %>% 
  dplyr::select(`comp 1`, `comp 2`,  subj, Dataset, Class) %>% 
  as.data.frame() %>% 
  rbind(., vardat) %>% 
  mutate(Dataset = factor(Dataset, c("Consensus", "Microbiome", "Metabolites", "Cytokines", "miRNAs"))) 

consensus.plot <- filter(vardat, Dataset == "Consensus") %>% 
  ggplot(aes(x = `comp 1`, y = `comp 2`, color = Class, shape=Class)) + 
  geom_point(size = 5) +
  stat_ellipse(data = filter(vardat, Dataset == "Consensus"), size = 1) +
  xlab("Component 1") + ylab("Component 2") + 
  theme(strip.text.x = element_text(size=26, face = "bold")) + 
  scale_color_manual(values=color.mixo(1:4)) + scale_shape_manual(values=c(1,2))

tiff('fullmetabolomics_consensusplot.tiff', res=600, height = 3000, width = 3000)
consensus.plot + 
  theme(axis.ticks = element_line(colour = "black", size=.6),axis.title = element_text(size = 13, face = "bold"), 
        axis.text = element_text(face = "bold", colour = "black", size =12), 
        plot.title = element_text(face = "bold", size = 15, hjust=.5), legend.text = element_text(size = 12, face = "bold"), 
        legend.title = element_text(size = 12, face = "bold"), legend.key = element_rect(fill = NA),
        legend.background = element_rect(fill = NA), axis.line = element_line(size = 0.5, linetype = "solid"),
        panel.background = element_rect(fill = NA)) +
  labs(title = "Consensus plot") + guides(color=guide_legend(override.aes=list(linetype=0)))
dev.off()
## Individual PLS-DA plots for each omics, to identify the contributions of each one to the whole separation

tiff('fullmetabolomics_individualcontributions.tiff', res=600, height = 4500, width=4500)
filter(vardat, Dataset != "Consensus") %>% 
  ggplot(aes(x = `comp 1`, y = `comp 2`, color = Class, shape=Class)) + 
  geom_point(size = 3) +
  facet_wrap(~Dataset, scales = "free", ncol = 2) + 
  stat_ellipse(data = filter(vardat, Dataset != "Consensus"), size = 1) +
  xlab("Component 1") + ylab("Component 2") + 
  scale_color_manual(values=color.mixo(1:4)) + scale_shape_manual(values=c(1,2)) + 
  guides(color=guide_legend(override.aes=list(linetype=0))) +
  theme(strip.text.x = element_text(size=15, face = "bold"), axis.line = element_line(size = 0.3,linetype = "solid"), 
        axis.ticks = element_line(colour = "black"), panel.grid.major = element_line(colour = "gray94"), 
        panel.grid.minor = element_line(colour = "gray94"), axis.text = element_text(size = 12, face = "bold",colour = "black"), 
        panel.background = element_rect(fill = NA), panel.border=element_rect(colour = 'black', fill=NA, size = 1.8),
        strip.background = element_blank(), legend.text = element_text(size = 12, face = "bold"), 
        legend.title = element_text(size = 13, face = "bold"), legend.key = element_rect(fill = NA), 
        legend.background = element_rect(fill = NA), axis.title = element_text(size = 15, face = "bold"), 
        plot.title = element_text(face = "bold"))
dev.off()
## ARROW plot: start of the arrow indicates the centroid between all datasets for a given sample and the tip of the arrow
## the location of the same sample in each block. This plot highlights the agreement between all datasets at the sample level
## when modeled with DIABLO
plotArrow(sgccda.res, ind.names = FALSE, legend = TRUE, title = 'DIABLO')
## Variable plot: correlation between variables and components. Most important variables should be closer to the large circle.
plotVar(sgccda.res, var.names = FALSE, style = 'graphics', legend = TRUE, 
        pch = c(15, 16, 17, 12), cex = c(2,2,2,2), col = c('darkgreen', 'blue4', 'chocolate1', 'darkorchid2'))
plotVar(sgccda.res, var.names = FALSE, style = 'graphics', legend = TRUE, 
        pch = c(15, 16, 17, 12), cex = c(2,2,2,2), col = c('darkgreen', 'blue4', 'chocolate1', 'darkorchid2'),
        comp.select = 1)
plotVar(sgccda.res, var.names = FALSE, style = 'graphics', legend = TRUE, 
        pch = c(15, 16, 17, 12), cex = c(2,2,2,2), col = c('darkgreen', 'blue4', 'chocolate1', 'darkorchid2'),
        comp.select = 2)
plotVar(sgccda.res, var.names = FALSE, style = 'ggplot2', legend = TRUE, 
        pch = c(15, 16, 17, 12), cex = c(3.5,3.5,3.5,3.5), col = c('darkgreen', 'blue4', 'chocolate1', 'darkorchid2'))
plotVar(sgccda.res, var.names = F, style = 'ggplot2', legend = T, blocks = c('Microbiome', 'Metabolites'))
## Circos plot also represents the correlations BETWEEN differents variables.
#tiff('circosfirstcomp.tiff', res=600, height = 4500, width = 4500)
pdf('fullmetabolomics_circosfirstcomp.pdf', height = 12, width = 12)
circosPlot(sgccda.res, cutoff=.7, line = T, 
           color.blocks= c('darkgreen', 'blue4', 'chocolate1', 'darkorchid2'),
           color.cor = c("cadetblue4","coral4"), showIntraLinks = F, size.variables = 0.6, comp=1)
dev.off()
## Relevant network: unites the correlated variables, indicating which kind of correlation is established between each
## one of them. Can be exported to Cytoscape.
my.network <- network(sgccda.res, blocks = c(1,2,3,4), color.node = c('darkgreen', 'blue4', 'chocolate1', 'darkorchid2'), cutoff = .5)
write.graph(my.network$gR, file='fullmetabolomics_networkcutoff5.gml', format='gml')
## Plot loadings of each block per each component.
pdf('fullmetabolomics_loadings.pdf', width = 15, height = 10)
plotLoadings(sgccda.res, contrib = 'max', method = 'median', size.name = .9)#, xlim=c(-.6,.6))
dev.off()
plotLoadings(sgccda.res, contrib = 'max', method = 'median', size.name = .9, comp=2)
## HEatmap of the combination of the distinct blocks per component.
cimDiablo(sgccda.res, color.blocks = c('darkgreen', 'blue4', 'chocolate1', 'darkorchid2'),
          margins=c(8,20), legend.position = "right", transpose = F, col.names = F) #+ 

## performance of the model computation ##################################################################################
set.seed(131)
t1 <- proc.time()
perf.diablo <- perf(sgccda.res, validation = 'Mfold', folds = 10, nrepeat = 10, dist = 'centroids.dist')
t2 <- proc.time()
running_time <- t2-t1; running_time
plot(perf.diablo)
perf.diablo$MajorityVote.error.rate
perf.diablo$WeightedVote.error.rate
## Extract the ROC curve for each OMICS block and plot the 4 of them together
micro.plot <- (auroc(sgccda.res)$graph.Microbiome)$comp1$data
metab.plot <- (auroc(sgccda.res, roc.block=2)$graph.Metabolites)$comp1$data
cytos.plot <- (auroc(sgccda.res, roc.block=3)$graph.Cytokines)$comp1$data
mirna.plot <- (auroc(sgccda.res, roc.block=4)$graph.miRNAs)$comp1$data
aucs <- melt(list(Microbiome=micro.plot, Metabolites=metab.plot, Cytokines=cytos.plot, miRNAs=mirna.plot), 
             id=c('Sensitivity', 'Specificity'))
aucs$L1 <- as.factor(aucs$L1)
aucs$L1 <- ordered(aucs$L1, levels=c('Microbiome', 'Metabolites', 'Cytokines', 'miRNAs'))
tiff('fullmetabolomics_aucs.tiff', res=600, height = 4500, width = 3500)
ggplot(aucs, aes(Specificity, Sensitivity, color=L1)) + geom_line(size=1) +#), linetype='dashed') + 
  scale_color_manual('Omics type', values=c('Cytokines' = 'chocolate1', 'Metabolites'='blue', 
                                            'Microbiome' = 'chartreuse4', 'miRNAs'='darkorchid2' )) +
  geom_abline(intercept=0, linetype="dashed") + 
  theme(axis.line = element_line(size = 0.5, linetype = "solid"), axis.title = element_text(size = 13, face = "bold"), 
        axis.text = element_text(size = 12, face = "bold", colour = "black"), legend.title = element_text(size = 12,face = "bold"), 
        panel.background = element_rect(fill = NA), legend.key = element_rect(fill = "white"), 
        legend.background = element_rect(fill = "white"), legend.text = element_text(face = "bold", size = 11), legend.position = c(0.87,0.25)) +
  labs(x = "100 - Specificity (%)", y = "Sensitivity (%)") +
  annotate('text', label='AUC = 0.797', colour = 'chocolate1', x=95, y=3, fontface=2) + 
  annotate('text', label='AUC = 0.833', colour='chartreuse4', x=95, y=9, fontface=2) + 
  annotate('text', label='AUC = 0.804', colour='blue', x=95, y=6, fontface=2) + 
  annotate('text', label='AUC = 0.712', colour='darkorchid2', x=95, y=0, fontface=2)
dev.off()

## ropls for whole data #################################################################################################
all.data.full <- data.frame(cbind(otus1m, full.metabolomicsm, cytosm, mirnasm))
plsda.all.full <- opls(all.data.full, Y)#, predI=3, permI=20)
plot(plsda.all.full, parAsColFcVn=Y)
scores <- as.data.frame(plsda.all.full@scoreMN)
scores$class <- Y
## scores plot
tiff('fullmetabolomics_scoresropls.tiff', height = 4300, width = 4300, res=600)
ggplot(scores, aes(p1, p2, color=class)) + geom_point() + stat_ellipse(geom = 'polygon', aes(fill=class), alpha=0.3) + 
  labs(x = "t1 (6%)", y = "t2 (4%)") + scale_fill_manual(values=c('cadetblue1', 'coral2')) + 
  scale_color_manual(values=c('blue2', 'red')) + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + 
  labs(title = "Scores (PLS-DA)") + 
  theme(axis.line = element_line(size = 0.6, linetype = "solid"), axis.ticks = element_line(colour = "gray0"), 
        axis.title = element_text(size = 13, face = "bold"), axis.text = element_text(face = "bold", colour = "gray0"), 
        panel.background = element_rect(fill = NA), plot.background = element_rect(colour = NA, linetype = "solid"), legend.text = element_text(face = "bold"), 
        legend.title = element_blank(), plot.title = element_text(hjust = .5), 
        panel.grid.major = element_line(colour = "gray89"), panel.grid.minor = element_line(linetype = "blank"), 
        panel.border=element_rect(colour = 'black', fill=NA, size = 3), legend.position = c(0.05,0.94))
dev.off()
## loadings plot -> we only display the names of those variables that are over ABS(#comp)>.15
loadings.plsda <- as.data.frame(cbind('label'=rownames(getLoadingMN(plsda.all.full)), getLoadingMN(plsda.all.full)))
loadings.plsda[,2:4] <- sapply(loadings.plsda[,2:4], function(x)as.numeric(as.character(x)))
tiff('fullmetabolomics_loadingsplsdaALL.tiff', res=600, height = 6500, width = 6500)
ggplot(loadings.plsda, aes(p1, p2)) + geom_point(aes(colour=ifelse(abs(loadings.plsda$p1)>.07|abs(loadings.plsda$p2)>.07, 'red', 'black'))) + 
  geom_hline(yintercept = 0) + geom_vline(xintercept = 0) +
  geom_text_repel(data=subset(loadings.plsda, abs(p1)>0.07 | abs(p2)>0.07), aes(label=label), cex=3, hjust=1.05, yjust=0) + 
  ylim(0.15, -0.15) + xlim(-0.15,0.15) + 
  theme(axis.line = element_line(size = 0.8, linetype = "solid"), axis.title = element_text(face = "bold"),
        axis.text = element_text(face = "bold"), panel.background = element_rect(fill = NA)) + scale_color_manual(values=c('black', 'red')) +
  guides(colour=F)
dev.off()
getSummaryDF(plsda.all.full)
## VIPs plot -> we generated two distinct plots, as there were 190 variables, plotting all of them in one and the
## 50 most contributing variables to allow and easier inspection.
VIPs <- as.data.frame(sort(getVipVn(plsda.all.full)))
names(VIPs) <- 'VIPscore'
VIPs$order <- factor(row.names(VIPs), levels=row.names(VIPs))

tiff('fullmetabolomics_VIPsALL.tiff', height = 90000, width = 4500, res=600)
ggplot(VIPs, aes(x=VIPscore, y=order)) + geom_point(size=3) + 
  theme(axis.line = element_line(size = 0.3, linetype = "solid"), axis.text = element_text(size = 11, face = "bold", colour = "black"), 
        panel.background = element_rect(fill = NA)) +labs(x = "VIP score", y = "Variable") +
  theme(axis.title = element_text(size = 14,  face = "bold"), axis.ticks.y = element_blank()) + geom_segment(aes(yend=order), xend=0) +labs(y = NULL) +
  scale_x_continuous(expand=c(0,0), limits = c(0,3))
dev.off()

VIPs.subset <- tail(VIPs, 50) ## plot only the 50 most contributing variables
tiff('fullmetabolomics_VIPs.tiff', height = 6000, width = 4500, res=600)
ggplot(VIPs.subset, aes(x=VIPscore, y=order)) + geom_point(size=3) + 
  theme(axis.line = element_line(size = 0.3, linetype = "solid"), axis.text = element_text(size = 11, face = "bold", colour = "black"), 
        panel.background = element_rect(fill = NA)) +labs(x = "VIP score", y = "Variable") +
  theme(axis.title = element_text(size = 14,  face = "bold"), axis.ticks.y = element_blank()) + geom_segment(aes(yend=order), xend=0) +labs(y = NULL) +
  scale_x_continuous(expand=c(0,0), limits = c(0,3))
dev.off()

## mixOmics analysis - multivariate rCCA correlation analysis between datasets -------------------------------------------

## Correlations between OTUs and full metabolomics dataset ---------------------------------------------------------------
otus <- read.csv2('level6.csv', header=T, sep='\t', row.names = 1)
otus[] <- sapply(otus, function(x)as.numeric(as.character(x)))
otus <- as.data.frame(t(otus))
otus <- otus[, colSums(otus!=0)>36]

# full.metabolomics <- readRDS('fullmetabolomics.rds')
# full.metabolomics <- readRDS('fullmetabolomics2.rds')
full.metabolomics <- readRDS('fullmetabolomics3.rds')

full.metabolomics <- full.metabolomics[!row.names(full.metabolomics)%in%'FIBR033C',]  ## we know FIBR033C did not work OK

otus <- otus[row.names(otus)%in%row.names(full.metabolomics),]
full.metabolomics <- full.metabolomics[row.names(full.metabolomics)%in%row.names(otus),]
row.names(full.metabolomics) == row.names(otus)
otus <- otus[match(row.names(full.metabolomics), row.names(otus)),]
row.names(full.metabolomics) == row.names(otus)

otus_deseq <- t(otus)
coldata <- data.frame('Status'=full.metabolomics$Status)
full.metabolomics$Status <- NULL
row.names(coldata) <- row.names(full.metabolomics)
ddstable <- DESeqDataSetFromMatrix(countData = otus_deseq, colData = coldata, design = ~Status)
gm_mean = function(x, na.rm=TRUE){  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
geoMeans = apply(counts(ddstable), 1, gm_mean)
statdds = estimateSizeFactors(ddstable, geoMeans = geoMeans)
statdds2 = DESeq(statdds, fitType="local")
otus_normal <- counts(statdds2, normalized = T)
otus_normal <- as.data.frame(t(otus_normal))

full.metabolomics[] <- sapply(full.metabolomics, NAfunc)
row.names(full.metabolomics) == row.names(otus_normal)
cor.dataset <- cbind(full.metabolomics, otus_normal)
cor.matrix <- cor(cor.dataset, method = 'spearman')
cor.result <- cor.mtest(cor.dataset)

tiff('spearman_cor_matrix1.tiff', height = 5000, width = 5000)
corrplot(cor.matrix, p.mat=cor.result$p, insig='blank', method = 'square', order='hclust') 
dev.off()

corrplot(cor.matrix, p.mat=cor.result$p, insig='blank', method = 'square', order='hclust') #,outline = T,
#         addgrid.col = 'darkgray', mar=c(4,0,4,0), cl.pos='b', tl.col='indianred4', tl.cex=.5, cl.cex=6)


tiff('spearman_cor_matrix_div.tiff', height = 5000, width = 7500)
corrplot(cor.matrix[1079:1153, 1:1078], p.mat=cor.result$p[1079:1153, 1:1078], insig='blank',  method = 'square', 
         outline = T, addgrid.col = 'darkgray', mar=c(4,0,4,0), cl.pos='b', tl.col='indianred4', tl.cex=.5, cl.cex=6) #order="hclust",
dev.off()


cor.heatmap <- cor.matrix
cor.heatmap <- cor.heatmap[, colnames(cor.heatmap)%in%colnames(otus_normal)]
cor.heatmap <- cor.heatmap[rownames(cor.heatmap)%in%colnames(full.metabolomics),]

 tiff('spearman_heatplot_corFIN.tiff', res=600, units='px', height = 15000, width = 7000)
#tiff('heatplot_cor2.tiff', res=600, units='px', height = 5000, width = 7000)
pdf('spearman_heatplot_corrFIN.pdf', width = 12, height = 20)
pdf('spearman_heatplot_corrFINLEGEND.pdf', width = 12, height = 10)
htplot <- heatplot(cor.heatmap, dend='both')
dev.off()

pdf('spearman_heatplot_corrFINfornames.pdf', width = 12, height = 150)
htplot <- heatplot(cor.heatmap, dend='both', cexRow=1)
dev.off()

mets <- c('X0.88_524.3714', 'X1.18_226.1062', 'X1.39_136.0490', 'X2.16_120.0665', 'X2.40_203.1509', 'X2.45_147.0769', 'X2.48_130.0508', 
          'X2.70_170.0931', 'X2.73_183.1109', 'X2.87_175.1196', 'X2.94_129.1029', 'X2.94_84.082', 'X3.00_155.0794')
bacteria <- c('f__Bifidobacteriaceae;g__Bifidobacterium', 'f__Paraprevotellaceae;g__Prevotella', 'f__Lachnospiraceae;g__Clostridium', 'f__Lachnospiraceae;g__Dorea',
              'f__Lachnospiraceae;g__Roseburia', 'f__Bacteroidaceae;g__Bacteroides', 'f__Rikenellaceae;g__Alistipes', 'f__Ruminococcaceae;g__Clostridium',
              'f__Ruminococcaceae;g__Papillibacter')
cor.heatmap2 <- cor.heatmap[rownames(cor.heatmap)%in%mets,]
cor.heatmap3 <- cor.heatmap[rownames(cor.heatmap)%in%mets, colnames(cor.heatmap)%in%bacteria]

tiff('spearman_heatplot_onlyIDmets.tiff', res=600, units='px', height=3000, width=6000)
pdf('spearman_heatplot_onlyIDmets.pdf', width = 15, height = 10)
heatplot(cor.heatmap2, dend='both', rowsep=c(7,9), colsep=c(16,28,47))
dev.off()

pdf('spearman_heatplot_MINI.pdf', width = 10, height = 10)
tiff('spearman_heatplot_MINI.tiff', res=600, units='px', height=3000, width=3000)
heatplot(cor.heatmap3, dend='both', rowsep=c(4, 7, 11), colsep=c(2,5))
dev.off()
