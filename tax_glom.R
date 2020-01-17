library(phyloseq); library(microbiome); library(VennDiagram); library(ggplot2); library(RColorBrewer); library(DESeq2)
setwd('closedref') # to use the CLOSED REFERENCE CLUSTERED OTUs

## upload data --------------------------------------------------------------------------------------------------------
otus <- read.csv2('otus_tax.csv', header = T, row.names = 1, dec='.', stringsAsFactors = F)
otus$Confidence <- NULL
taxonomy <- as.character(otus$Taxon)
names(taxonomy) <- row.names(otus)
otus$Taxon <- NULL

metadata <- read.csv2('metadata.csv', header = T, row.names=1)
metadata <- subset(metadata, metadata$Status2!='Related')
# metadata <- subset(metadata, metadata$Status!='Related')
metadata$Status2 <- droplevels(metadata$Status2)
#metadata$Status <- droplevels(metadata$Status)
remove <- c('FIBR025C', 'FIBR025F') #, 'FIBR013C', 'FIBR014C', 'FIBR021C', 'FIBR028C', 'FIBR030C', 'FIBR033C')
metadata <- metadata[!rownames(metadata) %in% remove, ]
otus <- as.data.frame(t(otus), stringsAsFactors = F)
otus <- otus[row.names(otus)%in%row.names(metadata),]
rm(remove)


## mount phyloseq object ----------------------------------------------------------------------------------------------
## otu_table: numeric matrix, specify species in row or column
rnames_otutb <- row.names(otus)
cnames_otutb <- names(otus)
otus_m <- t(as.matrix(otus))
rownames(otus_m) <- cnames_otutb
colnames(otus_m) <- rnames_otutb
#rm(otus, rnames_otutb, cnames_otutb)
# sample_data: data.frame, rownames must match sample names in otu_Table to combine with phyloseq
# tax_table: character matrix. Rownames must match OTU names (otu_table) of the otu_table.
# write.table(taxonomy, 'tax_to_physeq.csv') #modify table in order to fit the required format
taxmat <- read.csv2('tax_to_physeq.csv', header=T, row.names=1)
taxmat <- as.matrix(taxmat)
rownames(taxmat) <- rownames(otus_m)
colnames(taxmat) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")

# convert to phyloseq format
OTU<-otu_table(otus_m, taxa_are_rows = T)
TAX<-tax_table(taxmat)
#rm(otus_m, taxmat)

# add metadata
#metadata2<-metadata
#metadata2<-subset(metadata2, metadata2$Status!='Unsure')
physeq<-phyloseq(OTU, TAX, sample_data(metadata))


## summarize taxa to specific level -----------------------------------------------------------------------------------
## phylum level analysis ==============================================================================================
physeq.phylum <- tax_glom(physeq, taxrank = rank_names(physeq)[2])
phylum.ctl <- subset_samples(physeq.phylum, Status2%in%c('Control'))
phylum.fib <- subset_samples(physeq.phylum, Status2%in%c('Fibromyalgia'))
 ## relative frequency
phylum.rel <- microbiome::transform(physeq.phylum, 'compositional')
head(prevalence(phylum.rel, detection = 1/100, sort = TRUE))
head(prevalence(phylum.rel, detection = 1/100, sort = TRUE, count = TRUE))

 ## control phylum core
p.con <- prevalence(phylum.ctl, detection = 0, sort  = T) #population frequencies
p.con.compo<-microbiome::transform(phylum.ctl, 'compositional') #transform to relative abundances
ctl.core<-core(p.con.compo, detection=.2/100, prevalence = 50/100) #taxa with over 50%  prevalence at .2% relative abundance
ctl.core.otus<-taxa(ctl.core)
 ## fibro core
p.fib <- prevalence(phylum.fib, detection = 0, sort  = T) #population frequencies
p.fib.compo<-microbiome::transform(phylum.fib, 'compositional') #transform to relative abundances
fib.core<-core(p.fib.compo, detection=.2/100, prevalence = 50/100) #taxa with over 50%  prevalence at .2% relative abundance
fib.core.otus<-taxa(fib.core)

phylum.core <- merge_phyloseq(ctl.core, fib.core)
phylum.core

phylum.names <- data.frame(tax_table(phylum.core)[,2])
colorCount = 4

 ## composition plot
p <- plot_composition(microbiome::transform(phylum.core, "compositional"), plot.type = "barplot", sample.sort = "Status2") + 
  scale_fill_manual(values=getPalette(colorCount), labels=c('Proteobacteria', 'Firmicutes', 'Bacteroidetes', 'Actinobacteria')) +
  scale_y_continuous(expand = c(0,0), limits=c(0,1.01)) + coord_flip()
p + theme(axis.text.x = element_text(size=12, angle = 0,  face='bold'), axis.title.y = element_blank(), 
          axis.title.x =  element_text(face='bold', size=12), axis.text.y = element_text(size=12, face='bold'),
          panel.grid.major = element_blank(), legend.text = element_text(size=10, face='bold'), legend.title = element_blank()) + 
  ylab('Relative abundance')  

 ## DESeq2 analysis
phylum.dds <- phyloseq_to_deseq2(physeq.phylum, ~Status2)
gm_mean <- function(x, na.rm=TRUE){  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
geoMeans <-  apply(counts(phylum.dds), 1, gm_mean)
phylum.dds <-  estimateSizeFactors(phylum.dds, geoMeans = geoMeans)
phylum.dds = DESeq(phylum.dds, fitType="local")

res.phylum <- results(phylum.dds, cooksCutoff = F)
phy.sigtab <- res.phylum[which(res.phylum$padj<.05),]
## family level analysis ==============================================================================================
physeq.family <- tax_glom(physeq, taxrank = rank_names(physeq)[5])
family.ctl <- subset_samples(physeq.family, Status2%in%c('Control'))
family.fib <- subset_samples(physeq.family, Status2%in%c('Fibromyalgia'))
## relative frequency
family.rel <- microbiome::transform(physeq.family, 'compositional')
head(prevalence(family.rel, detection = 1/100, sort = TRUE))
head(prevalence(family.rel, detection = 1/100, sort = TRUE, count = TRUE))

## control family core
p.con <- prevalence(family.ctl, detection = 0, sort  = T) #population frequencies
p.con.compo<-microbiome::transform(family.ctl, 'compositional') #transform to relative abundances
ctl.core<-core(p.con.compo, detection=.2/100, prevalence = 50/100) #taxa with over 50%  prevalence at .2% relative abundance
ctl.core.otus<-taxa(ctl.core)
## fibro family core
p.fib <- prevalence(family.fib, detection = 0, sort  = T) #population frequencies
p.fib.compo<-microbiome::transform(family.fib, 'compositional') #transform to relative abundances
fib.core<-core(p.fib.compo, detection=.2/100, prevalence = 50/100) #taxa with over 50%  prevalence at .2% relative abundance
fib.core.otus<-taxa(fib.core)

family.core <- merge_phyloseq(ctl.core, fib.core)
family.core

family.names <- data.frame(tax_table(family.core)[,5])
colorCount = 12

## composition plot
tiff('family_core_summ.tiff', width = 12000, height =  9500, res = 600, units = 'px', compression = 'lzw')
p <- plot_composition(microbiome::transform(family.core, "compositional"), plot.type = "barplot", sample.sort = "Status2") + 
  scale_fill_manual(values = c('#A6CEE3', '#1F78B4', '#FB9A99', '#E31A1C', '#B2DF6A', '#33A02C', '#CAB2D6', '#6A3D9A','#FDBF6F',
                               '#FF7F00','#B3B6B7','#484949'), labels=tax_table(family.core)[,5]) +
  scale_y_continuous(expand = c(0,0), limits=c(0,1.01)) + coord_flip()
p + theme(axis.text.x = element_text(size=12, angle = 0,  face='bold'), axis.title.y = element_blank(), 
          axis.title.x =  element_text(face='bold', size=12), axis.text.y = element_text(size=12, face='bold'),
          panel.grid.major = element_blank(), legend.text = element_text(size=10, face='bold'), legend.title = element_blank()) + 
  ylab('Relative abundance') + ggtitle(label = 'Family level')  
dev.off()

## DESeq2 analysis
family.dds <- phyloseq_to_deseq2(physeq.family, ~Status2)
gm_mean <- function(x, na.rm=TRUE){  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
geoMeans <-  apply(counts(family.dds), 1, gm_mean)
family.dds <-  estimateSizeFactors(family.dds, geoMeans = geoMeans)
family.dds = DESeq(family.dds, fitType="local")

res.family <- results(family.dds, cooksCutoff = F)
fam.sigtab <- res.family[which(res.family$padj<.05),]


## genus level analysis ===============================================================================================
physeq.genus <- tax_glom(physeq, taxrank = rank_names(physeq)[6])
genus.ctl <- subset_samples(physeq.genus, Status2%in%c('Control'))
genus.fib <- subset_samples(physeq.genus, Status2%in%c('Fibromyalgia'))
## relative frequency
genus.rel <- microbiome::transform(physeq.genus, 'compositional')
head(prevalence(genus.rel, detection = 1/100, sort = TRUE))
head(prevalence(genus.rel, detection = 1/100, sort = TRUE, count = TRUE))

## control genus core
p.con <- prevalence(genus.ctl, detection = 0, sort  = T) #population frequencies
p.con.compo<-microbiome::transform(genus.ctl, 'compositional') #transform to relative abundances
ctl.core<-core(p.con.compo, detection=.2/100, prevalence = 50/100) #taxa with over 50%  prevalence at .2% relative abundance
ctl.core.otus<-taxa(ctl.core)
## fibro genus core
p.fib <- prevalence(genus.fib, detection = 0, sort  = T) #population frequencies
p.fib.compo<-microbiome::transform(genus.fib, 'compositional') #transform to relative abundances
fib.core<-core(p.fib.compo, detection=.2/100, prevalence = 50/100) #taxa with over 50%  prevalence at .2% relative abundance
fib.core.otus<-taxa(fib.core)

genus.core <- merge_phyloseq(ctl.core, fib.core)
genus.core

genus.names <- data.frame(tax_table(genus.core)[,6])
colorCount = 15

## composition plot
# p <- plot_composition(microbiome::transform(genus.core, "compositional"), plot.type = "barplot", sample.sort = "Status2") + 
#   scale_fill_manual(values=getPalette(colorCount), labels=genus.names$Genus) +
#   scale_y_continuous(expand = c(0,0), limits=c(0,1.01)) + coord_flip()
# p + theme(axis.text.x = element_text(size=12, angle = 0,  face='bold'), axis.title.y = element_blank(), 
#           axis.title.x =  element_text(face='bold', size=12), axis.text.y = element_text(size=12, face='bold'),
#           panel.grid.major = element_blank(), legend.text = element_text(size=10, face='bold'), legend.title = element_blank()) + 
#   ylab('Relative abundance') + ggtitle(label = 'genus level')  
tiff('GENUS_core_summ.tiff', width = 12000, height =  9500, res = 600, units = 'px', compression = 'lzw')
p <- plot_composition(microbiome::transform(genus.core, "compositional"), plot.type = "barplot", sample.sort = "Status2") + 
  scale_fill_manual(values = c('#A6CEE3', '#1F78B4', '#FB9A99', '#EFB26C', '#B2DF6A', '#6CEFBA', '#8612B7', '#CAB2D6', 
                               '#F7C27C', '#EFA3EA', '#B3B6B7', '#A3E3FA', '#484949', '#B04859', 'bisque', '#FCF891'), labels=genus.names$Genus) +
  scale_y_continuous(expand = c(0,0), limits=c(0,1.01)) + coord_flip()
p + theme(axis.text.x = element_text(size=12, angle = 0,  face='bold'), axis.title.y = element_blank(), 
          axis.title.x =  element_text(face='bold', size=12), axis.text.y = element_text(size=12, face='bold'),
          panel.grid.major = element_blank(), legend.text = element_text(size=10, face='bold'), legend.title = element_blank()) + 
  ylab('Relative abundance') + ggtitle(label = 'Genus level')  
dev.off()

## DESeq2 analysis 
otus.gen <- otu_table(physeq.genus)
otus.gen <- as.data.frame(t(otus.gen))
otus.gen.fib <- otus.gen[row.names(otus.gen)%in%row.names(subset(metadata, metadata$Status2=='Fibromyalgia')),]
otus.gen.ctl <- otus.gen[row.names(otus.gen)%in%row.names(subset(metadata, metadata$Status2=='Control')),]

otus.gen.ctl <- as.data.frame(t(otus.gen.ctl), stringsasFactors=F)
otus.gen.ctl[] <- sapply(otus.gen.ctl, function(x){as.numeric(as.character(x))})
otus.gen.ctl <- otus.gen.ctl[rowSums(otus.gen.ctl>0)>=dim(otus.gen.ctl)[2]/3,]
genus.ctl2 <- prune_taxa(row.names(otus.gen.ctl), genus.ctl)

otus.gen.fib <- as.data.frame(t(otus.gen.fib), stringsasFactors=F)
otus.gen.fib[] <- sapply(otus.gen.fib, function(x){as.numeric(as.character(x))})
otus.gen.fib <- otus.gen.fib[rowSums(otus.gen.fib>0)>=dim(otus.gen.fib)[2]/3,]
genus.fib2 <- prune_taxa(row.names(otus.gen.fib), genus.fib)

physeq.genus.cl <- merge_phyloseq(genus.ctl2, genus.fib2)




genus.dds <- phyloseq_to_deseq2(physeq.genus.cl, ~Status2)
gm_mean <- function(x, na.rm=TRUE){  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
geoMeans <-  apply(counts(genus.dds), 1, gm_mean)
genus.dds <-  estimateSizeFactors(genus.dds, geoMeans = geoMeans)
genus.dds = DESeq(genus.dds, fitType="local")

res.genus <- results(genus.dds, cooksCutoff = F)
gen.sigtab <- res.genus[which(res.genus$padj<.05),]

physeq.species <- tax_glom(physeq, taxrank = rank_names(physeq)[7])

## genus OTUs - pain correlations
genus.corr <- as.data.frame(otu_table(physeq.genus.cl))
genus.annot <- as.data.frame(tax_table(physeq.genus.cl))

pain.indexs <- read.csv2('pain_indices.csv', header=T, row.names=1)
genus.corr <- genus.corr[,names(genus.corr)%in%row.names(pain.indexs)]
pain.indexs <- pain.indexs[row.names(pain.indexs)%in%names(genus.corr), ]
genus.corr <- as.data.frame(t(genus.corr))
row.names(genus.corr) == row.names(pain.indexs)
genus.corr <- genus.corr[match(row.names(pain.indexs), row.names(genus.corr)),]
genus.annot[,5:6]
rnames.cor <- c()
for (i in 1:nrow(genus.annot)){
  rnames.cor[i] <- paste(genus.annot[i,5], genus.annot[i, 6], sep=' ')
}
names(genus.corr) <- rnames.cor
genus.corr <- data.frame(genus.corr, pain.indexs)
genus.corr$Erysipelotrichaceae.Eubacterium <- NULL
genus.corr[] <- sapply(genus.corr, NAfunc)


library(corrplot)
genus.pain.cor <- cor(genus.corr)
genus.pain.res <- cor.mtest(genus.corr)
tiff('genusvspain.tiff', width = 3500, height = 5500, res=600, units = 'px')
corrplot(genus.pain.cor[1:25,26:33], p.mat = genus.pain.res$p[1:25,26:33],insig='blank', 
         tl.co='black')
dev.off()
