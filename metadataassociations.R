## miRNAs 139-5p and 142-5p associations to any other variable ####
micros.fibro
pain
metadata.associations <- metadata[row.names(metadata)%in%row.names(micros.fibro),]

row.names(micros.fibro) == row.names(pain)
row.names(micros.fibro) == row.names(metadata.associations)

micros.associations <- micros.fibro[, c("hsa.mir.139.5p", "hsa.mir.142.5p")]
micros.associations <- data.frame(micros.associations, pain, metadata.associations)
micros.associations2 <- micros.associations

micros.associations$Age <- NULL
micros.associations$Status1 <- NULL
micros.associations$Status2 <- NULL
micros.associations$Subject <- NULL

pdf('associationplotstest2.pdf')
for (x in seq_along(micros.associations[,10:76])){
  if (length(levels(micros.associations[,x]))==2){
    #levs <- levels(micros.associations[,x])
    #comparisons = list(c(levs[1], levs[2]))
    print(ggplot(micros.associations, aes_string(x=names(micros.associations)[9+x], y='hsa.mir.139.5p')) + geom_boxplot() + 
           geom_beeswarm() + stat_compare_means())#(comparisons = comparisons))
  }
  else if ((length(levels(micros.associations[,x]))==3)){
    levs <- levels(micros.associations[,x])
    comparisons = list(c(levs[1], levs[2]),
                       c(levs[1], levs[3]),
                       c(levs[2], levs[3]))
    print(ggplot(micros.associations, aes_string(x=names(micros.associations)[9+x], y='hsa.mir.139.5p')) + geom_boxplot() + 
            geom_beeswarm() + stat_compare_means(comparisons = comparisons))
    
  }
  else {}
}
dev.off()


pdf('associationplotstest2.pdf')
for (x in seq_along(micros.associations[,11:76])){
    print(ggplot(micros.associations, aes_string(x=names(micros.associations)[10+x], y='hsa.mir.139.5p')) + geom_boxplot() + 
            geom_beeswarm() + stat_compare_means())#(comparisons = comparisons))
}
dev.off()

pdf('associationplotstest_1425p.pdf')
for (x in seq_along(micros.associations[,11:76])){
  print(ggplot(micros.associations, aes_string(x=names(micros.associations)[10+x], y='hsa.mir.142.5p')) + geom_boxplot() + 
          geom_beeswarm() + stat_compare_means())#(comparisons = comparisons))
}
dev.off()



for (i in seq_along(micros.associations)){
  levels(micros.associations[,i])[levels(micros.associations[,i])==''] <- NA
}


sink('miR-139-5p_associations.txt')
for (x in seq_along(micros.associations[,10:76])){
    if (length(levels(micros.associations[,x]))==2){
      tryCatch({print(paste0(names(micros.associations[x]), ' t-test: ', 
                             t.test(micros.associations$hsa.mir.139.5p~micros.associations[,x], na.action = na.omit)$p.val))
      }, error=function(a){cat(names(micros.associations[x]), ' ERROR: ', conditionMessage(a), '\n')})
  }
  else if (length(levels(micros.associations[,x]))==3){
    tryCatch({
      aovmodel <- aov(hsa.mir.139.5p~micros.associations[,x], data=micros.associations);
      print(paste0(names(micros.associations)[x], ' ANOVA: ', summary(aovmodel)[[1]][['Pr(>F)']][1]));
      tryCatch({
        print(TukeyHSD(aovmodel)[[1]][,4])
      }, error=function(b){cat('ERROR: ', conditionMessage(b), '\n')})
    }, error=function(a){cat('ERROR: ', conditionMessage(a), '\n')})
  }
  else {}
}
sink()

sink('miR-142-5p_associations.txt')
for (x in seq_along(micros.associations[,10:76])){
  if (length(levels(micros.associations[,x]))==2){
    tryCatch({print(paste0(names(micros.associations[x]), ' t-test: ', 
                           t.test(micros.associations$hsa.mir.142.5p~micros.associations[,x], na.action = na.omit)$p.val))
    }, error=function(a){cat(names(micros.associations[x]), ' ERROR: ', conditionMessage(a), '\n')})
  }
  else if (length(levels(micros.associations[,x]))==3){
    tryCatch({
      aovmodel <- aov(hsa.mir.142.5p~micros.associations[,x], data=micros.associations);
      print(paste0(names(micros.associations)[x], ' ANOVA: ', summary(aovmodel)[[1]][['Pr(>F)']][1]));
      tryCatch({
        print(TukeyHSD(aovmodel)[[1]][,4])
      }, error=function(b){cat('ERROR: ', conditionMessage(b), '\n')})
    }, error=function(a){cat('ERROR: ', conditionMessage(a), '\n')})
  }
  else {}
}
sink()

## boxplots for hsa-miR-139-5p: alcohol, CHCM, Glomerular Filtrate, Vegetables, Benzodiazepinas ####
tiff('mir139associations.tiff', height = 5500, width = 3000, res=600, units = 'px')
multiplot(
  ggplot(data=subset(micros.associations, !is.na(Alcohol)), aes(x=Alcohol, y=hsa.mir.139.5p)) + geom_boxplot() + geom_beeswarm() + 
    stat_compare_means(method='t.test', label.x = 1, label.y = 0.5, aes(label=paste0('p= ', round(..p.adj.., 3))),fontface='bold.italic' ) + 
    theme(axis.line = element_line(size = 0.6, linetype = "solid"), axis.ticks = element_line(colour = "gray16"),
          axis.title = element_text(size = 13, face = "bold"), axis.text = element_text(size = 11, face = "bold", colour = "black"),
          panel.background = element_rect(fill = NA)) +labs(y = "hsa-mir-139-5p"),
  ggplot(data=subset(micros.associations, !is.na(CHCM)), aes(x=CHCM, y=hsa.mir.139.5p)) + geom_boxplot() + geom_beeswarm() + 
    stat_compare_means(method='t.test', label.x = 1, label.y = 0.5, aes(label=paste0('p= ', round(..p.adj.., 3))),fontface='bold.italic') + 
    theme(axis.line = element_line(size = 0.6, linetype = "solid"), axis.ticks = element_line(colour = "gray16"), 
          axis.title = element_text(size = 13, face = "bold"), axis.text = element_text(size = 11, face = "bold", colour = "black"), 
          panel.background = element_rect(fill = NA)) +labs(y = "hsa-mir-139-5p"),
  ggplot(data=subset(micros.associations, !is.na(Vegetables)), aes(x=Vegetables, y=hsa.mir.139.5p)) + geom_boxplot() + geom_beeswarm() + 
    stat_compare_means(method='t.test', label.x = 1, label.y = 0.5, aes(label=paste0('p= ', round(..p.adj.., 3))),fontface='bold.italic') + 
    theme(axis.line = element_line(size = 0.6, linetype = "solid"), axis.ticks = element_line(colour = "gray16"), 
          axis.title = element_text(size = 13, face = "bold"), axis.text = element_text(size = 11, face = "bold", colour = "black"), 
          panel.background = element_rect(fill = NA)) +labs(y = "hsa-mir-139-5p"),
  ggplot(data=subset(micros.associations, !is.na(Benzodiazepina)), aes(x=Benzodiazepina, y=hsa.mir.139.5p)) + geom_boxplot() + geom_beeswarm() + 
    stat_compare_means(method='t.test', label.x = 1, label.y = 0.5, aes(label=paste0('p= ', round(..p.adj.., 3))),fontface='bold.italic') + 
    theme(axis.line = element_line(size = 0.6, linetype = "solid"), axis.ticks = element_line(colour = "gray16"), 
          axis.title = element_text(size = 13, face = "bold"), axis.text = element_text(size = 11, face = "bold", colour = "black"), 
          panel.background = element_rect(fill = NA)) +labs(y = "hsa-mir-139-5p"),
  ggplot(data=subset(micros.associations, !is.na(Glomerular.filtrate)), aes(x=Glomerular.filtrate, y=hsa.mir.139.5p)) + geom_boxplot() + geom_beeswarm() + 
    stat_compare_means(method='anova', label.y = 4, aes(label=paste0('p= ', round(..p.adj.., 3))),fontface='bold.italic') + 
    stat_compare_means(comparisons = list(c('<80', '>120'), c('<80', '>80'), c('>120', '>80')), aes(label=paste0('p= ', round(..p.adj.., 3))),fontface='bold.italic') + 
    theme(axis.line = element_line(size = 0.6, linetype = "solid"), axis.ticks = element_line(colour = "gray16"), 
          axis.title = element_text(size = 13, face = "bold"), axis.text = element_text(size = 11, face = "bold", colour = "black"), 
          panel.background = element_rect(fill = NA)) +labs(y = "hsa-mir-139-5p"), cols=2)
dev.off()

## boxplots for hsa-miR-142-5p: Excercise, Fiber, Hospital, Antiepileptico_GABA ####
tiff('mir142associations.tiff', height = 3500, width = 2500, res=600, units = 'px')
multiplot(
  ggplot(data=subset(micros.associations, !is.na(Excercise)), aes(x=Excercise, y=hsa.mir.142.5p)) + geom_boxplot() + geom_beeswarm() + 
            stat_compare_means(method='t.test', label.x = 2, label.y = 2.75, aes(label=paste0('p= ', round(..p.adj.., 3))),fontface='bold.italic' ) + 
            theme(axis.line = element_line(size = 0.6, linetype = "solid"), axis.ticks = element_line(colour = "gray16"),
                  axis.title = element_text(size = 13, face = "bold"), axis.text = element_text(size = 11, face = "bold", colour = "black"),
                  panel.background = element_rect(fill = NA)) +labs(y = "hsa-mir-142-5p"),
  ggplot(data=subset(micros.associations, !is.na(Fiber)), aes(x=Fiber, y=hsa.mir.142.5p)) + geom_boxplot() + geom_beeswarm() + 
    stat_compare_means(method='t.test', label.x = 2, label.y = 2.75, aes(label=paste0('p= ', round(..p.adj.., 3))),fontface='bold.italic' ) + 
    theme(axis.line = element_line(size = 0.6, linetype = "solid"), axis.ticks = element_line(colour = "gray16"),
          axis.title = element_text(size = 13, face = "bold"), axis.text = element_text(size = 11, face = "bold", colour = "black"),
          panel.background = element_rect(fill = NA)) +labs(y = "hsa-mir-142-5p"),
  ggplot(data=subset(micros.associations, !is.na(Hospital)), aes(x=Hospital, y=hsa.mir.142.5p)) + geom_boxplot() + geom_beeswarm() + 
    stat_compare_means(method='t.test', label.x = 2, label.y = 2.75, aes(label=paste0('p= ', round(..p.adj.., 3))),fontface='bold.italic' ) + 
    theme(axis.line = element_line(size = 0.6, linetype = "solid"), axis.ticks = element_line(colour = "gray16"),
          axis.title = element_text(size = 13, face = "bold"), axis.text = element_text(size = 11, face = "bold", colour = "black"),
          panel.background = element_rect(fill = NA)) +labs(y = "hsa-mir-142-5p"),
  ggplot(data=subset(micros.associations, !is.na(Antiepileptico_GABA)), aes(x=Antiepileptico_GABA, y=hsa.mir.142.5p)) + geom_boxplot() + geom_beeswarm() + 
    stat_compare_means(method='t.test', label.x = 2, label.y = 2.75, aes(label=paste0('p= ', round(..p.adj.., 3))),fontface='bold.italic' ) + 
    theme(axis.line = element_line(size = 0.6, linetype = "solid"), axis.ticks = element_line(colour = "gray16"),
          axis.title = element_text(size = 13, face = "bold"), axis.text = element_text(size = 11, face = "bold", colour = "black"),
          panel.background = element_rect(fill = NA)) +labs(y = "hsa-mir-142-5p"),
cols=2)
dev.off()
