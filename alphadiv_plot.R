library(xlsx); library(ggplot2)
## alpha diversity
alpha <- read.xlsx2('alpha_div.xlsx', header=T, sheetIndex = 1)
alpha$Subject <- NULL
alpha$Status <- NULL
alpha$Sex.1 <- NULL
row.names(alpha) <- alpha$sample
alpha$sample <- NULL
alpha[,1:6] <- sapply(alpha[,1:6], function(t)as.numeric(as.character(t)))
alpha1=alpha
alpha1[,7:22] <- sapply(alpha1[,7:22], function(f){is.na(f)<-which(f == '');f}) 


alpha_status <- alpha[, c( 'ace', 'chao1', 'observed_otus', 'shannon', 'simpson', 'faith_pd', 'Status2')]
alpha_status <- gather(alpha_status, key='tests', value = 'value', 1:6) 

alpha_status$value <- as.numeric(alpha_status$value)
ggplot(alpha_status, aes(Status2, value))+geom_boxplot()+geom_beeswarm()+theme_bw()+
  theme(axis.title.x = element_blank()) + facet_wrap(~tests, scales='free_y') +
  stat_compare_means(method='t.test', aes(label=paste0('p = ', round(..p.., digits = 2))), fontface=4, label.x = 2.25) + 
  theme(axis.text.y = element_text(colour = "black")) +labs(title = "Status")
                     
library(tidyverse)                     
                     
alpha1 <- gather(alpha, key='tests', value='value', 1:6)                     
plots <- list()
for (i in seq(1:16)){
  alp_plot <- data.frame('var'=alpha1[,i], 'tests'=alpha1$tests, 'value'=alpha1$value)
  alp_plot <- subset(alp_plot, alp_plot$var!="")
  plots[[i]] <- ggplot(alp_plot, aes(var, value))+geom_boxplot(na.rm = T)+geom_beeswarm()+theme_bw()+
    theme(axis.title.x = element_blank()) + facet_wrap(~tests, scales='free_y') +
   stat_compare_means(method='t.test', aes(label=paste0('p = ', round(..p.., digits = 2))), fontface=4, label.x = 2.25) + 
   theme(axis.text.y = element_text(colour = "black")) +labs(title = names(alpha1[i]))
  } 

pdf('diversity_plots.pdf', width = 12)
for (i in length(plots)){
  print(plots)
}
dev.off()  


#pdf('status_diversity.pdf', width = 12, height = 12)
tiff('status_diversity.tiff', res=600, units='px', width = 3750, height = 6000)
ggplot(alpha1, aes(Status2, value))+geom_boxplot(na.rm = T)+geom_beeswarm()+theme_bw()+
  theme(axis.title.x = element_blank()) + facet_wrap(~tests, scales='free_y', ncol = 2) +
  stat_compare_means(method='t.test', aes(label=paste0('p = ', round(..p.., digits = 2))), fontface=4, label.x = 2.15) + 
  theme(axis.text.y = element_text(colour = "black", size=13, face='bold'))  + theme(axis.title = element_text(size = 13, 
    face = "bold"), plot.title = element_text(face = "bold"), 
    strip.background = element_rect(fill = "gray98", 
        colour = "gray0"), strip.text = element_text(face = "bold", 
        colour = "gray0"),
    axis.text.x=element_text(colour = 'black', size=13, face='bold')) +labs(y = NULL)
dev.off()
