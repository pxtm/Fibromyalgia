summary(cytos.summ1)
pairs(cytos.summ2)
cytos.summ2 <- data.frame('Class'=cytos.summ1$Class, sapply(cytos.summ1[,2:53], NAfunc))
cytos.summ2$Class <- cytos.summ2$Class=='Fibromyalgia'
cor(cytos.summ2[,2:53])
model.full <- lm(Class~., data=cytos.summ2)
summary(model.full)
summary.aov(model.full)

step(model.full, direction = 'backward')
model.start <- lm(Class~IL.6, data=cytos.summ2)
step(model.start, direction='forward')
confint(model.full)

model.stats <- lm(Class~IL.6Ra+Axl+BST2+IL.1RI+IL.8, data=cytos.summ2)

par(mfrow=c(2,2))
plot(model.stats)
par(mfrow=c(1,1))

library(ROCR)


cytos.pca <- prcomp(cytos.summ2[,2:53], center=T, scale=T)
print(cytos.pca)
plot(cytos.pca, type='l')
summary(cytos.pca)
library(ggbiplot)
ggbiplot(cytos.pca, obs.scale = 1, var.scale = 1,
         groups = cytos.summ2$Class, ellipse = TRUE, circle = TRUE) +
  scale_color_discrete(name = '') +
  theme(legend.direction = 'horizontal', legend.position = 'top')
library(ropls)
cytos.pca2 <- opls(cytos.summ2[,2:53])
plot(cytos.pca2, typeVc='x-score', parAsColFcVn=cytos.summ2$Class)
cytos.pls <- opls(cytos.summ2[,2:53], cytos.summ2$Class, predI=3)
loadings.pls <- as.data.frame(cbind('label'=rownames(getLoadingMN(cytos.pls)), getLoadingMN(cytos.pls)))
loadings.pls[,2:4] <- sapply(loadings.pls[,2:4], function(x)as.numeric(as.character(x)))
ggplot(loadings.pls, aes(p1, p2)) + geom_point() + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) +
  geom_text_repel(aes(label=label), hjust=1.05, yjust=0) + ylim(0.3, -0.3) + xlim(-0.4,0.4) + theme(axis.line = element_line(size = 0.8, 
    linetype = "solid"), axis.title = element_text(face = "bold"), 
    axis.text = element_text(face = "bold"), 
    panel.background = element_rect(fill = NA))
  
getSummaryDF(cytos.pls)
VIPs <- as.data.frame(sort(getVipVn(cytos.pls)))
names(VIPs) <- 'VIPscore'
VIPs$order <- factor(row.names(VIPs), levels=row.names(VIPs))

ggplot(VIPs, aes(x=VIPscore, y=order)) + geom_point(size=3) + 
  theme(axis.line = element_line(size = 0.3, linetype = "solid"), axis.text = element_text(size = 11, face = "bold", colour = "black"), 
        panel.background = element_rect(fill = NA)) +labs(x = "VIP score", y = "Variable") +
  theme(axis.title = element_text(size = 14,  face = "bold"), axis.ticks.y = element_blank()) + geom_segment(aes(yend=order), xend=0) +labs(y = NULL)
                                                                             