library(Seurat)
library(ggplot2)

stmerged <- readRDS('gloms_manuscript_020824.RDS')

coltable <- read.csv('../identify_glom_cell_types/colors_new_types.csv')
coltable <- rbind(coltable,
                  list(New.types = 'Other',Colors = '#DDDDDD'))
coltable[3,1] <- 'VSM/P'
coltable[6,1] <- 'EC'
coltable <- coltable[c(1,3,6,9),]
rownames(coltable) <- coltable[,1]


pred <- as.data.frame(reshape2::melt(stmerged@assays$predsubclassl1@data))
pred <- pred[pred$Var1 != 'max',]
pred$celltype <- pred$Var1
pred$celltype <- as.character(pred$celltype)
pred[!pred$celltype %in% coltable$New.types,"celltype"] <- 'Other'
pred$condition <- stmerged@meta.data[pred$Var2,'condition']

fisher <- as.data.frame(matrix(0,ncol=2,nrow = 3,
                               dimnames = list(coltable$New.types[c(1,3,2)],c('p.value','odds.raio'))))
forest <- as.data.frame(matrix(0,ncol=3,nrow = 3,
                               dimnames = list(coltable$New.types[c(1,3,2)],c('low','est','high'))))
for (ct in coltable$New.types[c(1,3,2)]){
  fishmat <- matrix(c(sum(pred[pred$celltype==ct & pred$condition == 'DKD',"value"]),
                      sum(pred[pred$celltype==ct & pred$condition == 'Ref',"value"]),
                      sum(pred[pred$celltype!=ct & pred$condition == 'DKD',"value"]),
                      sum(pred[pred$celltype!=ct & pred$condition == 'Ref',"value"])),ncol=2)
  test <- fisher.test(fishmat)
  fisher[ct,"p.value"] <- test$p.value
  fisher[ct,"odds.raio"] <- test$estimate
  forest[ct,'est'] <- log(test$estimate)
  forest[ct,'low'] <- log(test$estimate) - 1.96*sqrt(1/fishmat[1,1] + 1/fishmat[2,1] + 
                                                      1/fishmat[1,2] + 1/fishmat[2,2])
  forest[ct,'high'] <- log(test$estimate) + 1.96*sqrt(1/fishmat[1,1] + 1/fishmat[2,1] + 
                                                      1/fishmat[1,2] + 1/fishmat[2,2])
}

forest$cell <- factor(rownames(forest),levels=c('EC','VSM/P','POD'))
pdf('forest_l1_proportions.pdf',width = 4,height = 3)
ggplot(data=forest, aes(x=cell, y=est, ymin=low, ymax=high,color=cell)) +
  geom_pointrange(size = .7,linewidth = 1,shape = 2) + 
  geom_hline(yintercept=0, lty=2) +  # add a dotted line at x=1 after flip
  coord_flip() +  # flip coordinates (puts labels on y axis)
  xlab("Cell type") + ylab("Ln(OR) (95% CI)")+
  scale_color_manual(values = c('#ff8c00','#483d8b','#db7295'),
                     name = element_blank(),
                     guide = guide_legend(reverse = TRUE))+
  theme_classic()+
  theme(panel.border = element_rect(color = "black", fill=NA, linewidth = 1.5),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10))
dev.off()

pred <- aggregate(pred$value,by=list(spot=pred$Var2,celltype=pred$celltype),FUN=sum)
pred$condition <- stmerged@meta.data[pred$spot,'condition']
pred <- aggregate(pred$x,by=list(celltype=pred$celltype,Condition=pred$condition),FUN=mean)
pred <- pred[pred$celltype %in% coltable$New.types[1:3],]
pred$celltype <- factor(pred$celltype,level=c(coltable$New.types[c(1,3,2)]))
pred$Condition <- factor(pred$Condition,level=c('Ref','DKD'))

pdf('barplot_l1_proportions.pdf',width = 6,height = 3)
ggplot(pred,aes(x=celltype,y=x,fill=Condition))+
  geom_bar(stat = 'identity',position = 'dodge')+
  scale_fill_manual(values = c('#8DCFFF','#FFA829'))+
  theme_classic()+
  theme(panel.border = element_rect(color = "black", fill=NA, linewidth = 1.5),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10))+
  ylab('Proportion')+
  xlab('Condition')
dev.off()

write.csv(fisher,'fishers_l1_proportion.csv',quote = F,row.names = T)


coltable <- read.csv('../identify_glom_cell_types/colors_new_types.csv')
coltable <- rbind(coltable,
                  list(New.types = 'Other',Colors = '#DDDDDD'))
rownames(coltable) <- coltable[,1]

pred <- as.data.frame(reshape2::melt(stmerged@assays$pred.new.types@data))
pred <- pred[pred$Var1 != 'max',]
pred$celltype <- pred$Var1
pred$celltype <- as.character(pred$celltype)
pred[!pred$celltype %in% coltable$New.types,"celltype"] <- 'Other'
pred$condition <- stmerged@meta.data[pred$Var2,'condition']

fisher <- as.data.frame(matrix(0,ncol=2,nrow = 9,
                               dimnames = list(coltable$New.types,c('p.value','odds.raio'))))
forest <- as.data.frame(matrix(0,ncol=3,nrow = 9,
                               dimnames = list(coltable$New.types,c('low','est','high'))))
for (ct in coltable$New.types){
  fishmat <- matrix(c(sum(pred[pred$celltype==ct & pred$condition == 'DKD',"value"]),
                      sum(pred[pred$celltype==ct & pred$condition == 'Ref',"value"]),
                      sum(pred[pred$celltype!=ct & pred$condition == 'DKD',"value"]),
                      sum(pred[pred$celltype!=ct & pred$condition == 'Ref',"value"])),ncol=2)
  test <- fisher.test(fishmat)
  fisher[ct,"p.value"] <- test$p.value
  fisher[ct,"odds.raio"] <- test$estimate
  forest[ct,'est'] <- log(test$estimate)
  forest[ct,'low'] <- log(test$estimate) - 1.96*sqrt(1/fishmat[1,1] + 1/fishmat[2,1] + 
                                                       1/fishmat[1,2] + 1/fishmat[2,2])
  forest[ct,'high'] <- log(test$estimate) + 1.96*sqrt(1/fishmat[1,1] + 1/fishmat[2,1] + 
                                                        1/fishmat[1,2] + 1/fishmat[2,2])
}


forest$cell <- factor(rownames(forest),levels=rev(coltable$New.types))
pdf('forest_newtypes_proportions.pdf',width = 4,height = 3)
ggplot(data=forest, aes(x=cell, y=est, ymin=low, ymax=high,color=cell)) +
  geom_pointrange(size = .7,linewidth = 1,shape = 2) + 
  geom_hline(yintercept=0, lty=2) +  # add a dotted line at x=1 after flip
  coord_flip() +  # flip coordinates (puts labels on y axis)
  xlab("Cell type") + ylab("Ln(OR) (95% CI)")+
  scale_x_discrete(limits=rev(coltable$New.types[c(1,3:8)]))+
  scale_color_manual(values = rev(coltable$Colors[c(1,3:8)]),
                     name = element_blank(),
                     guide = guide_legend(reverse = TRUE),
                     limits=rev(coltable$New.types[c(1,3:8)]))+
  theme_classic()+
  theme(panel.border = element_rect(color = "black", fill=NA, linewidth = 1.5),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10))
dev.off()


pred <- aggregate(pred$value,by=list(spot=pred$Var2,celltype=pred$celltype),FUN=sum)
pred$condition <- stmerged@meta.data[pred$spot,'condition']
pred <- aggregate(pred$x,by=list(celltype=pred$celltype,Condition=pred$condition),FUN=mean)
pred <- pred[pred$celltype %in% coltable$New.types[1:8],]
pred$celltype <- factor(pred$celltype,level=c(coltable$New.types))
pred$Condition <- factor(pred$Condition,level=c('Ref','DKD'))

pdf('barplot_newtypes_proportions.pdf',width = 6,height = 3)
ggplot(pred,aes(x=celltype,y=x,fill=Condition))+
  geom_bar(stat = 'identity',position = 'dodge')+
  scale_fill_manual(values = c('#8DCFFF','#FFA829'))+
  theme_classic()+
  theme(panel.border = element_rect(color = "black", fill=NA, linewidth = 1.5),
        axis.text.x = element_text(size = 10,angle = 90,vjust = .5),
        axis.text.y = element_text(size = 10))+
  ylab('Proportion')+
  xlab('Condition')
dev.off()

write.csv(fisher,'fishers_newtypes_proportion.csv',quote = F,row.names = T)



### Check
stmerged <- readRDS('merged_manuscript_020824.RDS')

pred <- as.data.frame(reshape2::melt(stmerged@assays$pred.new.types@data))
pred <- pred[pred$Var1 != 'max',]
pred$celltype <- pred$Var1
pred$celltype <- as.character(pred$celltype)
pred[!pred$celltype %in% coltable$New.types,"celltype"] <- 'Other'
pred$gloms <- stmerged@meta.data[pred$Var2,'Gloms']

fisher <- as.data.frame(matrix(0,ncol=2,nrow = 9,
                               dimnames = list(coltable$New.types,c('p.value','odds.raio'))))
forest <- as.data.frame(matrix(0,ncol=3,nrow = 9,
                               dimnames = list(coltable$New.types,c('low','est','high'))))
for (ct in coltable$New.types){
  fishmat <- matrix(c(sum(pred[pred$celltype==ct & pred$gloms == 'glom',"value"]),
                      sum(pred[pred$celltype==ct & pred$gloms == 'no_glom',"value"]),
                      sum(pred[pred$celltype!=ct & pred$gloms == 'glom',"value"]),
                      sum(pred[pred$celltype!=ct & pred$gloms == 'no_glom',"value"])),ncol=2)
  test <- fisher.test(fishmat)
  fisher[ct,"p.value"] <- test$p.value
  fisher[ct,"odds.raio"] <- test$estimate
  forest[ct,'est'] <- log(test$estimate)
  forest[ct,'low'] <- log(test$estimate) - 1.96*sqrt(1/fishmat[1,1] + 1/fishmat[2,1] + 
                                                       1/fishmat[1,2] + 1/fishmat[2,2])
  forest[ct,'high'] <- log(test$estimate) + 1.96*sqrt(1/fishmat[1,1] + 1/fishmat[2,1] + 
                                                        1/fishmat[1,2] + 1/fishmat[2,2])
}

forest$cell <- factor(rownames(forest),levels=rev(coltable$New.types))
pdf('forest_region_proportions.pdf',width = 4,height = 3)
ggplot(data=forest, aes(x=cell, y=est, ymin=low, ymax=high,color=cell)) +
  geom_pointrange(size = .7,linewidth = 1,shape = 2) + 
  geom_hline(yintercept=0, lty=2) +  # add a dotted line at x=1 after flip
  coord_flip() +  # flip coordinates (puts labels on y axis)
  xlab("Cell type") + ylab("Ln(OR) (95% CI)")+
  scale_x_discrete(limits=rev(coltable$New.types[c(1,3:8)]))+
  scale_color_manual(values = rev(coltable$Colors[c(1,3:8)]),
                     name = element_blank(),
                     guide = guide_legend(reverse = TRUE),
                     limits=rev(coltable$New.types[c(1,3:8)]))+
  theme_classic()+
  theme(panel.border = element_rect(color = "black", fill=NA, linewidth = 1.5),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10))
dev.off()


pred <- aggregate(pred$value,by=list(spot=pred$Var2,celltype=pred$celltype),FUN=sum)
pred$condition <- stmerged@meta.data[pred$spot,'Gloms']
pred <- aggregate(pred$x,by=list(celltype=pred$celltype,Condition=pred$condition),FUN=mean)
pred <- pred[pred$celltype %in% coltable$New.types[1:8],]
pred$celltype <- factor(pred$celltype,level=c(coltable$New.types))
pred$Condition <- factor(pred$Condition,level=c('glom','no_glom'))

pdf('barplot_region_proportions.pdf',width = 6,height = 3)
ggplot(pred,aes(x=celltype,y=x,fill=Condition))+
  geom_bar(stat = 'identity',position = 'dodge')+
  scale_fill_manual(values = c('#AC3CE6','#CFE679'))+
  theme_classic()+
  theme(panel.border = element_rect(color = "black", fill=NA, linewidth = 1.5),
        axis.text.x = element_text(size = 10,angle = 90,vjust = .5),
        axis.text.y = element_text(size = 10))+
  ylab('Proportion')+
  xlab('Condition')
dev.off()

write.csv(fisher,'fishers_region_proportion.csv',quote = F,row.names = T)


ecs <- c("EC-DVR","EC-AVR","EC-AEA","EC-PTC","EC-GC","dEC-PTC","cycEC","EC-LYM","dEC")
pred <- as.data.frame(reshape2::melt(stmerged@assays$predsubclassl2@data))
pred <- pred[pred$Var1 != 'max',]
pred$celltype <- pred$Var1
pred$celltype <- as.character(pred$celltype)
pred <- pred[pred$celltype %in% ecs,] 
pred$region <- stmerged@meta.data[pred$Var2,'Gloms']
pred <- aggregate(pred$value,by=list(celltype=pred$celltype,region=pred$region),FUN=sum)
pred$proportion <- NA
for (i in 1:nrow(pred)){
  pred[i,'proportion'] <- pred[i,'x'] / sum(pred[pred$region == pred[i,'region'],'x'])
}
prop <- reshape2::dcast(pred,celltype~region,value.var = 'proportion')
prop[6,2:3]/prop[5,2:3]
