library(Seurat)
library(sp)
library(ggplot2)

selec_files <- list.files('../glom_selections_xenium/','coordinates')
selec_files <- selec_files[grep('Selection',selec_files,invert = T)]
samples <- stringr::str_split(selec_files,'_',simplify = T)[,1]

xenium_list <- readRDS('list_xenium_noNS.RDS')
xenium_list <- xenium_list[samples]
xenium_list[['3723']]@meta.data$condition = 'Ref'
xenium_list[['3775']]@meta.data$condition = 'Ref'
xenium_list[['f59']]@meta.data$condition = 'Ref'
xenium_list[['40440']]@meta.data$condition = 'DKD'
xenium_list[['40610']]@meta.data$condition = 'DKD'
xenium_list[['40775']]@meta.data$condition = 'DKD'
xenium_list[['5582']]@meta.data$condition = 'DKD'
for (smp in names(xenium_list)[1]){
  xenium_list[[smp]]@meta.data$orig.ident <- smp
}

merged.anchors <- FindIntegrationAnchors(object.list = xenium_list,
                                         reduction = 'rpca')
merged_obj <- IntegrateData(anchorset = merged.anchors)

merged_obj <- SCTransform(merged_obj,assay = 'Xenium')
merged_obj <- ScaleData(merged_obj)
merged_obj <- RunPCA(merged_obj)
saveRDS(merged_obj,'xenium_merged_noNS.RDS')
write.table(unique(merged_obj@meta.data[,c("orig.ident","condition")]),
            '../samples_list/xenium_samples.csv',
            quote = F,row.names = F,col.names = F,sep = ',')

names(merged_obj@images) <- names(xenium_list)

merged_obj@meta.data$segment <- 'TI'
for (smp in names(merged_obj@images)){
  selections <- read.csv(paste0('../glom_selections_xenium/',smp,'_coordinates..csv'),
                         row.names = NULL,skip = 2)
  
  coords <- as.data.frame(merged_obj@images[[smp]]@boundaries$centroids@coords)
  rownames(coords) <- merged_obj@images[[smp]]@boundaries$centroids@cells
  
  cells <- c()
  for (gl in unique(selections$Selection)){
    in_selec1 <- point.in.polygon(coords$x,coords$y,
                                  selections[selections$Selection == gl,'X'],
                                  selections[selections$Selection == gl,'Y'])
    
    cells <- rownames(coords)[in_selec1==1]
    merged_obj@meta.data[cells,"segment"] <- paste0(smp,'_G',
                                                   stringr::str_split(gl,' ',simplify = T)[,2])
  }

}

Idents(merged_obj) <- merged_obj@meta.data$segment

saveRDS(merged_obj,'xenium_merged_noNS.RDS')
# merged_obj <- readRDS('xenium_merged.RDS')

merged_obj@meta.data$level1 <- merged_obj@meta.data$anchors_predictions
merged_obj@meta.data$level1 <- ifelse(merged_obj@meta.data$anchors_predictions %in% c('pEC','EC-GC','EC-NOS','I-EC'),
                       'EC',merged_obj@meta.data$level1)
merged_obj@meta.data$level1 <- ifelse(merged_obj@meta.data$anchors_predictions %in% c('MC','VSMC-NOS','I-VSMC'),
                       'VSM/P',merged_obj@meta.data$level1)
merged_obj@meta.data$level1 <- ifelse(merged_obj@meta.data$anchors_predictions %in% c('POD'),
                        'POD',merged_obj@meta.data$level1)


for (smp in unique(merged_obj@meta.data$orig.ident)){
  write.csv(data.frame(list(cell_id=rownames(merged_obj@meta.data[merged_obj$orig.ident == smp,]),
                            group=merged_obj@meta.data[merged_obj$orig.ident == smp,'level1'])),
            paste0('level1_cell_types/',smp,'_l1_anchors.csv'),row.names = F,quote = F)
}

Idents(merged_obj) <- merged_obj@meta.data$segment

gloms <- subset(merged_obj,idents = 'TI',invert=T)
ImageDimPlot(gloms,fov = 'f59')
saveRDS(gloms,'xenium_gloms_noNS.RDS')
# gloms <- readRDS('xenium_gloms.RDS')


endotypes <- c('pEC','EC-GC','I-EC')
df <- gloms@meta.data[,c('segment',"condition")]
df <- unique(df)
df$num.pEC <- 0
df$num.IEC <- 0
df$num.EC <- 0
df$num.ECGC <- 0
df$num.cells <- 0
for(i in 1:nrow(df)){
  df[i,'num.pEC'] <- nrow(gloms@meta.data[gloms$segment==df[i,"segment"] &
                                            gloms$anchors_predictions=='pEC',])
  df[i,'num.IEC'] <- nrow(gloms@meta.data[gloms$segment==df[i,"segment"] &
                                            gloms$anchors_predictions=='I-EC',])
  df[i,'num.EC'] <- nrow(gloms@meta.data[gloms$segment==df[i,"segment"] &
                                           gloms$anchors_predictions %in% endotypes,])
  df[i,'num.ECGC'] <- nrow(gloms@meta.data[gloms$segment==df[i,"segment"] &
                                           gloms$anchors_predictions=='EC-GC',])
  df[i,'num.cells'] <- nrow(gloms@meta.data[gloms$segment==df[i,"segment"],])
}


ggplot(df,aes(x=condition,y=num.pEC/num.cells))+
  geom_dotplot(binaxis = 'y',stackdir = 'center')
pdf('violin_prec_all_ratio_noNS.pdf',width = 3,height = 3)
ggplot(df,aes(x=condition,y=num.pEC/num.cells,fill=condition))+
  geom_violin(draw_quantiles = 0.5)+
  coord_flip()+
  scale_fill_manual(values = c('#FFA829','#8DCFFF'))+
  theme_classic()
dev.off()

t.test(df[df$condition=='DKD',"num.pEC"]/df[df$condition=='DKD',"num.cells"],
       df[df$condition=='Ref',"num.pEC"]/df[df$condition=='Ref',"num.cells"])

mean(df[df$condition=='DKD',"num.pEC"]/df[df$condition=='DKD',"num.cells"]) /
  mean(df[df$condition=='Ref',"num.pEC"]/df[df$condition=='Ref',"num.cells"])

sd(df[df$condition=='DKD',"num.pEC"]/df[df$condition=='DKD',"num.cells"])
sd(df[df$condition=='Ref',"num.pEC"]/df[df$condition=='Ref',"num.cells"])

ggplot(df,aes(x=condition,y=num.pEC/num.EC))+
  geom_dotplot(binaxis = 'y',stackdir = 'center')

ggplot(df,aes(x=condition,y=num.pEC/num.ECGC))+
  geom_dotplot(binaxis = 'y',stackdir = 'center')
ggplot(df,aes(x=condition,y=num.pEC/num.ECGC))+
  geom_jitter()
pdf('violin_prec_ecgc_ratio.pdf',width = 3,height = 3)
ggplot(df,aes(x=condition,y=num.pEC/num.ECGC,fill = condition))+
  geom_violin(draw_quantiles = 0.5)+
  coord_flip()+
  scale_fill_manual(values = c('#FFA829','#8DCFFF'))
dev.off()

t.test(df[df$condition=='DKD',"num.pEC"]/df[df$condition=='DKD',"num.ECGC"],
       df[df$condition=='Ref',"num.pEC"]/df[df$condition=='Ref',"num.ECGC"],)

ggplot(df,aes(x=condition,y=num.pEC/num.EC))+
  geom_violin(draw_quantiles = 0.5)

ggplot(df,aes(x=condition,y=num.IEC/num.EC))+
  geom_violin(draw_quantiles = 0.5)

ggplot(df,aes(x=condition,y=num.pEC/num.cells))+
  geom_violin()

df$ratio <- df$num.pEC / df$num.EC
(sum(df[df$condition=='DKD',"ratio"]) / nrow(df[df$condition=='DKD',]))/ 
  (sum(df[df$condition=='Ref',"ratio"]) / nrow(df[df$condition=='Ref',]))
#2.608359
t.test(df[df$condition=='DKD',"ratio"],df[df$condition=='Ref',"ratio"])

df$condition=factor(df$condition,levels = c('Ref','DKD'))
pdf('violin_ratio_pecs.pdf',width = 3.75,height = 3)
ggplot(df,aes(x=condition,y=ratio,fill=condition))+
  geom_violin(draw_quantiles = 0.5)+
  scale_fill_manual(values = c('#8DCFFF','#FFA829'))+
  theme_classic()
dev.off()

DotPlot(gloms,c('SEMA6A','PLXNA2'),group.by = 'anchors_predictions',idents = endotypes)


coltable <- read.csv('../identify_glom_cell_types/colors_new_types.csv')
coltable <- rbind(coltable,
                  list(New.types = 'Other',Colors = '#DDDDDD'))
coltable[3,1] <- 'VSM/P'
coltable[6,1] <- 'EC'
coltable <- coltable[c(1,3,6,9),]
rownames(coltable) <- coltable[,1]


fisher <- as.data.frame(matrix(0,ncol=2,nrow = 3,
                               dimnames = list(coltable$New.types[c(1,3,2)],c('p.value','odds.raio'))))
forest <- as.data.frame(matrix(0,ncol=3,nrow = 3,
                               dimnames = list(coltable$New.types[c(1,3,2)],c('low','est','high'))))
for (ct in coltable$New.types[c(1,3,2)]){
  fishmat <- matrix(c(nrow(gloms@meta.data[gloms$level1==ct & gloms$condition == 'DKD',]),
                      nrow(gloms@meta.data[gloms$level1==ct & gloms$condition == 'Ref',]),
                      nrow(gloms@meta.data[gloms$level1!=ct & gloms$condition == 'DKD',]),
                      nrow(gloms@meta.data[gloms$level1!=ct & gloms$condition == 'Ref',])),ncol=2)
  test <- fisher.test(fishmat)
  fisher[ct,"p.value"] <- test$p.value
  fisher[ct,"odds.raio"] <- test$estimate
  forest[ct,'est'] <- test$estimate#log(test$estimate)
  forest[ct,'low'] <- exp(log(test$estimate) - 1.96*sqrt(1/fishmat[1,1] + 1/fishmat[2,1] + 
                                                       1/fishmat[1,2] + 1/fishmat[2,2]))
  forest[ct,'high'] <- exp(log(test$estimate) + 1.96*sqrt(1/fishmat[1,1] + 1/fishmat[2,1] + 
                                                        1/fishmat[1,2] + 1/fishmat[2,2]))
}

forest$cell <- factor(rownames(forest),levels=c('EC','VSM/P','POD'))
pdf('forest_l1_proportions_linear.pdf',width = 4,height = 3)
ggplot(data=forest, aes(x=cell, y=est, ymin=low, ymax=high,color=cell)) +
  geom_pointrange(size = .7,linewidth = 1,shape = 2) + 
  geom_hline(yintercept=1, lty=2) +  # add a dotted line at x=1 after flip
  coord_flip() +  # flip coordinates (puts labels on y axis)
  xlab("Cell type") + ylab("OR (95% CI)")+
  scale_color_manual(values = c('#ff8c00','#483d8b','#db7295'),
                     name = element_blank(),
                     guide = guide_legend(reverse = TRUE))+
  theme_classic()+
  theme(panel.border = element_rect(color = "black", fill=NA, linewidth = 1.5),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10))
dev.off()


fisher <- as.data.frame(matrix(0,ncol=2,nrow = 3,
                               dimnames = list(coltable$New.types[c(1,3,2)],c('p.value','odds.raio'))))
forest <- as.data.frame(matrix(0,ncol=3,nrow = 3,
                               dimnames = list(coltable$New.types[c(1,3,2)],c('low','est','high'))))
for (ct in coltable$New.types[c(1,3,2)]){
  fishmat <- matrix(c(nrow(gloms@meta.data[gloms$level1==ct & gloms$condition == 'DKD',]),
                      nrow(gloms@meta.data[gloms$level1==ct & gloms$condition == 'Ref',]),
                      nrow(gloms@meta.data[gloms$level1!=ct & gloms$condition == 'DKD',]),
                      nrow(gloms@meta.data[gloms$level1!=ct & gloms$condition == 'Ref',])),ncol=2)
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
pdf('forest_l1_proportions_log.pdf',width = 4,height = 3)
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


coltable <- read.csv('../identify_glom_cell_types/colors_new_types.csv')
coltable <- rbind(coltable,
                  list(New.types = 'MYOF',Colors = '#dbda59'))
fisher <- as.data.frame(matrix(0,ncol=2,nrow = 9,
                               dimnames = list(coltable$New.types,c('p.value','odds.raio'))))
forest <- as.data.frame(matrix(0,ncol=3,nrow = 9,
                               dimnames = list(coltable$New.types,c('low','est','high'))))
for (ct in coltable$New.types){
  fishmat <- matrix(c(nrow(gloms@meta.data[gloms$anchors_predictions ==ct & gloms$condition == 'DKD',]),
                      nrow(gloms@meta.data[gloms$anchors_predictions==ct & gloms$condition == 'Ref',]),
                      nrow(gloms@meta.data[gloms$anchors_predictions!=ct & gloms$condition == 'DKD',]),
                      nrow(gloms@meta.data[gloms$anchors_predictions!=ct & gloms$condition == 'Ref',])),ncol=2)
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
pdf('forest_anchors_proportions_log.pdf',width = 4,height = 3)
ggplot(data=forest, aes(x=cell, y=est, ymin=low, ymax=high,color=cell)) +
  geom_pointrange(size = .7,linewidth = 1,shape = 2) + 
  geom_hline(yintercept=0, lty=2) +  # add a dotted line at x=1 after flip
  coord_flip() +  # flip coordinates (puts labels on y axis)
  xlab("Cell type") + ylab("Ln(OR) (95% CI)")+
  scale_color_manual(values = rev(coltable$Colors),
                     name = element_blank(),
                     guide = guide_legend(reverse = TRUE))+
  theme_classic()+
  theme(panel.border = element_rect(color = "black", fill=NA, linewidth = 1.5),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10))
dev.off()
