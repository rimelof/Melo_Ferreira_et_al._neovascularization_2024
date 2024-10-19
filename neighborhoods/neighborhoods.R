library(Seurat)
library(ggplot2)
library(patchwork)
library(cluster)

coltable <- read.csv('../identify_glom_cell_types/colors_new_types.csv')
coltable <- rbind(coltable,
                  list(New.types = 'Other',Colors = '#DDDDDD'))
rownames(coltable) <- coltable[,1]


stmerged <- readRDS('../create_glom_ST/gloms_manuscript_020824.RDS')

SCT <- stmerged[['SCT']]
Spatial <- stmerged[['Spatial']]
predsub2 <- stmerged[['predsubclassl2']]
predsub1 <- stmerged[['predsubclassl1']]
DefaultAssay(stmerged) <- 'pred.new.types'
stmerged <- subset(stmerged,features = rownames(stmerged@assays$pred.new.types)[1:(nrow(stmerged@assays$pred.new.types)-1)]) 
VariableFeatures(stmerged) <- rownames(stmerged@assays$pred.new.types)
stmerged[['SCT']] <- SCT
stmerged[['Spatial']] <- Spatial
stmerged[['predsubclassl2']] <- predsub2
stmerged[['predsubclassl1']] <- predsub1
stmerged  <- ScaleData(stmerged)
stmerged <- RunPCA(stmerged, verbose = FALSE) #should I be setting approx = FALSE to run a full PCA instead
stmerged <- FindNeighbors(stmerged, dims = 1:50)
stmerged <- FindClusters(stmerged, verbose = FALSE,resolution = 2)
stmerged <- RunUMAP(stmerged, dims = 1:50)

pdf('umap_res2.pdf',height = 4,width = 6)
DimPlot(stmerged,label = T)
dev.off()
DimPlot(stmerged,group.by = 'orig.ident')
pdf('umap_cond_res2.pdf',height = 4,width = 5)
DimPlot(stmerged,group.by = 'condition')
dev.off()

saveRDS(stmerged,'gloms_neighborhood_res2_020524.RDS')
#stmerged <- readRDS('gloms_neighborhood_res2_020524.RDS')

stmerged@meta.data$cond_glom <- paste0(stmerged@meta.data$condition,stmerged@meta.data$glom)
cluster_tab <- table(stmerged@meta.data$cond_glom, stmerged@meta.data$seurat_clusters) # creates a table by condition and seurat cluster
#cluster_tab <- rbind(cluster_tab, colSums(cluster_tab[c(1,3),]))
sum_tab <- data.frame(rowSums(cluster_tab)) #table of all the conditions summed for all the clusters

fisher_results <- data.frame(matrix(ncol=2, nrow = 21), row.names = 0:20) #empty df to hold the results
colnames(fisher_results) <- c('p.value', 'odds.ratio') #name of the columns
fish_tab <- data.frame(matrix(ncol = 2, nrow = 2), row.names = c('glom_SLE', 'glom_control')) #table to hold the fisher values when iterating
colnames(fish_tab) <- c('cluster_x', 'other') 

cluster_tab

forest <- as.data.frame(matrix(0,ncol=3,nrow = 21,
                               dimnames = list(0:20,c('low','est','high'))))

for(i in  1:21){
  fish_tab$cluster_x[1] <- cluster_tab[1, i] # assigning value to freq of Cluster x in SLE
  fish_tab$cluster_x[2] <- cluster_tab[2, i] #assigning value to freq of cluster x in control
  fish_tab$other[1] <- (sum_tab[1, 1]- cluster_tab[1, i]) #assigning value  to freq of other clusters in SLE
  fish_tab$other[2] <- (sum_tab[2,1] - cluster_tab[2, i]) #assigning value to freq of other clusters in control
  
  test <- fisher.test(fish_tab)
  
  fisher_results$p.value[i] <-test$p.value #assign p-val to the fisher result table
  fisher_results$odds.ratio[i] <-data.frame(test$estimate)[1,1] #assign odd ratio to the fisher result table
  forest[i,'est'] <- test$estimate#log(test$estimate)
  forest[i,'low'] <- exp(log(test$estimate) - 1.96*sqrt(1/fish_tab[1,1] + 1/fish_tab[2,1] + 
                                                           1/fish_tab[1,2] + 1/fish_tab[2,2]))
  forest[i,'high'] <- exp(log(test$estimate) + 1.96*sqrt(1/fish_tab[1,1] + 1/fish_tab[2,1] + 
                                                            1/fish_tab[1,2] + 1/fish_tab[2,2]))
  
}

write.csv(fisher_results,'fishers_res2.csv')
#fisher_results <- read.csv('fishers_res2.csv',row.names = 1)
write.csv(forest,'forest_res2.csv')
#forest <- read.csv('forest_res2.csv',row.names = 1)

###########
# Barplot #
###########

coltable <- rbind(coltable,
                  data.frame(list(New.types=c('FIB','IMM'),Colors=c('#DBDA59','#5EDB3A'))))
rownames(coltable) <- coltable$New.types
coltable <- coltable[c(1,3:8,10,11,9),]

pred <- as.data.frame(reshape2::melt(stmerged@assays$pred.new.types@data))
pred$celltype <- pred$Var1
pred$celltype <- as.character(pred$celltype)
pred[pred$Var1 %in% c('MYOF','cycMYOF','FIB','M-FIB','dM-FIB','aFIB','dFIB'),'celltype'] <- 'FIB'
pred[pred$Var1 %in% c('B','PL','T','NKT','MAST','MAC-M2','cycMNP','MDC','cDC','pDC','ncMON','N'),'celltype'] <- 'IMM'
pred[!pred$celltype %in% coltable$New.types,"celltype"] <- 'Other'
pred <- aggregate(pred$value,by=list(spot=pred$Var2,celltype=pred$celltype),FUN=sum)
pred$cluster <- stmerged@meta.data[pred$spot,'seurat_clusters']
pred <- aggregate(pred$x,by=list(celltype=pred$celltype,cluster=pred$cluster),FUN=mean)
pred$celltype <- factor(pred$celltype,level=c(coltable$New.types))

fisher_results <- fisher_results[order(fisher_results$odds.ratio,decreasing = T),]
pred$cluster <- factor(pred$cluster,levels = as.numeric(rownames(fisher_results)))


pdf('barplot_all_res2.pdf',width = 8,height = 5)
ggplot(pred,aes(x=cluster,y=x,fill=celltype))+
  geom_bar(stat = 'identity',position = 'stack')+
  scale_fill_manual(values = coltable$Colors)+
  theme_classic()
dev.off()

sigclus <- rownames(fisher_results[abs(fisher_results$p.value) < .01,])
pdf('barplot_signif_res2.pdf',width = 5,height = 5)
ggplot(pred,aes(x=cluster,y=x,fill=celltype))+
  geom_bar(stat = 'identity',position = 'stack')+
  scale_fill_manual(values = coltable$Colors)+
  scale_x_discrete(limits=sigclus)+
  theme_classic()
dev.off()


fisher_results$cluster <- factor(rownames(fisher_results),levels = rownames(fisher_results))
fisher_results$condition <- ifelse(fisher_results$odds.ratio >1,'DKD','Ref')
fisher_results$log2odds <- log2(fisher_results$odds.ratio)
fisher_results['19',"log2odds"] <- fisher_results['18',"log2odds"] - .4
ggplot(fisher_results,aes(x=cluster,y=odds.ratio-1))+
  geom_bar(stat = "identity")
pdf('barods_all_res2.pdf',width = 8,height = 2)
ggplot(fisher_results,aes(x=cluster,y=log2odds,fill=condition))+
  geom_bar(stat = "identity")+
  ylim(-3.9,3)+
  scale_fill_manual(values = c('#FFA829','#8DCFFF'))+
  theme_classic()
dev.off()
pdf('barods_signif_res2.pdf',width = 4,height = 2)
ggplot(fisher_results,aes(x=cluster,y=log2odds,fill=condition))+
  geom_bar(stat = "identity")+
  ylim(-3.9,3)+
  scale_fill_manual(values = c('#FFA829','#8DCFFF'))+
  scale_x_discrete(limits=sigclus)+
  theme_classic()
dev.off()

forest$cell <- factor(rownames(forest),levels=rev(rownames(fisher_results)))
pdf('forest_all_res2.pdf',width = 6,height = 2)
ggplot(data=forest, aes(x=cell, y=est, ymin=low, ymax=high,color=cell)) +
  geom_pointrange(size = .7,linewidth = 1,shape = 2) + 
  geom_hline(yintercept=1, lty=2) +  # add a dotted line at x=1 after flip
  #coord_flip() +  # flip coordinates (puts labels on y axis)
  scale_x_discrete(limits=rownames(fisher_results))+
  xlab("Cell type") + ylab("OR (95% CI)")+
  scale_color_manual(values = c(rep('#8DCFFF',10),rep('#dddddd',6),rep('#FFA829',5)),
                     name = element_blank(),
                     guide = guide_legend(reverse = TRUE))+
  theme_classic()+
  theme(panel.border = element_rect(color = "black", fill=NA, linewidth = 1.5),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10))
dev.off()

forest$cell <- factor(rownames(forest),levels=rev(rownames(fisher_results)))
pdf('forest_signif_res2.pdf',width = 5,height = 2)
ggplot(data=forest, aes(x=cell, y=est, ymin=low, ymax=high,color=cell)) +
  geom_pointrange(size = .7,linewidth = 1,shape = 2) + 
  geom_hline(yintercept=1, lty=2) +  # add a dotted line at x=1 after flip
#  coord_flip() +  # flip coordinates (puts labels on y axis)
  xlab("Cell type") + ylab("OR (95% CI)")+
  scale_color_manual(values = c(rep('#8DCFFF',10),rep('#dddddd',6),rep('#FFA829',5)),
                     name = element_blank(),
                     guide = guide_legend(reverse = TRUE))+
#  scale_x_discrete(limits=rev(sigclus))+
  scale_x_discrete(limits=sigclus)+
  theme_classic()+
  theme(panel.border = element_rect(color = "black", fill=NA, linewidth = 1.5),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10))
dev.off()

#########################
# Explore neighborhoods #
#########################

#stmerged <- readRDS('gloms_neighborhood_res2_020524.RDS')

DefaultAssay(stmerged) <- 'SCT'
markers_1vs2018 <- FindMarkers(stmerged,ident.1 = 1,ident.2 = c(20,18))
markers_1vs2018 <- markers_1vs2018[markers_1vs2018$p_val_adj < 0.01,]
markers_1vs2018$genes <- rownames(markers_1vs2018)
markers_1vs2018 <- markers_1vs2018[order(markers_1vs2018$avg_log2FC,decreasing = T),]

write.csv(markers_1vs2018,'neigborhood_1vs2018_dkd_spatial_markers.csv',quote = F,row.names = F)

markers_8vs2018 <- FindMarkers(stmerged,ident.1 = 8,ident.2 = c(20,18))
markers_8vs2018 <- markers_8vs2018[markers_8vs2018$p_val_adj < 0.01,]
markers_8vs2018$genes <- rownames(markers_8vs2018)
markers_8vs2018 <- markers_8vs2018[order(markers_8vs2018$avg_log2FC,decreasing = T),]

write.csv(markers_8vs2018,'neigborhood_8vs2018_dkd_spatial_markers.csv',quote = F,row.names = F)

markers_1511vs31419 <- FindMarkers(stmerged,ident.1 = c(15,11),ident.2 = c(3,14,19))
markers_1511vs31419 <- markers_1511vs31419[markers_1511vs31419$p_val_adj < 0.01,]
markers_1511vs31419$genes <- rownames(markers_1511vs31419)
markers_1511vs31419 <- markers_1511vs31419[order(markers_1511vs31419$avg_log2FC,decreasing = T),]

write.csv(markers_1511vs31419,'neigborhood_1511vs31419_dkd_spatial_markers.csv',quote = F,row.names = F)



ecgenes <- rev(c('IGFBP7','AEBP1','TCIM','APLNR','FBLN5',
                'MMP2','MT1G','FRZB','IGLC2','IGLC3',
                'CCN1','SOST','DEFB1',
                'CLIC5','SPP1','PTPRM','DDIT4','CDKN1C','SPOCK2','SPOCK1','TYRO3',
                'NEDD9','KDR','CEACAM1','NTN4','FABP1','C1QTNF1','NKD1','AGT','C1QL1','SOCS3','PKHD1L1'))
podgenes <- rev(c('APLNR','ARHGEF15','SOST','TBX3','TEK','IGFBP5','EGFL7','POSTN','ETS1','FLT1','PTPRB','TGFBR2',
                  'COX7C','NDUFS5','MIF','CRYAB'))
iecgenes <- c()


# DotPlot(stmerged,idents = sigclus,features = angiogenes)+
#   theme(axis.text.x = element_text(angle = 90,vjust = .5))
pdf('dotplot_ec_neighborhoods.pdf',width=4,height=8)
DotPlot(stmerged,idents = c('1','20','18'),features = ecgenes)+
#  theme(axis.text.x = element_text(angle = 90,vjust = .5))+
  coord_flip()+
  scale_y_discrete(limits = c('1','20','18'))
dev.off()

DotPlot(stmerged,idents = c('15','11','3','14','19'),features = podgenes)+
  #  theme(axis.text.x = element_text(angle = 90,vjust = .5))+
  coord_flip()+
  scale_y_discrete(limits = c('15','11','3','14','19'))

table(stmerged@meta.data[,c("orig.ident","seurat_clusters")])
imgs <- stringr::str_replace_all(c('V10S14-087_XY01_21-0061'),'-','.')
SpatialDimPlot(stmerged,images = imgs,crop = F)+
  scale_fill_manual(values=randomcoloR::distinctColorPalette(21))
