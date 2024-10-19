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
stmerged <- FindClusters(stmerged, verbose = FALSE,resolution = c(0.4,0.6,0.8,1,1.2,1.4,1.6,1.8,2,2.2,2.4,2.6,2.8,3))
stmerged <- RunUMAP(stmerged, dims = 1:50)

pdf('featureplot_test.pdf',height = 12)
FeaturePlot(stmerged,c('POD',
                              'EC-GC','I-EC','EC-NOS',
                              'VSMC-NOS','MC','I-VSMC'),ncol = 2)
dev.off()
pdf('dotplot_test.pdf')
for(res in c(0.4,0.6,0.8,1,1.2,1.4,1.6,1.8,2,2.2,2.4,2.6,2.8,3)){
  plot(DotPlot(stmerged,c('POD',
                   'EC-GC','I-EC','EC-NOS',
                   'VSMC-NOS','MC','I-VSMC',
                   'B', 'PL', 'T', 'N', 'ncMON', "MDC", 'MAST',
                   'FIB', 'M-FIB', 'MYOF', "dFIB" ),
        group.by = paste0('pred.new.types_snn_res.',res))+
         theme(axis.text.x = element_text(angle=90,vjust = .5))+
         ggtitle(paste0('res = ',res)))
}
dev.off()

pdf('dimplot_test.pdf',width = 12)
for(res in c(0.4,0.6,0.8,1,1.2,1.4,1.6,1.8,2,2.2,2.4,2.6,2.8,3)){
  a <- DimPlot(stmerged,label = T,
               group.by = paste0('pred.new.types_snn_res.',res))
  b <- DimPlot(stmerged,label = T,
               group.by = 'condition')
  plot(a + b)
}
dev.off()

df_sil <- NULL
dist.matrix <- dist(x = Embeddings(object = stmerged[['pca']])[, 1:50])
for(res in c(0.4,0.6,0.8,1,1.2,1.4,1.6,1.8,2,2.2,2.4,2.6,2.8,3)){
  clusters <- stmerged@meta.data[,paste0('pred.new.types_snn_res.',res)]
  sil <- silhouette(x = as.numeric(x = as.factor(x = clusters)), dist = dist.matrix)
  df_sil <- cbind(df_sil,sil[,3])
  colnames(df_sil)[ncol(df_sil)] <- paste0('res_',res)
}
df_sil <- reshape2::melt(df_sil)
ggplot(df_sil,aes(x=Var2,y=value))+
  geom_boxplot()

#Fishers Exact Test
stmerged@meta.data$cond_glom <- paste0(stmerged@meta.data$condition,stmerged@meta.data$glom)
cluster_tab <- table(stmerged@meta.data$cond_glom, stmerged@meta.data$pred.new.types_snn_res.1.2) # creates a table by condition and seurat cluster
#cluster_tab <- rbind(cluster_tab, colSums(cluster_tab[c(1,3),]))
sum_tab <- data.frame(rowSums(cluster_tab)) #table of all the conditions summed for all the clusters

fisher_results <- data.frame(matrix(ncol=2, nrow = 13), row.names = 0:12) #empty df to hold the results
colnames(fisher_results) <- c('p.value', 'odds.ratio') #name of the columns
fish_tab <- data.frame(matrix(ncol = 2, nrow = 2), row.names = c('glom_SLE', 'glom_control')) #table to hold the fisher values when iterating
colnames(fish_tab) <- c('cluster_x', 'other') 

cluster_tab

for(i in  1:13){
  fish_tab$cluster_x[1] <- cluster_tab[1, i] # assigning value to freq of Cluster x in SLE
  fish_tab$cluster_x[2] <- cluster_tab[2, i] #assigning value to freq of cluster x in control
  fish_tab$other[1] <- (sum_tab[1, 1]- cluster_tab[1, i]) #assigning value  to freq of other clusters in SLE
  fish_tab$other[2] <- (sum_tab[2,1] - cluster_tab[2, i]) #assigning value to freq of other clusters in control
  
  test <- fisher.test(fish_tab)
  
  fisher_results$p.value[i] <-test$p.value #assign p-val to the fisher result table
  fisher_results$odds.ratio[i] <-data.frame(test$estimate)[1,1] #assign odd ratio to the fisher result table
}


#pdf('umap_73_res3.pdf',height = 4,width = 6)
DimPlot(stmerged,label = T)
#dev.off()
DimPlot(stmerged,group.by = 'orig.ident')
#pdf('umap_73_cond_res3.pdf',height = 4,width = 5)
DimPlot(stmerged,group.by = 'condition')
#dev.off()


