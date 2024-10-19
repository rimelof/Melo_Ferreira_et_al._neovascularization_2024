library(Seurat)
library(rio)
library(ggplot2)
library(ggpubr)

ffpe <- readRDS('integrated_ffpe_04242024.RDS')
write.table(unique(ffpe@meta.data[,c("orig.ident")]),
            '../samples_list/ffpe_samples.csv',
            quote = F,row.names = F,col.names = F,sep = ',')
Idents(ffpe) <- ffpe$Gloms

gloms <- subset(ffpe,idents = 'no_glom',invert=T)
gloms <- FindVariableFeatures(gloms)
#gloms <- SCTransform(gloms,assay = 'spatial')
gloms <- ScaleData(gloms)
gloms <- RunPCA(gloms)
gloms <- RunUMAP(gloms,dims = 1:30)
DimPlot(gloms,group.by = 'orig.ident',label = T)

DefaultAssay(gloms) <- 'SCT'


ST <- readRDS("../neighborhoods/gloms_neighborhood_res2_020524.RDS")
ST <- UpdateSeuratObject(ST)
DefaultAssay(ST) <- 'SCT'
Idents(ST) <- 'seurat_clusters' # new

DefaultAssay(ST)
DefaultAssay(gloms)
# find anchors
gloms <- PrepSCTFindMarkers(gloms)
anchors <- FindTransferAnchors(reference = ST, query = gloms)
# transfer labels
new_types <- TransferData(
  anchorset = anchors,
  refdata = ST$seurat_clusters,prediction.assay = T)
gloms[['pred.neighborhood']] <- new_types

new_types <- TransferData(
  anchorset = anchors,
  refdata = ST$seurat_clusters)
gloms <- AddMetaData(object = gloms, metadata = new_types)
gloms@meta.data$predneighborhood <- gloms$predicted.id

saveRDS(gloms,'gloms_ffpe.RDS')
#gloms <- readRDS('gloms_ffpe.RDS')

DimPlot(gloms,group.by = 'predneighborhood',label = T)


scoring <- import('scoring.xlsx')
scoring$glom_index <- paste0(stringr::str_split(scoring$Sample,'_',simplify = T)[,3],
                             '_G',scoring$Number)
scoring[scoring == 'na'] <- NA
colnames(scoring)[4:8] <- c('Expansion','KW','Obsolescence','FSGS','drop')

gloms$glom_index <- paste0(gloms$orig.ident,'_',gloms$Gloms)
df <- as.data.frame(t(gloms@assays$pred.neighborhood@data))
colnames(df) <- paste0('nb_',colnames(df))
df$glom_index <- gloms@meta.data[rownames(df),"glom_index"]

fulldf <- merge.data.frame(df,scoring,by='glom_index',all=F)
plotdf <- reshape2::melt(fulldf[,c("nb_1",'nb_18','Neovasc',"glom_index")])
plotdf$variable <- factor(plotdf$variable,levels=c("nb_18",'nb_1'))
plotdf <- plotdf[!is.na(plotdf$Neovasc),]
pdf('violin_neaovasc.pdf',width = 2.5,height = 3)
ggplot(plotdf,aes(x=Neovasc,y=value,fill=variable))+
  geom_violin(position = 'dodge',scale = 'width',draw_quantiles = c(0.5))+
  scale_fill_manual(values = c('#8DCFFF','#FFA829'))+
  theme_classic()
dev.off()
fmatrix <- matrix(c(sum(plotdf[plotdf$variable == 'nb_1' & plotdf$Neovasc == 'Y',"value"]),
                    sum(plotdf[plotdf$variable == 'nb_18' & plotdf$Neovasc == 'Y',"value"]),
                    sum(plotdf[plotdf$variable == 'nb_1' & plotdf$Neovasc == 'N',"value"]),
                    sum(plotdf[plotdf$variable == 'nb_18' & plotdf$Neovasc == 'N',"value"])),ncol = 2)
ftest <- fisher.test(fmatrix)
ftest$estimate
ftest$p.value

plotdf <- reshape2::melt(fulldf[,c("nb_6",'FSGS',"glom_index")])
#plotdf$variable <- factor(plotdf$variable,levels=c("nb_18",'nb_1'))
pdf('violin_FSGS.pdf',width = 2,height = 3)
ggplot(plotdf[!is.na(plotdf$FSGS),],aes(x=FSGS,y=value,fill=variable))+
  geom_violin(position = 'dodge',scale = 'width',draw_quantiles = c(0.5))+
  scale_fill_manual(values = c('#dbda59'))+
  theme_classic()
dev.off()

assoc <- matrix(0,ncol=6,nrow=11)
colnames(assoc) <- colnames(scoring)[c(3:8)]
rownames(assoc) <- paste0('nb_',c(1,15,2,6,11,12,3,9,14,18,19))
pval <- matrix(0,ncol=6,nrow=11)
colnames(pval) <- colnames(scoring)[c(3:8)]
rownames(pval) <- paste0('nb_',c(1,15,2,6,11,12,3,9,14,18,19))
testdf <- reshape2::melt(fulldf[,c(colnames(assoc),row.names(assoc))],)
for (rr in rownames(assoc)){
  for (cc in colnames(assoc)){
    plotdf <- testdf[!is.na(testdf[,cc]),]
    fmatrix <- matrix(c(sum(plotdf[plotdf$variable == rr & plotdf[,cc] == 'Y',"value"]),
                        sum(plotdf[plotdf$variable != rr & plotdf[,cc] == 'Y',"value"]),
                        sum(plotdf[plotdf$variable == rr & plotdf[,cc] == 'N',"value"]),
                        sum(plotdf[plotdf$variable != rr & plotdf[,cc] == 'N',"value"])),ncol = 2)
    ftest <- fisher.test(fmatrix)
    print(c(rr,cc,ftest$estimate,ftest$p.value))
    assoc[rr,cc] <- ftest$estimate
    pval[rr,cc] <- ftest$p.value
  }
}

row.names(assoc) <- c('DKD EC','DKD POD 1', 'DKD PG','DKD FIB','DKD POD 2',
                      'Ref I-EC','Ref POD 1','Ref PG','Ref POD 2','Ref EC 2', 'Ref POD 3')
row.names(pval) <- c('DKD EC','DKD POD 1', 'DKD PG','DKD FIB','DKD POD 2',
                      'Ref I-EC','Ref POD 1','Ref PG','Ref POD 2','Ref EC 2', 'Ref POD 3')
