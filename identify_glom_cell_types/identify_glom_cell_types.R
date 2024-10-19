library(Seurat)
library(harmony)
library(ggplot2)
library(SeuratDisk)
#library(SeuratWrappers)

gloms <- readRDS('../neighborhoods/gloms_ref_DKD_w_glom_index.RDS')

Idents(gloms) <- gloms$condition
DefaultAssay(gloms) <- 'SCT'
dkdref_glom <- FindMarkers(gloms,ident.1 = 'DKD',ident.2 = 'Ref')

stmerged <- readRDS('../neighborhoods/merged_41_samples_12122023.RDS')

Idents(stmerged) <- stmerged$glom
DefaultAssay(stmerged) <- 'SCT'
glom_noglom <- FindMarkers(stmerged,ident.1='glom',ident.2 = 'no_glom')

mcatlas <- readRDS('../vascular_atlas.RDS')
mcatlas <- FindVariableFeatures(mcatlas)
mcatlas <- RunPCA(FindVariableFeatures(mcatlas))
mcatlas <- RunHarmony(mcatlas,group.by.vars = 'experiment')
mcatlas <- RunUMAP(mcatlas,reduction = 'harmony', dims = 1:30)

DimPlot(mcatlas,group.by = 'subclass.l2')

dkdref_glom <- dkdref_glom[dkdref_glom$avg_log2FC > 0 & dkdref_glom$p_val_adj < .01,]
dkdref_glom <- dkdref_glom[order(dkdref_glom$avg_log2FC),]
glom_noglom <- glom_noglom[glom_noglom$avg_log2FC > 0 & glom_noglom$p_val_adj < .01,]
glom_noglom <- glom_noglom[order(glom_noglom$avg_log2FC),]

mcatlas <- AddModuleScore(mcatlas,
                          list(dkdref = rownames(dkdref_glom),
                               glomti = rownames(glom_noglom)
))

DimPlot(mcatlas,group.by = 'subclass.l2')
FeaturePlot(mcatlas,c('Cluster1','Cluster2'))

FeaturePlot(mcatlas,rownames(dkdref_glom)[1:15])
FeaturePlot(mcatlas,rownames(glom_noglom)[1:15])

mcatlas <- FindNeighbors(object = mcatlas, reduction = "harmony")
mcatlas <- FindClusters(mcatlas, resolution = c(0.2, 0.4, 0.6, 0.8, 1.0, 1.2))
DimPlot(mcatlas,group.by = "RNA_snn_res.1.2",label = T)

DotPlot(mcatlas,rownames(glom_noglom)[1:200],group.by = "RNA_snn_res.1")

genes <- c('ROBO1','PIEZO2','DAAM2','PHTF2','GATA3','POSTN',#'PIP5K1B',
           'NOTCH3','PDGFRB','ITGA8',
           'TAGLN','FLNA','MYL9','ACTA2','DSTN')
DotPlot(mcatlas,genes,group.by = "RNA_snn_res.0.8")
FeaturePlot(mcatlas,genes)


all_markers <- c('NOTCH3','PDGFRB','ITGA8',# All VSMC
                 'PIP5K1B','ROBO1','PIEZO2','DAAM2','PHTF2','GATA3','POSTN',# MC
                 'REN','PDE10A','ABCC8','COL13A1','GRID2',#'PIP5K1B','ROBO1',# REN
                 'NTRK3','MYH11','RGS6','ADRA1A','LDB3','MCAM',# VSMC
                 'RGS5','ABCC9','ADCY3','ADGRB3',#'NTRK3','CCDC102B',# VSMC/P
                 'TAGLN','FLNA','MYL9','ACTA2','DSTN')# dVSMC

Idents(mcatlas) <- mcatlas$RNA_snn_res.1.2
DimPlot(mcatlas,label = T)
DotPlot(mcatlas,all_markers)+
  theme(axis.text.x = element_text(angle = 90,vjust = .5))+
  scale_y_discrete(limits=factor(c(6,7,8,1,4,5,9,0,2,3,11,10)))

FeaturePlot(mcatlas,c('ROBO1','TAGLN'),blend = T,order = T)

DimPlot(mcatlas,label = T,split.by = "condition.l2")


markers <- list()
for (res in as.character(c(0.2, 0.4, 0.6, 0.8, 1.0, 1.2))){
  Idents(mcatlas) <- paste0('RNA_snn_res.',res)
  markers[[res]] <- FindAllMarkers(mcatlas)
}

# mk67v0 <- FindMarkers(mcatlas,ident.1 = 0,ident.2 = c(6,7),
#                       logfc.threshold = 0,min.pct = 0)
# mk67v0$genes <- rownames(mk67v0)
# 
# write.csv(mk67v0,'markers_iVSMC_vs_MC.csv',row.names = F,quote = F)


write.csv(mcatlas@meta.data,'mcatlas_cluster_inglom_metadata_dec23.csv')

mcatlas@meta.data$new_types <- ifelse(mcatlas@meta.data$RNA_snn_res.1.2 %in% c(6,7),
                                      'MC','VSMC-NOS')
mcatlas@meta.data$new_types <- ifelse(mcatlas@meta.data$RNA_snn_res.1.2 %in% c(0),
                                      'I-VSMC',mcatlas@meta.data$new_types)
DimPlot(mcatlas,label = T,group.by = "new_types")

Idents(mcatlas) <- mcatlas@meta.data$new_types
mk67v0 <- FindMarkers(mcatlas,ident.1 = 'I-VSMC',ident.2 = 'MC',
                      logfc.threshold = 0,min.pct = 0)
mk67v0$genes <- rownames(mk67v0)
ivsmc_cells <- WhichCells(mcatlas, idents = "I-VSMC")
mc_cells <- WhichCells(mcatlas, idents = "MC")
avg_exp <- data.frame(list(
  avg_ivsmc = apply(mcatlas@assays$RNA@data[,ivsmc_cells], 1, mean),
  std_ivsmc = apply(mcatlas@assays$RNA@data[,ivsmc_cells], 1, sd),
  avg_mc = apply(mcatlas@assays$RNA@data[,mc_cells], 1, mean),
  std_mc = apply(mcatlas@assays$RNA@data[,mc_cells], 1, sd)
),row.names = rownames(mcatlas@assays$RNA@data))
avg_exp$signal2noise <- (avg_exp$avg_ivsmc - avg_exp$avg_mc) / 
  (avg_exp$std_ivsmc + avg_exp$std_mc)

mk67v0$signal2noise <- avg_exp[rownames(mk67v0),"signal2noise"] 

write.csv(mk67v0,'markers_iVSMC_vs_MC.csv',row.names = F,quote = F)


saveRDS(mcatlas,'mcatlas_reclustered.RDS')
#mcatlas <- readRDS('mcatlas_reclustered.RDS')

ecatlas <- readRDS('../endothelial_atlas.RDS')
ecatlas <- FindVariableFeatures(ecatlas)
ecatlas <- RunPCA(FindVariableFeatures(ecatlas))
ecatlas <- RunHarmony(ecatlas,group.by.vars = 'experiment')
ecatlas <- RunUMAP(ecatlas,reduction = 'harmony', dims = 1:30)

DimPlot(ecatlas,group.by = 'subclass.l2')

ecatlas <- FindNeighbors(object = ecatlas, reduction = "harmony")
ecatlas <- FindClusters(ecatlas, resolution = c(0.2, 0.4, 0.6, 0.8, 1.0, 1.2))
DimPlot(ecatlas,group.by = "RNA_snn_res.1.2",label = T)

write.csv(ecatlas@meta.data,'ecatlas_cluster_inglom_metadata_dec23.csv')

saveRDS(ecatlas,'ecatlas_reclustered.RDS')
#ecatlas <- readRDS('ecatlas_reclustered.RDS')
stecgenes <- rev(c('IGFBP7','RPS10','AEBP1','TCIM','APLNR','CXCL12','FBLN5',
                 'TAGLN','POSTN','LUM','MMP2','MT1G','ERAP2','FRZB','IGLC2','IGLC3',
                 'SOST','CLIC5','SPP1','PTPRM','DDIT4','CDKN1C','SPOCK2','SPOCK1','TYRO3',
                 'NEDD9','KDR','CEACAM1','NTN4','FABP1','C1QTNF1','NKD1','AGT','C1QL1','SOCS3','PKHD1L1'))
DotPlot(ecatlas,stecgenes)+
  coord_flip()
DimPlot(ecatlas,label = T)
prolifGC <- FindMarkers(ecatlas,ident.1 = c('1','2','8'),ident2=c('4','6','14'))
prolifmarkers <- FindMarkers(ecatlas,ident.1 = c(1,2,8))
testmarkes <- FindMarkers(ecatlas,ident.1 = c('1','2','8'),ident2=c('12','9'))
prolifmarkers2 <- FindMarkers(ecatlas,ident.1 = c('1'))
allmarkers <- FindAllMarkers(ecatlas)
nrow(allmarkers[allmarkers$cluster == '1',])
prolifGC2 <- FindMarkers(ecatlas,ident.1 = c('1'),ident2=c('4','6','14'))

pdf('feature_ec_st_genes.pdf')
for (g in stecgenes){
  plot(FeaturePlot(ecatlas,g,order = T))
}
dev.off()

all_markers <- c('CD34','PECAM1','PTPRB','MEIS2','FLT1','EMCN',# EC
                 'HECW2','PLAT','ITGA8','EHD3','KDR','SOST',# 'EMCN',# EC-GC
                 'BTNL9','ADAMTS6','PALMD','AQP1','TM4SF1',
                 'VEGFC','CCDC3','CDH5','SERPINE2','FBLN5','CXCL12','SOX17',# EC-AEA
                 'MCTP1','SLC14A1','ENPP2','LYPD6B',#'BTNL9','ADAMTS6','PALMD','AQP1','TM4SF1', # EC-DVR
                 'CEACAM1','DNASE1L3','PLVAP','PITPNC1','GRB10','SLCO2A1','RAPGEF4',# EC-PTC
                 'GPM6A','EDIL3','TLL1','ZNF385D','NR2F2',#'CEACAM1','DNASE1L3','PLVAP', # EC-AVR
                 'MMRN1','CD36','TBX1','PKHD1L1','PROX1',# EC-LYM
                 'B2M','TMSB4X','TMSB10','HLA-B','IGFBP5',# dEC
                 'FTL','DHFR','SPP1','GPX3','PEBP1')# dEC-PTC

Idents(ecatlas) <- ecatlas$RNA_snn_res.1.2
DimPlot(ecatlas,label = T)
DotPlot(ecatlas,all_markers)+
  theme(axis.text.x = element_text(angle = 90,vjust = .5))+
  scale_y_discrete(limits=factor(c(4,6,14,9,12,7,0,2,3,11,13,1,10,5,8)))

DimPlot(ecatlas,label = T,group.by = 'subclass.l2')

markers <- list()
for (res in as.character(c(0.2, 0.4, 0.6, 0.8, 1.0, 1.2))){
  Idents(ecatlas) <- paste0('RNA_snn_res.',res)
  markers[[res]] <- FindAllMarkers(ecatlas)
}

saveRDS(markers,'markers_ec_resolutions.RDS')

# mk_5.8.10_4.6.14 <- FindMarkers(ecatlas, ident.1 = c(5,8,10),
#                                 ident.2 = c(4,6,14))
# mk_5.8.10_4.6.14$genes <- rownames(mk_5.8.10_4.6.14)
# write.csv(mk_5.8.10_4.6.14,'markers_iEC_vs_EC-GC.csv',row.names = F,quote = F)

ecatlas@meta.data$new_types <- ifelse(ecatlas@meta.data$RNA_snn_res.1.2 %in% c(4,6,14),
                                      'EC-GC','EC-NOS')
ecatlas@meta.data$new_types <- ifelse(ecatlas@meta.data$RNA_snn_res.1.2 %in% c(5,8,10),
                                      'I-EC',ecatlas@meta.data$new_types)
DimPlot(ecatlas,label = T,group.by = "new_types")

ecatlas@meta.data$new_prolif <- ifelse(ecatlas@meta.data$RNA_snn_res.1.2 %in% c(4,6,14),
                                      'EC-GC',ecatlas@meta.data$RNA_snn_res.1.2)
ecatlas@meta.data$new_prolif <- ifelse(ecatlas@meta.data$RNA_snn_res.1.2 %in% c(5,8,10),
                                      'I-EC',ecatlas@meta.data$new_prolif)
DimPlot(ecatlas,label = T,group.by = "new_prolif")


Idents(ecatlas) <- ecatlas@meta.data$new_types
mk_5.8.10_4.6.14 <- FindMarkers(ecatlas, ident.1 = 'I-EC',
                                ident.2 = 'EC-GC')
mk_5.8.10_4.6.14$genes <- rownames(mk_5.8.10_4.6.14)
iec_cells <- WhichCells(ecatlas, idents = "I-EC")
ecgc_cells <- WhichCells(ecatlas, idents = "EC-GC")
avg_exp <- data.frame(list(
  avg_iec = apply(ecatlas@assays$RNA@data[,iec_cells], 1, mean),
  std_iec = apply(ecatlas@assays$RNA@data[,iec_cells], 1, sd),
  avg_ecgc = apply(ecatlas@assays$RNA@data[,ecgc_cells], 1, mean),
  std_ecgc = apply(ecatlas@assays$RNA@data[,ecgc_cells], 1, sd)
),row.names = rownames(ecatlas@assays$RNA@data))
avg_exp$signal2noise <- (avg_exp$avg_iec - avg_exp$avg_ecgc) / 
  (avg_exp$std_iec + avg_exp$std_ecgc)

mk_5.8.10_4.6.14$signal2noise <- avg_exp[rownames(mk_5.8.10_4.6.14),"signal2noise"] 


write.csv(mk_5.8.10_4.6.14,'markers_iEC_vs_EC-GC.csv',row.names = F,quote = F)



saveRDS(ecatlas,'ecatlas_reclustered.RDS')
#ecatlas <- readRDS('ecatlas_reclustered.RDS')

podatlas <- readRDS('../podocyte_atlas.RDS')
podatlas <- FindVariableFeatures(podatlas)
podatlas <- RunPCA(FindVariableFeatures(podatlas))
podatlas <- RunHarmony(podatlas,group.by.vars = 'experiment')
podatlas <- RunUMAP(podatlas,reduction = 'harmony', dims = 1:30)

DimPlot(podatlas,group.by = 'subclass.l2')

podatlas <- FindNeighbors(object = podatlas, reduction = "harmony")
podatlas <- FindClusters(podatlas, resolution = c(0.2, 0.4, 0.6, 0.8, 1.0, 1.2))
DimPlot(podatlas,group.by = "RNA_snn_res.1.2",label = T)

write.csv(podatlas@meta.data,'podatlas_cluster_inglom_metadata_dec23.csv')

saveRDS(podatlas,'podatlas_reclustered.RDS')
#podatlas <- readRDS('podatlas_reclustered.RDS')

all_markers <- c('PTPRQ','WT1','NTNG1','NPHS1','NPHS2','CLIC5','PODXL',# POD
                 'CDKN1C','SPOCK2','B2M','CD81','S100A6','IGFBP2')# dPOD


Idents(podatlas) <- podatlas$RNA_snn_res.1.2
DimPlot(podatlas,label = T)
DotPlot(podatlas,all_markers)+
  theme(axis.text.x = element_text(angle = 90,vjust = .5))+
 scale_y_discrete(limits=factor(c(0:1,3:6,8:11,2,7)))

DimPlot(podatlas,label = T,group.by = 'subclass.l2')

markers <- list()
for (res in as.character(c(0.2, 0.4, 0.6, 0.8, 1.0, 1.2))){
  Idents(podatlas) <- paste0('RNA_snn_res.',res)
  markers[[res]] <- FindAllMarkers(podatlas)
}

saveRDS(markers,'markers_pod_resolutions.RDS')


# mk_5.8.10_4.6.14 <- FindMarkers(podatlas, ident.1 = c(5,8,10),
#                                 ident.2 = c(4,6,14))
# mk_5.8.10_4.6.14$genes <- rownames(mk_5.8.10_4.6.14)


podatlas@meta.data$new_types <- ifelse(podatlas@meta.data$RNA_snn_res.1.2 %in% c(7),
                                      'I-POD','POD')
DimPlot(podatlas,label = T,group.by = "new_types")
saveRDS(podatlas,'podatlas_reclustered.RDS')
#podatlas <- readRDS('podatlas_reclustered.RDS')


glom <- LoadH5Seurat('../../atlas_v1_paper/Kidney_Healthy-Injury_Cell_Atlas_snCv3_Seurat_noNA_04222023.h5Seurat')
glom@meta.data$new_types <- 'other'
glom@meta.data[rownames(podatlas@meta.data),"new_types"] <- podatlas@meta.data$new_types
glom@meta.data[rownames(mcatlas@meta.data),"new_types"] <- mcatlas@meta.data$new_types
glom@meta.data[rownames(ecatlas@meta.data),"new_types"] <- ecatlas@meta.data$new_types

Idents(glom) <- glom@meta.data$new_types
glom <- subset(glom,idents = 'other',invert=T)

glom <- FindVariableFeatures(glom)
glom <- RunPCA(FindVariableFeatures(glom))
glom <- RunHarmony(glom,group.by.vars = 'experiment')
glom <- RunUMAP(glom,reduction = 'harmony', dims = 1:30)


new_colors <- as.data.frame(list(`New types` = c('POD','I-POD','MC','I-VSMC','VSMC-NOS',
                                                 'EC-GC','I-EC','EC-NOS'),
                                   Colors = c("#DB7295", "#9370DB",
                                              "#483D8B", "#48E7EB", "#ADD8E6",
                                              "#FF8C00", "#F04681","#FFDAB9"
                                              )))

write.csv(new_colors,'colors_new_types.csv',row.names = F,quote = F)
new_colors <- read.csv('colors_new_types.csv')

glom@meta.data[glom@meta.data$new_types == 'I-POD','new_types'] <- 'POD'
new_colors <- new_colors[new_colors$New.types != 'I-POD',]

Idents(glom) <- factor(glom@meta.data$new_types,
                       levels = new_colors$New.types)
pdf('umap_glom_new_types.pdf',width = 6,height = 4.5)                       
DimPlot(glom,label = T,cols = new_colors$Colors)
dev.off()

genes <- c('PTPRQ','WT1','NTNG1','NPHS1','NPHS2','CLIC5','PODXL',
          'CDKN1C','SPOCK2','B2M','CD81','S100A6','IGFBP2',
          'NOTCH3','PDGFRB','ITGA8',
          'PIP5K1B','ROBO1','PIEZO2','DAAM2','PHTF2','GATA3','POSTN',
          'SLIT3','DACH1','PDE1C','IGFBP7','LAMC3','PLXDC1','PLXNA2',
          'CD34','PECAM1','PTPRB','MEIS2','FLT1','EMCN',
          'HECW2','PLAT','EHD3','KDR','SOST',
          'PKHD1','ESRRG','SPP1','HNF1B','PAX2')
pdf('dotplot_glom_new_types.pdf',width = 11,height = 3)
DotPlot(glom,genes)+
  theme(axis.text.x = element_text(angle=90,vjust = .5))+
  scale_y_discrete(limits=rev)
dev.off()

saveRDS(glom,'glom_reclustered.RDS')
glom <- readRDS('glom_reclustered.RDS')

write.csv(glom@meta.data,'glom_new_types_metadata_jan24.csv')

glom@meta.data$new_prolif <- glom@meta.data$new_types
glom@meta.data[rownames(ecatlas@meta.data),"new_prolif"] <- ecatlas@meta.data$new_prolif

DimPlot(glom,label = T,group.by = 'new_prolif')

write.csv(glom@meta.data,'glom_new_prolif_metadata_feb24.csv')


