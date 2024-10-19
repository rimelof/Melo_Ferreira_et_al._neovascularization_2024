library(Seurat) # Version 5.0.3
library(SeuratDisk)



atlas <- LoadH5Seurat('../../atlas_v1_paper/Kidney_Healthy-Injury_Cell_Atlas_snCv3_Seurat_noNA_04222023.h5Seurat')

Idents(atlas) <- 'region.l2'
atlas2 <- subset(atlas,idents = c('P','M'),invert=T)


mapped <- readRDS('../multiome_analysis/mesangial_snRNAApr_01_2024.RDS')

atlas2$predsnRNA0.5 <- atlas2$subclass.l2
atlas2@meta.data[rownames(mapped@meta.data),"predsnRNA0.5"] <- mapped$predsnRNA0.5

Idents(atlas2) <- atlas2$predsnRNA0.5
atlas2 <- SCTransform(atlas2)

atlas2 <- subset(atlas2,
                 idents = c('dPT','dPT/DTL','dM-TAL','dM-PC','dFIB','dC-TAL',
                            'dATL','dC-IC-A','dCNT','dIMCD','dM-FIB','dDTL3','dDCT'),
                 invert=T)

# unique(atlas2$predsnRNA0.5)
# 'pEC' %in% unique(atlas2$predsnRNA0.5)
# rm(atlas)


dirs <- list.dirs('../../xenium/Xenium_1_02122024',full.names = F,recursive = F)
samples <- stringr::str_split(dirs,'__',simplify = T)[,3] #double _
names(dirs) <- samples
list_xenium <- list()
for (smp in samples){
  
  xenium <- LoadXenium(paste0('../../xenium/Xenium_1_02122024/',dirs[smp]))
  xenium <- subset(xenium,subset = nCount_Xenium > 0)
  
  # VlnPlot(xenium, features = c("nFeature_Xenium", "nCount_Xenium"), ncol = 2, pt.size = 0)
  # 
  # ImageDimPlot(xenium, fov = "fov", molecules = c("NPHS2","PECAM1","TAGLN","SEMA6A","PLXNA2"), nmols = 20000)
  # ImageFeaturePlot(xenium, features = c("NPHS2","PECAM1","TAGLN","SEMA6A","PLXNA2"))
  
  xenium <- SCTransform(xenium, assay = "Xenium")
  xenium <- RunPCA(xenium, npcs = 50, features = rownames(xenium))
  xenium <- RunUMAP(xenium, dims = 1:50)
  xenium <- FindNeighbors(xenium, reduction = "pca", dims = 1:50)
  xenium <- FindClusters(xenium, resolution = 0.8)
  
  # DimPlot(xenium)
  # FeaturePlot(xenium, features = c("NPHS2","PECAM1","TAGLN","SEMA6A","PLXNA2"))
  
  
  ####################
  # anchors pipeline #
  ####################
  
  anchors <- FindTransferAnchors(reference = atlas2, query = xenium, normalization.method = "SCT")
  predictions.assay <- TransferData(anchorset = anchors, refdata = atlas2$predsnRNA0.5, prediction.assay = TRUE,
                                    weight.reduction = xenium[["pca"]], dims = 1:50)
  xenium[["predictions"]] <- predictions.assay
  predictions <- TransferData(anchorset = anchors, refdata = atlas2$predsnRNA0.5,
                              weight.reduction = xenium[["pca"]], dims = 1:50)
  xenium@meta.data$anchors_predictions <- predictions[rownames(xenium@meta.data),"predicted.id"]
  
  sort(table(xenium$anchors_predictions))
  cell_groups <- data.frame(list(cell_id=rownames(xenium@meta.data),group=xenium$anchors_predictions))
  write.csv(cell_groups,paste0(smp,'_anchors.csv'),row.names = F,quote = F)
  
  list_xenium[[smp]] <- xenium
}

saveRDS(list_xenium,'list_xenium.RDS')
