library(Seurat) # Version 5.0.3
library(SeuratDisk)

f59 <- LoadXenium('../../xenium/Xenium_1_02122024/output-XETG00126__0010200__f59__20240214__210015/')
f59 <- subset(f59,subset = nCount_Xenium > 0)

VlnPlot(f59, features = c("nFeature_Xenium", "nCount_Xenium"), ncol = 2, pt.size = 0)

ImageDimPlot(f59, fov = "fov", molecules = c("NPHS2","PECAM1","TAGLN","SEMA6A","PLXNA2"), nmols = 20000)
ImageFeaturePlot(f59, features = c("NPHS2","PECAM1","TAGLN","SEMA6A","PLXNA2"))

f59 <- SCTransform(f59, assay = "Xenium")
f59 <- RunPCA(f59, npcs = 50, features = rownames(f59))
f59 <- RunUMAP(f59, dims = 1:50)
f59 <- FindNeighbors(f59, reduction = "pca", dims = 1:50)
f59 <- FindClusters(f59, resolution = 0.8)

DimPlot(f59)
FeaturePlot(f59, features = c("NPHS2","PECAM1","TAGLN","SEMA6A","PLXNA2"))

atlas <- LoadH5Seurat('../../atlas_v1_paper/Kidney_Healthy-Injury_Cell_Atlas_snCv3_Seurat_noNA_04222023.h5Seurat')

Idents(atlas) <- 'region.l2'
atlas2 <- subset(atlas,idents = c('P','M'),invert=T)

Idents(atlas2) <- 'state.l2'
atlas2 <- subset(atlas2,idents = c('degenerative'),invert=T)


mapped <- readRDS('../multiome_analysis/mesangial_snRNAApr_01_2024.RDS')

atlas2$predsnRNA0.5 <- atlas2$subclass.l2
atlas2@meta.data[rownames(mapped@meta.data),"predsnRNA0.5"] <- mapped$predsnRNA0.5
# unique(atlas2$predsnRNA0.5)
# 'pEC' %in% unique(atlas2$predsnRNA0.5)
rm(atlas)


####################
# anchors pipeline #
####################
atlas2 <- SCTransform(atlas2)
anchors <- FindTransferAnchors(reference = atlas2, query = f59, normalization.method = "SCT")
predictions.assay <- TransferData(anchorset = anchors, refdata = atlas2$predsnRNA0.5, prediction.assay = TRUE,
                                  weight.reduction = f59[["pca"]], dims = 1:50)
f59[["predictions"]] <- predictions.assay
predictions <- TransferData(anchorset = anchors, refdata = atlas2$predsnRNA0.5,
                            weight.reduction = f59[["pca"]], dims = 1:50)
f59@meta.data$anchors_predictions <- predictions[rownames(f59@meta.data),"predicted.id"]

sort(table(f59$anchors_predictions))
cell_groups <- data.frame(list(cell_id=rownames(f59@meta.data),group=f59$anchors_predictions))
write.csv(cell_groups,'f59_anchors.csv',row.names = F,quote = F)

f59@assays$predictions@data[,'lgeokdbp-1']

#################
# RCTD pipeline #
#################


library(spacexr)

query.counts <- GetAssayData(f59, assay = "Xenium", slot = "counts")[, Cells(f59)]
coords <- GetTissueCoordinates(f59, which = "centroids")
rownames(coords) <- coords$cell
coords$cell <- NULL
query <- SpatialRNA(coords, query.counts, colSums(query.counts))


sort(table(atlas2$predsnRNA0.5))
Idents(atlas2) <- atlas2$predsnRNA0.5
atlas2 <- subset(atlas2,idents = c('dDTL3','PapE'),invert=T)

counts <- GetAssayData(atlas2, assay = "RNA", slot = "counts")
cluster <- as.factor(atlas2$predsnRNA0.5)
names(cluster) <- colnames(atlas2)
nUMI <- atlas2$nCount_RNA
names(nUMI) <- colnames(atlas2)
nUMI <- colSums(counts)
levels(cluster) <- gsub("/", "-", levels(cluster))
reference <- Reference(counts, cluster, nUMI)

# run RCTD with many cores
RCTD <- create.RCTD(query, reference, max_cores = 8)
RCTD <- run.RCTD(RCTD, doublet_mode = "doublet")

annotations.df <- RCTD@results$results_df
annotations <- annotations.df$first_type
names(annotations) <- rownames(annotations.df)
f59$predicted.celltype <- annotations
keep.cells <- Cells(f59)[!is.na(f59$predicted.celltype)]
f59 <- subset(f59, cells = keep.cells)

head(f59@meta.data)
sort(table(f59$predicted.celltype))
cell_groups <- data.frame(list(cell_id=rownames(f59@meta.data),group=f59$predicted.celltype))
write.csv(cell_groups,'f59_RCTD.csv',row.names = F,quote = F)


####################

