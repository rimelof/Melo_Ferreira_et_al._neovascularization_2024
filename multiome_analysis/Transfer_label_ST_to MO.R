library(Seurat)
library(Signac)
#B4-Transfer label ST to Multiome
#MO Multiome
#st spatial
st <- readRDS("data/gloms_neighborhood_res2_020524.RDS")
ast <- UpdateSeuratObject(st)
st <- ast
ast <- NULL
DefaultAssay(st) <- 'SCT'
Idents(st) <- 'seurat_clusters' # new


obj0 <- readRDS('../../plot_ATAC/data/blue_mypathfrags_Oct2022.RDS')

DefaultAssay(MO) <- "RNA"
MO <- SCTransform(MO)
DefaultAssay(MO)
DefaultAssay(st)
# find anchors
anchors <- FindTransferAnchors(reference = st, query = MO)
# transfer labels
new_types <- TransferData(
  anchorset = anchors,
  refdata = st$seurat_clusters,prediction.assay = T)
MO[['pred.neighborhood']] <- new_types

new_types <- TransferData(
  anchorset = anchors,
  refdata = st$seurat_clusters)
MO <- AddMetaData(object = MO, metadata = new_types)

Idents(MO) <- MO$predicted.id
aux <- Idents(MO)
aux <- levels(aux[[1]])
aux <- as.numeric(aux)
aux <- max(aux)
seq(0,aux,1)
MO$predicted.id <- factor(MO$predicted.id,levels = seq(0,aux,1))
Idents(MO) <- MO$predicted.id
MO@meta.data$predneighborhood <- MO@meta.data$predicted.id
saveRDS(MO,paste0("MELOrds/B4_",date,".RDS"))
```