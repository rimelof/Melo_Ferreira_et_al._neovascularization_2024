library(Seurat)
library(SeuratDisk)

atlas <- LoadH5Seurat('../../atlas_v1_paper/Kidney_Healthy-Injury_Cell_Atlas_snCv3_Seurat_noNA_04222023.h5Seurat')
atlas <- UpdateSeuratObject(atlas)
atlas[['RNA']] <- as(object = atlas[['RNA']], Class = 'SCTAssay')
DefaultAssay(atlas) <- 'RNA'

atlas@meta.data$glomtypes <- ifelse(atlas@meta.data$subclass.l1 %in% c('POD','VSM/P','EC'),
                                    atlas@meta.data$subclass.l1,'Other')
Idents(atlas) <- factor(atlas@meta.data$glomtypes,levels=c('POD','VSM/P','EC','Other'))
cols <- c("#DB7295","#483D8B","#FF8C00","#BBBBBB")
pdf('umap_glom_level1.pdf',width = 6,height = 4.5)                       
DimPlot(atlas,label = T,cols = cols)
dev.off()


glom <- subset(atlas,idents=c('POD','VSM/P','EC'))
saveRDS(glom,'glom_atlas_manuscript_012724.RDS')
