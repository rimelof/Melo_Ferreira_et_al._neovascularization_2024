library(ggplot2)
library(Seurat)
library(EnhancedVolcano)
library(clusterProfiler)


gloms <- readRDS('../neighborhoods/gloms_ref_DKD_w_glom_index.RDS')

gtab <- as.data.frame(table(gloms@meta.data[,c("glom_index","orig.ident")]))
gtab <- gtab[gtab$Freq > 0,]

table(gtab$orig.ident)

Idents(gloms) <- gloms$orig.ident

img <- stringr::str_replace_all('V13J17-332_XY03_23-0123',pattern = '-','.')
pdf('glom.pdf')
SpatialDimPlot(gloms,images = img,crop = F)
dev.off()
