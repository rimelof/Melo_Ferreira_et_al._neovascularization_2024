library(ggplot2)
library(Seurat)
library(EnhancedVolcano)


coltable <- read.csv('../identify_glom_cell_types/colors_new_types.csv')
coltable <- rbind(coltable,
                  list(New.types = 'Other',Colors = '#DDDDDD'))
coltable[3,1] <- 'VSM-P'
coltable[6,1] <- 'EC'
coltable <- coltable[c(1,3,6,9),]
rownames(coltable) <- coltable$New.types

gloms <- readRDS('../figure1_analysis/glom_atlas_manuscript_012724.RDS')

comparisons <- matrix(c('I-POD','I-EC','I-VSMC','POD','EC-GC','MC','POD','EC','VSM-P'),ncol=3)

glommeta <- read.csv('../identify_glom_cell_types/glom_new_types_metadata_jan24.csv')
glommeta$seurat_clusters
rownames(glommeta) <- glommeta$X
gloms@meta.data$new_types <- gloms@meta.data$subclass.l2
gloms@meta.data[rownames(glommeta),'new_types'] <- as.character(glommeta$new_types)



Idents(gloms) <- gloms$new_types
deglist <- list()
for (ct in 1:3){
  degglom <- FindMarkers(gloms,ident.1 = comparisons[ct,1],ident.2 = comparisons[ct,2])
  degglom$genes <- rownames(degglom)

  deglist[[comparisons[ct,3]]] <- degglom
  write.csv(degglom,paste0('degglom_',comparisons[ct,3],'_.csv'),quote = F,row.names = F)
    
  pdf(paste0('volcano_',comparisons[ct,1],'_vs_',comparisons[ct,2],'.pdf'))
  plot(  EnhancedVolcano(degglom,x = "avg_log2FC",y = "p_val_adj",lab = degglom$genes,
                         title = NULL,subtitle = NULL,legendPosition = 'none',
                         pCutoff = 0.05,FCcutoff = .5,
                         col=c("grey30","grey30","grey30",coltable[comparisons[ct,3],"Colors"])))
  dev.off()
}
saveRDS(deglist,'deglist.rds')
