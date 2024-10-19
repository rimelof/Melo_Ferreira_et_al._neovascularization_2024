library(ggplot2)
library(Seurat)
library(EnhancedVolcano)


coltable <- read.csv('../identify_glom_cell_types/colors_new_types.csv')
coltable <- rbind(coltable,
                  list(New.types = 'Other',Colors = '#DDDDDD'))
coltable[3,1] <- 'VSM/P'
coltable[6,1] <- 'EC'
coltable <- coltable[c(1,3,6,9),]
rownames(coltable) <- coltable$New.types

gloms <- readRDS('glom_atlas_manuscript_012724.RDS')

Idents(gloms) <- gloms$condition.l2

degglom <- FindMarkers(gloms,ident.1 = 'DKD',ident.2 = 'Ref')
degglom$genes <- rownames(degglom)

pdf('volcano_DKD_vs_Ref.pdf')
plot(  EnhancedVolcano(degglom,x = "avg_log2FC",y = "p_val_adj",lab = degglom$genes,
                       title = NULL,subtitle = NULL,legendPosition = 'none',
                       pCutoff = 0.05,FCcutoff = .5))
dev.off()



Idents(gloms) <- gloms$subclass.l1
for (ct in c('POD','EC','VSM/P')){
  obj <- subset(gloms,idents = ct)
  Idents(obj) <- obj$condition.l2
  
  degglom <- FindMarkers(obj,ident.1 = 'DKD',ident.2 = 'Ref')
  degglom$genes <- rownames(degglom)
  
  pdf(paste0('volcano_',stringr::str_replace(ct,'/','-'),'_DKD_vs_Ref.pdf'))
  plot(  EnhancedVolcano(degglom,x = "avg_log2FC",y = "p_val_adj",lab = degglom$genes,
                         title = NULL,subtitle = NULL,legendPosition = 'none',
                         pCutoff = 0.05,FCcutoff = .5,
                         col=c("grey30","grey30","grey30",coltable[ct,"Colors"])))
  dev.off()
}
