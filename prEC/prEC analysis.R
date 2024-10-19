library(Seurat)
library(rio)
library(ggplot2)

atlas <- readRDS('../multiome_analysis/mesangial_snRNAApr_01_2024.RDS')

degs <- import('DEGS_atlas_Multiome_pECxEC-GC_EC-NS_I-EC.xlsx')
degs[is.na(degs$...8),"...8"] <- ' '

genes <- degs[degs$...8 == 'x',"gene"]
genes <- c(genes,'PLAT','SOST','ESRRG','PKHD1','SPP1')
genes <- genes[c(1,8,12,2:7,9:11,13,14,16,21,23,24,
                 15,17,25,26,28,22,18,27,20,29:33,19)]

ecs <- subset(atlas,idents = c('pEC','EC-GC','I-EC','EC-NOS'))
Idents(ecs) <- factor(ecs$predsnRNA0.5,levels = rev(c('pEC','EC-GC','I-EC','EC-NOS')))

pdf('dotplot_prEC.pdf',width = 8,height = 2)
DotPlot(ecs,genes,idents = c('pEC','EC-GC','I-EC','EC-NOS'))+
  theme(axis.text.x = element_text(angle = 90,vjust = .5))
dev.off()

kos <- import('MEF2A_MEF2C_TRPS1_DEGs_lognorm_KO_and_NO_KO_pseudbulk.csv')
rownames(kos) <- kos$KO_names
kos <- kos[genes,]

kos$KO_names <- factor(kos$KO_names,levels = genes)
pdf('dotplot_KO_FC.pdf',width = 8,height = 2)
ggplot(kos,aes(x=KO_names,y=KO_logfoldchanges,size = -log10(KO_pvals_adj)))+
  geom_point()+
  geom_hline(yintercept=0, lty=2)+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90,vjust = .5))
dev.off()
