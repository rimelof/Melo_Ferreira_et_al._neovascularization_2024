library(CellChat)
library(NMF)
library(ggalluvial)



cond <- 'DKD'

load(paste0('sn_cellchat_pEC/',cond,'_glom_all/cellchat.RData'))
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
pdf(paste0('netanalys_',cond,'_glom_all.pdf'))
netAnalysis_signalingRole_network(cellchat, signaling = cellchat@netP$pathways, width = 8, height = 2.5, font.size = 10)
dev.off()

pdf(paste0('netanalys_path_',cond,'_glom_all.pdf'))
ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing")
ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming")
ht1 + ht2
dev.off()

ht <- netAnalysis_signalingRole_heatmap(cellchat, signaling = c("SEMA6"))
ht


# library(NMF)
# library(ggalluvial)
selectK(cellchat, pattern = "outgoing")

nPatterns = 7

pdf(paste0('patterns_out',cond,'_glom_all.pdf'))
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "outgoing", k = nPatterns)
dev.off()


pdf(paste0('patterns_out_river',cond,'_glom_all.pdf'))
netAnalysis_river(cellchat, pattern = "outgoing")
dev.off()

pdf(paste0('patterns_out_dots',cond,'_glom_all.pdf'))
netAnalysis_dot(cellchat, pattern = "outgoing")
dev.off()

pdf(paste0('patterns_in',cond,'_glom_all.pdf'))
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "incoming", k = nPatterns)
dev.off()


pdf(paste0('patterns_in_river',cond,'_glom_all.pdf'))
netAnalysis_river(cellchat, pattern = "incoming")
dev.off()

pdf(paste0('patterns_in_dots',cond,'_glom_all.pdf'))
netAnalysis_dot(cellchat, pattern = "incoming")
dev.off()

#Manifold and classification learning analysis of signaling networks
cellchat <- computeNetSimilarity(cellchat, type = "functional")
cellchat <- netEmbedding(cellchat, type = "functional")
#> Manifold learning of the signaling networks for a single dataset
cellchat <- netClustering(cellchat, type = "functional")
#> Classification learning of the signaling networks for a single dataset
# Visualization in 2D-space
pdf(paste0('manifold_functional_',cond,'_glom_all.pdf'))
netVisual_embedding(cellchat, type = "functional", label.size = 3.5)
dev.off()

saveRDS(cellchat,paste0('sn_cellchat_pEC/',cond,'_glom_all/cellchat_after_visualization.RDS'))


cond <- 'Ref'

load(paste0('sn_cellchat_pEC/',cond,'_glom_all/cellchat.RData'))
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
pdf(paste0('netanalys_',cond,'_glom_all.pdf'))
netAnalysis_signalingRole_network(cellchat, signaling = cellchat@netP$pathways, width = 8, height = 2.5, font.size = 10)
dev.off()

pdf(paste0('netanalys_path_',cond,'_glom_all.pdf'))
ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing")
ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming")
ht1 + ht2
dev.off()

ht <- netAnalysis_signalingRole_heatmap(cellchat, signaling = c("SEMA6"))
ht


# library(NMF)
# library(ggalluvial)
selectK(cellchat, pattern = "outgoing")

nPatterns = 7

pdf(paste0('patterns_out',cond,'_glom_all.pdf'))
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "outgoing", k = nPatterns)
dev.off()


pdf(paste0('patterns_out_river',cond,'_glom_all.pdf'))
netAnalysis_river(cellchat, pattern = "outgoing")
dev.off()

pdf(paste0('patterns_out_dots',cond,'_glom_all.pdf'))
netAnalysis_dot(cellchat, pattern = "outgoing")
dev.off()

pdf(paste0('patterns_in',cond,'_glom_all.pdf'))
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "incoming", k = nPatterns)
dev.off()


pdf(paste0('patterns_in_river',cond,'_glom_all.pdf'))
netAnalysis_river(cellchat, pattern = "incoming")
dev.off()

pdf(paste0('patterns_in_dots',cond,'_glom_all.pdf'))
netAnalysis_dot(cellchat, pattern = "incoming")
dev.off()

#Manifold and classification learning analysis of signaling networks
cellchat <- computeNetSimilarity(cellchat, type = "functional")
cellchat <- netEmbedding(cellchat, type = "functional")
#> Manifold learning of the signaling networks for a single dataset
cellchat <- netClustering(cellchat, type = "functional")
#> Classification learning of the signaling networks for a single dataset
# Visualization in 2D-space
pdf(paste0('manifold_functional_',cond,'_glom_all.pdf'))
netVisual_embedding(cellchat, type = "functional", label.size = 3.5)
dev.off()

saveRDS(cellchat,paste0('sn_cellchat_pEC/',cond,'_glom_all/cellchat_after_visualization.RDS'))

