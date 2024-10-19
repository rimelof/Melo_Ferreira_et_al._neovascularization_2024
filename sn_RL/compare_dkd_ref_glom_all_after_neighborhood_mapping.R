library(CellChat)
library(patchwork)

dkdcellchat <- readRDS('sn_cellchat_pEC/DKD_glom_all/cellchat_after_visualization.RDS')
refcellchat <- readRDS('sn_cellchat_pEC/ref_glom_all/cellchat_after_visualization.RDS')

data.dir <- './sn_cellchat_pEC/comparison_glom_all'
dir.create(data.dir)
setwd(data.dir)

object.list <- list(DKD = dkdcellchat, Ref = refcellchat)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))
cellchat

gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
gg1 + gg2


par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat, weight.scale = T)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")

gg1 <- netVisual_heatmap(cellchat)
#> Do heatmap based on a merged object
gg2 <- netVisual_heatmap(cellchat, measure = "weight")
#> Do heatmap based on a merged object
gg1 + gg2


weight.max <- getMaxWeight(object.list, attribute = c("idents","count"))
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
}


num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
gg <- list()
for (i in 1:length(object.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], title = names(object.list)[i], weight.MinMax = weight.MinMax)
}
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
patchwork::wrap_plots(plots = gg)


gg1 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "EC-GC")#, signaling.exclude = "MIF")
#> Visualizing differential outgoing and incoming signaling changes from NL to LS
#> The following `from` values were not present in `x`: 0
#> The following `from` values were not present in `x`: 0, -1
gg2 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "EC-NOS")#, signaling.exclude = c("MIF"))
#> Visualizing differential outgoing and incoming signaling changes from NL to LS
#> The following `from` values were not present in `x`: 0, 2
#> The following `from` values were not present in `x`: 0, -1
patchwork::wrap_plots(plots = list(gg1,gg2))

gg1 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "MC")#, signaling.exclude = "MIF")
#> Visualizing differential outgoing and incoming signaling changes from NL to LS
#> The following `from` values were not present in `x`: 0
#> The following `from` values were not present in `x`: 0, -1
gg2 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "VSMC-NOS")#, signaling.exclude = c("MIF"))
#> Visualizing differential outgoing and incoming signaling changes from NL to LS
#> The following `from` values were not present in `x`: 0, 2
#> The following `from` values were not present in `x`: 0, -1
patchwork::wrap_plots(plots = list(gg1,gg2))

gg1 <- rankNet(cellchat, mode = "comparison", measure = "weight", sources.use = NULL, targets.use = NULL, stacked = T, do.stat = TRUE)
gg2 <- rankNet(cellchat, mode = "comparison", measure = "weight", sources.use = NULL, targets.use = NULL, stacked = F, do.stat = TRUE)

gg1 + gg2


### Part III: Identify the up-gulated and down-regulated signaling ligand-receptor pairs
netVisual_bubble(cellchat, sources.use = c('EC-GC','EC-NOS','MC','MC-NOS'), 
                 targets.use = 1:7,  comparison = c(1, 2), angle.x = 45)
netVisual_bubble(cellchat, sources.use = 1:7, 
                 targets.use = c('EC-GC','EC-NOS','MC','MC-NOS'),  comparison = c(1, 2), angle.x = 45)

gg1 <- netVisual_bubble(cellchat, sources.use = c('EC-GC','EC-NOS','MC','MC-NOS'), comparison = c(1, 2), 
                        max.dataset = 1, title.name = "Increased signaling in DKD", angle.x = 45, remove.isolate = T)
#> Comparing communications on a merged object
gg2 <- netVisual_bubble(cellchat, sources.use = c('EC-GC','EC-NOS','MC','MC-NOS'), comparison = c(1, 2), 
                        max.dataset = 2, title.name = "Decreased signaling in DKD", angle.x = 45, remove.isolate = T)
#> Comparing communications on a merged object
gg1 + gg2

gg1 <- netVisual_bubble(cellchat, targets.use = c('EC-GC','EC-NOS','MC','MC-NOS'), comparison = c(1, 2), 
                        max.dataset = 1, title.name = "Increased signaling in DKD", angle.x = 45, remove.isolate = T)
#> Comparing communications on a merged object
gg2 <- netVisual_bubble(cellchat, targets.use = c('EC-GC','EC-NOS','MC','MC-NOS'), comparison = c(1, 2), 
                        max.dataset = 2, title.name = "Decreased signaling in DKD", angle.x = 45, remove.isolate = T)
#> Comparing communications on a merged object
gg1 + gg2



### Identify dysfunctional signaling by using differential expression analysis
# define a positive dataset, i.e., the dataset with positive fold change against the other dataset
pos.dataset = "DKD"
# define a char name used for storing the results of differential expression analysis
features.name = paste0(pos.dataset, ".merged")

# perform differential expression analysis 
# Of note, compared to CellChat version < v2, CellChat v2 now performs an ultra-fast Wilcoxon test using the presto package, which gives smaller values of logFC. Thus we here set a smaller value of thresh.fc compared to the original one (thresh.fc = 0.1). Users can also provide a vector and dataframe of customized DEGs by modifying the cellchat@var.features$LS.merged and cellchat@var.features$LS.merged.info. 

cellchat <- identifyOverExpressedGenes(cellchat, group.dataset = "datasets", pos.dataset = pos.dataset, features.name = features.name, only.pos = FALSE, thresh.pc = 0.1, thresh.fc = 0.05,thresh.p = 0.05, group.DE.combined = FALSE) 
#> Use the joint cell labels from the merged CellChat object

# map the results of differential expression analysis onto the inferred cell-cell communications to easily manage/subset the ligand-receptor pairs of interest
net <- netMappingDEG(cellchat, features.name = features.name, variable.all = TRUE)
# extract the ligand-receptor pairs with upregulated ligands in LS
net.up <- subsetCommunication(cellchat, net = net, datasets = "DKD",ligand.logFC = 0.05, receptor.logFC = NULL)
# extract the ligand-receptor pairs with upregulated ligands and upregulated receptors in NL, i.e.,downregulated in LS
net.down <- subsetCommunication(cellchat, net = net, datasets = "Ref",ligand.logFC = -0.05, receptor.logFC = NULL)

write.csv(net,'communication_all.csv',row.names = F,quote = F)
write.csv(net.up,'communication_up_in_dkd.csv',row.names = F,quote = F)
write.csv(net.down,'communication_down_in_dkd.csv',row.names = F,quote = F)

gene.up <- extractGeneSubsetFromPair(net.up, cellchat)
gene.down <- extractGeneSubsetFromPair(net.down, cellchat)

pairLR.use.up = net.up[, "interaction_name", drop = F]
gg1 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.up, sources.use = c('EC-GC','EC-NOS','MC','MC-NOS'), comparison = c(1, 2),  angle.x = 90, remove.isolate = T,title.name = "Up-regulated signaling in DKD")
#> Comparing communications on a merged object
pairLR.use.down = net.down[, "interaction_name", drop = F]
gg2 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.down, sources.use = c('EC-GC','EC-NOS','MC','MC-NOS'), comparison = c(1, 2),  angle.x = 90, remove.isolate = T,title.name = "Down-regulated signaling in DKD")
#> Comparing communications on a merged object
gg1 + gg2


lr_up <- paste0(net.up$interaction_name)
terms <- c(unlist(stringr::str_split(unique(CellChatDB.human$interaction[lr_up,'ligand.keyword']),', ')),
           unlist(stringr::str_split(unique(CellChatDB.human$interaction[lr_up,'receptor.keyword']),', ')))
#sort(table(terms))

lr_down <- paste0(net.down$interaction_name)
terms <- c(terms,
           unlist(stringr::str_split(unique(CellChatDB.human$interaction[lr_down,'ligand.keyword']),', ')),
           unlist(stringr::str_split(unique(CellChatDB.human$interaction[lr_down,'receptor.keyword']),', ')))
sort(table(terms))

target_terms <- c('Developmental protein','Cell junction','Cell projection','Differentiation','Angiogenesis',
                  'EGF-like domain','Myogenesis','Laminin EGF-like domain','Apoptosis','Notch signaling pathway',
                  'Growth factor')

net.up$ligand.keyword <- CellChatDB.human$interaction[lr_up,'ligand.keyword']
net.up$receptor.keyword <- CellChatDB.human$interaction[lr_up,'receptor.keyword']
net.down$ligand.keyword <- CellChatDB.human$interaction[lr_down,'ligand.keyword']
net.down$receptor.keyword <- CellChatDB.human$interaction[lr_down,'receptor.keyword']

ind <- unique(c(unlist(lapply(target_terms,function(x) {grep(x,net.up$ligand.keyword)})),
                unlist(lapply(target_terms,function(x) {grep(x,net.up$receptor.keyword)}))))
interactions <- net.up[ind,'interaction_name']
intup <- as.character(net.up[ind,])
ind <- unique(c(unlist(lapply(target_terms,function(x) {grep(x,net.down$ligand.keyword)})),
                unlist(lapply(target_terms,function(x) {grep(x,net.down$receptor.keyword)}))))
interactions <- rbind(interactions,net.down[ind,])
intdown <- as.character(net.down[ind,'interaction_name'])

pdf('buble_target_interactions_all_dkd_ref.pdf',height = 16,width = 10)
netVisual_bubble(cellchat, pairLR.use = interactions[, "interaction_name", drop = F], comparison = c(1, 2),
                 angle.x = 90, remove.isolate = T,title.name = "Target interactions")
dev.off()

pdf('buble_target_interactions_sourcesEC_dkd_ref.pdf',height = 16,width = 10)
netVisual_bubble(cellchat, pairLR.use = interactions[, "interaction_name", drop = F], 
                 comparison = c(1, 2), sources.use = c('EC-GC','pEC','EC-NOS','I-EC'),
                 angle.x = 90, remove.isolate = T,title.name = "Target interactions")
dev.off()

pdf('buble_target_interactions_targetsEC_dkd_ref.pdf',height = 16,width = 10)
netVisual_bubble(cellchat, pairLR.use = interactions[, "interaction_name", drop = F], 
                 comparison = c(1, 2),targets.use = c('EC-GC','pEC','EC-NOS','I-EC'),
                 angle.x = 90, remove.isolate = T,title.name = "Target interactions")
dev.off()

ec_deg <- read.csv('../../../neighborhoods/neigborhood_8vs2018_dkd_spatial_markers.csv')
ec_deg[ec_deg$genes %in% c(net.up$ligand,net.up$receptor,net.down$ligand,net.down$receptor),"genes"]
ec_deg[ec_deg$genes %in% c(net$ligand,net$receptor),"genes"]

pod_deg <- read.csv('../../../neighborhoods/neigborhood_1511vs31419_dkd_spatial_markers.csv')
pod_deg[pod_deg$genes %in% c(net.up$ligand,net.up$receptor,net.down$ligand,net.down$receptor),"genes"]
pod_deg[pod_deg$genes %in% c(net$ligand,net$receptor),"genes"]

st_rl <- read.csv('../../../spatial_RL_commot_155/boruta_signif_rlpairs_02202024.csv')
st_interactions <- net[net$ligand %in% st_rl$V1 | net$receptor %in% st_rl$V2,]
pdf('buble_ST_interactions_sourcesEC_dkd_ref.pdf',height = 10,width = 10)
netVisual_bubble(cellchat, pairLR.use = st_interactions[, "interaction_name", drop = F], 
                 comparison = c(1, 2), sources.use = c('EC-GC','pEC','EC-NOS','I-EC'),
                 angle.x = 90, remove.isolate = T,title.name = "Target interactions")
dev.off()
pdf('buble_ST_interactions_tagetsEC_dkd_ref.pdf',height = 10,width = 10)
netVisual_bubble(cellchat, pairLR.use = st_interactions[, "interaction_name", drop = F], 
                 comparison = c(1, 2), targets.use = c('EC-GC','pEC','EC-NOS','I-EC'),
                 angle.x = 90, remove.isolate = T,title.name = "Target interactions")
dev.off()
pdf('buble_ST_interactions_sourcesVSMC_dkd_ref.pdf',height = 10,width = 10)
netVisual_bubble(cellchat, pairLR.use = st_interactions[, "interaction_name", drop = F], 
                 comparison = c(1, 2), sources.use = c('MC','VSMC-NOS','I-VSMC'),
                 angle.x = 90, remove.isolate = T,title.name = "Target interactions")
dev.off()
pdf('buble_ST_interactions_tagetsVSMC_dkd_ref.pdf',height = 10,width = 10)
netVisual_bubble(cellchat, pairLR.use = st_interactions[, "interaction_name", drop = F], 
                 comparison = c(1, 2), targets.use = c('MC','VSMC-NOS','I-VSMC'),
                 angle.x = 90, remove.isolate = T,title.name = "Target interactions")
dev.off()

pairs <- c('VEGFA','SEMA6A', 'LRP6','FN1')#,'PDGF','PLXNA2','GAS6',
inters <- CellChatDB.human$interaction[, "interaction_name", drop = F]
pairs <- inters[grep(paste(pairs,collapse = '|'),inters$interaction_name),,drop = F]
rownames(pairs) <- 1:nrow(pairs)
#pairs <- pairs[c(1,17,21,29,26,24,32,34,37),,drop=F]
pdf('buble_target_selected_interactions_all_dkd_ref.pdf',height = 6,width = 10)
netVisual_bubble(cellchat, pairLR.use = pairs, comparison = c(1, 2),
                 sources.use = c('POD','EC-GC','pEC','I-EC','EC-NOS'),#,'MC', 'I-VSMC','VSMC-NOS'),
                 targets.use = c('POD','EC-GC','pEC','I-EC','EC-NOS','MC', 'I-VSMC','VSMC-NOS'),
                 sort.by.source = T,sort.by.target = T,
                 angle.x = 90, remove.isolate = F,title.name = "Target interactions")
dev.off()
viz <- netVisual_bubble(cellchat, pairLR.use = pairs, comparison = c(1, 2),
                        angle.x = 90, remove.isolate = F,title.name = "Target interactions",return.data = T)

write.table(unique(c(stringr::str_split(pairs$interaction_name,'_',simplify = T))),'selected_rl.csv',
            quote = F,row.names = F,col.names = F)


pairs <- CellChatDB.human$interaction[, "interaction_name", drop = F]
rownames(pairs) <- 1:nrow(pairs)
#pairs <- pairs[c(1,17,21,29,26,24,32,34,37),,drop=F]

pdf('test_all_pEC_source.pdf',height = 20)
netVisual_bubble(cellchat, pairLR.use = pairs, comparison = c(1, 2),
                 sources.use = c('pEC'),
                 targets.use = c('POD','EC-GC','pEC','I-EC','EC-NOS','MC', 'I-VSMC','VSMC-NOS'),
                 sort.by.source = T,sort.by.target = T,
                 angle.x = 90, remove.isolate = F,title.name = "Target interactions")
dev.off()

pdf('test_all_pEC_target.pdf',height = 8)
netVisual_bubble(cellchat, pairLR.use = pairs, comparison = c(1, 2),
                 sources.use = c('POD','EC-GC','pEC','I-EC','EC-NOS','MC', 'I-VSMC','VSMC-NOS'),
                 targets.use = c('pEC'),
                 sort.by.source = T,sort.by.target = T,
                 angle.x = 90, remove.isolate = F,title.name = "Target interactions")
dev.off()


pairs <- c('VEGFA','PLXNA2')#,'PDGF','PLXNA2','GAS6',
inters <- CellChatDB.human$interaction[, "interaction_name", drop = F]
pairs <- inters[grep(paste(pairs,collapse = '|'),inters$interaction_name),,drop = F]
rownames(pairs) <- 1:nrow(pairs)

pdf('bubble_iEC_story.pdf',height = 4,width = 10)
netVisual_bubble(cellchat, pairLR.use = pairs, comparison = c(1, 2),
                 sources.use = c('POD','EC-GC','I-EC','EC-NOS'),
                 targets.use = c('EC-GC','pEC','I-EC','EC-NOS','MC', 'VSMC-NOS'),
                 sort.by.source = T,sort.by.target = T,
                 angle.x = 90, remove.isolate = F,title.name = "Target interactions")
dev.off()

pairs <- c('NOTCH')#,'PDGF','PLXNA2','GAS6',
inters <- CellChatDB.human$interaction[, "interaction_name", drop = F]
pairs <- inters[grep(paste(pairs,collapse = '|'),inters$interaction_name),,drop = F]
rownames(pairs) <- 1:nrow(pairs)

netVisual_bubble(cellchat, pairLR.use = pairs, comparison = c(1, 2),
                 sources.use = c('POD','EC-GC','I-EC','EC-NOS','pEC'),
                 targets.use = c('POD','EC-GC','pEC','I-EC','EC-NOS','MC', 'VSMC-NOS','I-VSMC'),
                 sort.by.source = T,sort.by.target = T,
                 angle.x = 90, remove.isolate = F,title.name = "Target interactions")

pdf('bubble_pec_notch_story.pdf',height = 4,width = 6)
netVisual_bubble(cellchat, pairLR.use = pairs, comparison = c(1, 2),
                 sources.use = c('pEC'),
                 targets.use = c('POD','EC-GC','pEC','I-EC','EC-NOS','MC', 'VSMC-NOS','I-VSMC'),
                 sort.by.source = T,sort.by.target = T,
                 angle.x = 90, remove.isolate = F,title.name = "Target interactions")
dev.off()


pairs <- c('WNT2B','FN1')
inters <- CellChatDB.human$interaction[, "interaction_name", drop = F]
pairs <- inters[grep(paste(pairs,collapse = '|'),inters$interaction_name),,drop = F]
rownames(pairs) <- 1:nrow(pairs)

pdf('bubble_pec_FN1_WNT2B_story.pdf',height = 4,width = 6)
netVisual_bubble(cellchat, pairLR.use = pairs, comparison = c(1, 2),
                 sources.use = c('pEC'),
                 targets.use = c('POD','EC-GC','pEC','I-EC','EC-NOS','MC', 'VSMC-NOS','I-VSMC'),
                 sort.by.source = T,sort.by.target = T,
                 angle.x = 90, remove.isolate = F,title.name = "Target interactions")
dev.off()
