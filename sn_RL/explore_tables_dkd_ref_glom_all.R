library(CellChat)
library(ggplot2)


dkdcellchat <- readRDS('sn_cellchat/DKD_glom_all/cellchat_after_visualization.RDS')
refcellchat <- readRDS('sn_cellchat/ref_glom_all/cellchat_after_visualization.RDS')

ct <- 'EC-GC'

dkdnames <- dkdcellchat@LR$LRsig$interaction_name
dkd1p <- dkdcellchat@net$pval[ct,,]
dkd1prob <- dkdcellchat@net$pval[ct,,]

keep <- apply(dkd1p,2,function(x) {any(x < .01)}) & apply(dkd1prob,2,function(x) {any(x > 0)})
dkdnames <- dkdnames[keep]
dkd1p <- dkd1p[,keep]
dkd1prob <- dkd1prob[,keep]


refnames <- refcellchat@LR$LRsig$interaction_name
ref1p <- refcellchat@net$pval[ct,,]
ref1prob <- refcellchat@net$pval[ct,,]

keep <- apply(ref1p,2,function(x) {any(x < .01)}) & apply(ref1prob,2,function(x) {any(x > 0)})
refnames <- refnames[keep]
ref1p <- ref1p[,keep]
ref1prob <- ref1prob[,keep]

refnames[!refnames%in%dkdnames]

dkd1prob <- reshape2::melt(dkd1prob)
ref1prob <- reshape2::melt(ref1prob)

probdif <- merge(dkd1prob,ref1prob,by = c('Var1','Var2'),all = T)
colnames(probdif)[3:4] <- c('dkd','ref')
probdif[is.na(probdif)] <- 0

grep('angiogenesis',CellChatDB.human$interaction$ligand.keyword,ignore.case = T)
grep('angiogenesis',CellChatDB.human$interaction$receptor,ignore.case = T)
grep('vessel',CellChatDB.human$interaction$ligandreceptor.keyword)
grep('sema6a',rownames(CellChatDB.human$interaction),ignore.case = T)
CellChatDB.human$interaction[91,]
CellChatDB.human$interaction[c(1908,1910),]
unique(unlist(stringr::str_split(unique(CellChatDB.human$interaction$receptor.keyword),', ')))

grep('sema',rownames(CellChatDB.human$interaction[grep('angiogenesis',CellChatDB.human$interaction$ligand.keyword,ignore.case = T),]),ignore.case = T)



