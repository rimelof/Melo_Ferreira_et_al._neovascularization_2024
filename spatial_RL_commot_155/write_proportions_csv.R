library(Seurat)

dir.create('proportions')

samples <- read.csv('../Samples_Table.csv')


for (r in 1:nrow(samples)){
  smp <- samples[r,"Sample"]
  spatial <- readRDS(paste0('../../spatial/share_',smp,'/',smp,'.RDS'))
#  spatial <- subset(stmerge,idents = smp)
  prop <- spatial@assays$predsubclassl2@data[1:74,]
  prop <- t(t(prop)/colSums(prop))
  rownames(prop) <- stringr::str_replace_all(rownames(prop),'/','_')
  write.csv(t(prop),paste0('proportions/',smp,'_proportion.csv'))
}

write.csv(samples$Sample,'sample_list.csv',quote = F,row.names = F)
