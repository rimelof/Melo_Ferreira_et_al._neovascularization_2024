library(Seurat)
library(ggplot2)
library(rio)
library(readxl)

list_sample <- list.files('../../gloms/','*csv')
list_sample <- stringr::str_remove(list_sample,'.csv')

samples <- import('../../create_glom_ST/Processing_checklist_all.xlsx',range=cell_limits(c(2, 1), c(NA, NA)))
samples <- samples[samples$Sample %in% list_sample,]

spatial_list <- list()
for (r in 1:nrow(samples)){
  smp <- samples[r,"Sample"]
  spatial <- readRDS(paste0(smp,'_ec.RDS'))
  spatial@meta.data$orig.ident <- smp
  spatial@meta.data$condition <- samples[r,"Condition"]
  gloms <- read.csv(paste0('../../gloms/',smp,'.csv'))
  rownames(gloms) <- gloms$Barcode
  spatial@meta.data$Gloms <- gloms[rownames(spatial@meta.data),'Gloms']
  print(smp)
  print(any(is.na(spatial$Gloms)))
  print(sort(unique(spatial$Gloms)))
  spatial_list[[samples[r,"Sample"]]] <- spatial
}

stmerged <- merge(spatial_list[[1]],spatial_list[2:33],
                  add.cell.ids=stringr::str_split(samples$Sample,'_',simplify = T)[,3])
print(c(length(colnames(stmerged)),length(Idents(stmerged))))

names(stmerged@images) <- stringr::str_replace_all(samples,'-','.')

#stmerged <- SCTransform(stmerged,assay = 'Spatial')
saveRDS(stmerged,'merged_33_after_ec.RDS')
#stmerged <- readRDS('merged_33_after_ec.RDS')



DefaultAssay(stmerged) <- 'pred.ec1.2'

# pdf('vln_glom_ti_res1.2.pdf')
# VlnPlot(stmerged,features = as.character(0:10),group.by = 'gloms')
# dev.off()

sort(rownames(stmerged))

map_prop_data <- as.data.frame(matrix(0,nrow=15,ncol=6,
                                      dimnames = list(as.character(0:14),
                                                      c('In Glom','In TI',
                                                        'G Rel spot','T Rel spot',
                                                        'G Rel cl','T Rel cl'))))

glspots <- rownames(stmerged@meta.data[stmerged$Gloms == 'glom',])
tispots <- rownames(stmerged@meta.data[stmerged$Gloms != 'glom',])
for (cl in as.character(0:14)){
  map_prop_data[cl,'In Glom'] <- sum(stmerged@assays$pred.ec1.2@data[cl,glspots])
  map_prop_data[cl,'In TI'] <- sum(stmerged@assays$pred.ec1.2@data[cl,tispots])
  map_prop_data[cl,"G Rel spot"] <- map_prop_data[cl,'In Glom'] / length(glspots)
  map_prop_data[cl,"T Rel spot"] <- map_prop_data[cl,'In TI'] / length(tispots)
  map_prop_data[cl,"G Rel cl"] <- map_prop_data[cl,'In Glom'] / sum(stmerged@assays$pred.ec1.2@data[cl,])
  map_prop_data[cl,"T Rel cl"] <- map_prop_data[cl,'In TI'] / sum(stmerged@assays$pred.ec1.2@data[cl,])
}
map_prop_data[,'Rel spot ratio'] <- map_prop_data[,"G Rel spot"] / map_prop_data[,"T Rel spot"]

write.csv(map_prop_data,'map_prop_data_1.2.csv',quote = F,row.names = T)

map_prop_data <- as.data.frame(matrix(0,nrow=12,ncol=6,
                                      dimnames = list(as.character(0:11),
                                                      c('In DKD','In Ref',
                                                        'DKD Rel spot','Ref Rel spot',
                                                        'DKD Rel cl','Ref Rel cl'))))

dkdspots <- rownames(stmerged@meta.data[stmerged$Gloms == 'glom' & stmerged$condition=='DKD',])
refspots <- rownames(stmerged@meta.data[stmerged$Gloms == 'glom' & stmerged$condition=='Ref',])
for (cl in as.character(0:11)){
  map_prop_data[cl,'In DKD'] <- sum(stmerged@assays$pred.ec1.2@data[cl,dkdspots])
  map_prop_data[cl,'In Ref'] <- sum(stmerged@assays$pred.ec1.2@data[cl,refspots])
  map_prop_data[cl,"DKD Rel spot"] <- map_prop_data[cl,'In DKD'] / length(dkdspots)
  map_prop_data[cl,"Ref Rel spot"] <- map_prop_data[cl,'In Ref'] / length(refspots)
  map_prop_data[cl,"DKD Rel cl"] <- map_prop_data[cl,'In DKD'] / sum(stmerged@assays$pred.ec1.2@data[cl,])
  map_prop_data[cl,"Ref Rel cl"] <- map_prop_data[cl,'In Ref'] / sum(stmerged@assays$pred.ec1.2@data[cl,])
}
map_prop_data[,'Rel spot ratio'] <- map_prop_data[,"DKD Rel spot"] / map_prop_data[,"Ref Rel spot"]

write.csv(map_prop_data,'map_prop_dkd_data_1.2.csv',quote = F,row.names = T)

