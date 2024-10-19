library(Seurat)
library(ggplot2)

samples <- list.files('.','*RDS')
samples <- samples[!grepl('merged',samples)]
samples <- stringr::str_remove(samples,'_ec.RDS')

spatial_list <- list()
for (smp in samples){
  spatial <- readRDS(paste0(smp,'_ec.RDS'))
  spatial@meta.data$orig.ident <- smp
  spatial_list[[smp]] <- spatial
  print(c(smp,length(colnames(spatial)),length(Idents(spatial))))
}

stmerged <- merge(spatial_list[[1]],spatial_list[2:41],
                  add.cell.ids=stringr::str_split(samples,'_',simplify = T)[,3])
print(c(length(colnames(stmerged)),length(Idents(stmerged))))

names(stmerged@images) <- stringr::str_replace_all(samples,'-','.')

#stmerged <- SCTransform(stmerged,assay = 'Spatial')
saveRDS(stmerged,'merged_41_after_ec.RDS')
#stmerged <- readRDS('merged_41_after_ec.RDS')

gloms <- readRDS('../../neighborhoods/gloms_ref_DKD_w_glom_index.RDS')

stmerged@meta.data$gloms <- 'no-glom'
stmerged@meta.data$condition <- 'NA'
smp <- unique(gloms$orig.ident)[1]
for (smp in unique(gloms$orig.ident)){
  df <- gloms@meta.data[gloms$orig.ident == smp,]
  pref <- unique(stringr::str_split(rownames(stmerged@meta.data[stmerged$orig.ident == smp,]),'_',simplify = T)[,1])
  suf <- unique(stringr::str_split(rownames(stmerged@meta.data[stmerged$orig.ident == smp,]),'_',simplify = T)[,3])
  rownames(df) <- paste0(pref[1],'_',stringr::str_split(rownames(df),'_',simplify = T)[,2],'_',suf[1])
  print(smp)
  print(any(!rownames(df)%in%rownames(stmerged@meta.data)))
  print(rownames(df)[1:6])
  stmerged@meta.data[rownames(df),'gloms'] <- 'glom'  
  stmerged@meta.data[stmerged@meta.data$orig.ident == smp,'condition'] <- unique(df$condition)  
  print(length(colnames(stmerged)) - length(Idents(stmerged)))
}
stmerged@meta.data[stmerged@meta.data$orig.ident == 'V10S21-388_XY03_20-0072','condition'] <- 'DKD'
stmerged@meta.data[stmerged@meta.data$orig.ident == 'V12N16-374_XY01_23-0108','condition'] <- 'DKD'


DefaultAssay(stmerged) <- 'pred.ec1.2'
# glomid <- factor(stmerged@meta.data$gloms,levels=as.character(0:10))
# names(glomid) <- colnames(glomid)
# stmerged@active.ident <- glomid
Idents(stmerged) <- stmerged@meta.data$gloms
DotPlot(stmerged,features = as.character(0:14),group.by = 'gloms')
DotPlot(stmerged,features = as.character(0:14),group.by = 'condition')


all_markers <- c('CD34','PECAM1','PTPRB','MEIS2','FLT1','EMCN',# EC
                 'HECW2','PLAT','ITGA8','EHD3','KDR','SOST',# 'EMCN',# EC-GC
                 'BTNL9','ADAMTS6','PALMD','AQP1','TM4SF1',
                 'VEGFC','CCDC3','CDH5','SERPINE2','FBLN5','CXCL12','SOX17',# EC-AEA
                 'MCTP1','SLC14A1','ENPP2','LYPD6B',#'BTNL9','ADAMTS6','PALMD','AQP1','TM4SF1', # EC-DVR
                 'CEACAM1','DNASE1L3','PLVAP','PITPNC1','GRB10','SLCO2A1','RAPGEF4',# EC-PTC
                 'GPM6A','EDIL3','TLL1','ZNF385D','NR2F2',#'CEACAM1','DNASE1L3','PLVAP', # EC-AVR
                 'MMRN1','CD36','TBX1','PKHD1L1','PROX1',# EC-LYM
                 'B2M','TMSB4X','TMSB10','HLA-B','IGFBP5',# dEC
                 'FTL','DHFR','SPP1','GPX3','PEBP1')# dEC-PTC

DefaultAssay(stmerged) <- 'SCT'
DotPlot(stmerged,features = all_markers,group.by = 'gloms')+
  theme(axis.text.x = element_text(angle = 90,vjust = .5))


DefaultAssay(stmerged) <- 'pred.vsmc1.2'

pdf('vln_glom_ti_res1.2.pdf')
VlnPlot(stmerged,features = as.character(0:10),group.by = 'gloms')
dev.off()

map_prop_data <- as.data.frame(matrix(0,nrow=11,ncol=6,
                                      dimnames = list(as.character(0:11),
                                                      c('In Glom','In TI',
                                                        'G Rel spot','T Rel spot',
                                                        'G Rel cl','T Rel cl'))))

glspots <- rownames(stmerged@meta.data[stmerged$gloms == 'glom',])
tispots <- rownames(stmerged@meta.data[stmerged$gloms != 'glom',])
for (cl in as.character(0:11)){
  map_prop_data[cl,'In Glom'] <- sum(stmerged@assays$pred.vsmc1.2@data[cl,glspots])
  map_prop_data[cl,'In TI'] <- sum(stmerged@assays$pred.vsmc1.2@data[cl,tispots])
  map_prop_data[cl,"G Rel spot"] <- map_prop_data[cl,'In Glom'] / length(glspots)
  map_prop_data[cl,"T Rel spot"] <- map_prop_data[cl,'In TI'] / length(tispots)
  map_prop_data[cl,"G Rel cl"] <- map_prop_data[cl,'In Glom'] / sum(stmerged@assays$pred.vsmc1.2@data[cl,])
  map_prop_data[cl,"T Rel cl"] <- map_prop_data[cl,'In TI'] / sum(stmerged@assays$pred.vsmc1.2@data[cl,])
}

write.csv(map_prop_data,'map_prop_data_1.2.csv',quote = F,row.names = T)

map_prop_data <- as.data.frame(matrix(0,nrow=9,ncol=6,
                                      dimnames = list(as.character(0:8),
                                                      c('In Glom','In TI',
                                                        'G Rel spot','T Rel spot',
                                                        'G Rel cl','T Rel cl'))))

glspots <- rownames(stmerged@meta.data[stmerged$gloms == 'glom',])
tispots <- rownames(stmerged@meta.data[stmerged$gloms != 'glom',])
for (cl in as.character(0:8)){
  map_prop_data[cl,'In Glom'] <- sum(stmerged@assays$pred.vsmc.8@data[cl,glspots])
  map_prop_data[cl,'In TI'] <- sum(stmerged@assays$pred.vsmc.8@data[cl,tispots])
  map_prop_data[cl,"G Rel spot"] <- map_prop_data[cl,'In Glom'] / length(glspots)
  map_prop_data[cl,"T Rel spot"] <- map_prop_data[cl,'In TI'] / length(tispots)
  map_prop_data[cl,"G Rel cl"] <- map_prop_data[cl,'In Glom'] / sum(stmerged@assays$pred.vsmc.8@data[cl,])
  map_prop_data[cl,"T Rel cl"] <- map_prop_data[cl,'In TI'] / sum(stmerged@assays$pred.vsmc.8@data[cl,])
}

