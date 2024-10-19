library(Seurat)
library(rio)
library(readxl)

list_sample <- list.files('../gloms/','*csv')
list_sample <- stringr::str_remove(list_sample,'.csv')

samples <- import('Processing_checklist_all.xlsx',range=cell_limits(c(2, 1), c(NA, NA)))
samples <- samples[samples$Sample %in% list_sample,]

spatial_list <- list()
for (r in 1:nrow(samples)){
  smp <- samples[r,"Sample"]
  spatial <- readRDS(paste0('../map_glom_new_types_pod_feb24/',smp,'_pod.RDS'))
  spatial@meta.data$orig.ident <- smp
  spatial@meta.data$condition <- samples[r,"Condition"]
  gloms <- read.csv(paste0('../gloms/',smp,'.csv'))
  rownames(gloms) <- gloms$Barcode
  spatial@meta.data$Gloms <- gloms[rownames(spatial@meta.data),'Gloms']
  print(smp)
  print(any(is.na(spatial$Gloms)))
  print(sort(unique(spatial$Gloms)))
  spatial_list[[samples[r,"Sample"]]] <- spatial
}

length(spatial_list)
stmerged <- merge(spatial_list[[1]],spatial_list[2:36],
                  add.cell.ids=stringr::str_split(samples$Sample,'_',simplify = T)[,3])
names(stmerged@images) <- stringr::str_replace_all(samples$Sample,'-','.')

stmerged <- SCTransform(stmerged,assay = 'Spatial')


kernel <- c()
edge <- c()
for (smp in names(stmerged@images)){
  coords <- stmerged@images[[smp]]@coordinates
  
  for (sp in 1:nrow(coords)){
    #sp <- 595
    spot <- rownames(coords)[sp]
    i <- coords[sp,"row"]
    j <- coords[sp,"col"]
    neigh <- rownames(coords[coords$row==i-1 & coords$col == j-1,])
    neigh <- c(neigh,rownames(coords[coords$row==i-1 & coords$col == j+1,]))
    neigh <- c(neigh,rownames(coords[coords$row==i & coords$col == j-2,]))
    neigh <- c(neigh,rownames(coords[coords$row==i & coords$col == j+2,]))
    neigh <- c(neigh,rownames(coords[coords$row==i+1 & coords$col == j-1,]))
    neigh <- c(neigh,rownames(coords[coords$row==i+1 & coords$col == j+1,]))
    
    if (length(neigh) == 6){
      kernel <- c(kernel,spot)
    } else {
      edge <- c(edge,spot)
    }
  }
  
}

stmerged <- subset(stmerged,cells = kernel)
SpatialDimPlot(stmerged,images = smp)




saveRDS(stmerged,'merged_manuscript_020824.RDS')
#stmerged <- readRDS('merged_manuscript_020824.RDS')
write.table(unique(stmerged@meta.data[,c("orig.ident","condition")]),
            '../samples_list/st_samples.csv',
            row.names = F,col.names = F,quote = F,sep = ',')

VlnPlot(stmerged,'AEBP1',group.by = 'orig.ident',split.by = 'condition')

Idents(stmerged) <- stmerged$Gloms
gloms <- subset(stmerged,idents='glom')
saveRDS(gloms,'gloms_manuscript_020824.RDS')
# gloms <- readRDS('gloms_manuscript_020824.RDS')

