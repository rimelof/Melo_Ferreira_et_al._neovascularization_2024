library(CellChat)
library(patchwork)
library(Seurat)
library(SeuratDisk)
options(stringsAsFactors = FALSE)

atlas <- LoadH5Seurat('../../atlas_v1_paper/Kidney_Healthy-Injury_Cell_Atlas_snCv3_Seurat_noNA_04222023.h5Seurat')


glommeta <- read.csv('../identify_glom_cell_types/glom_new_types_metadata_jan24.csv')
glommeta$seurat_clusters
rownames(glommeta) <- glommeta$X
atlas@meta.data$new_types <- atlas@meta.data$subclass.l2
atlas@meta.data[rownames(glommeta),'new_types'] <- as.character(glommeta$new_types)
atlas@meta.data[atlas@meta.data$new_types == 'I-POD','new_types'] <- 'POD'

Idents(atlas) <- 'new_types'


conditions = c('Ref','DKD')
modellist <- list(glom_only=c('POD','MC','EC-GC'),
                  glom_all=ctypes <- c('POD','MC','EC-GC','EC-NOS', 'I-EC', 'I-VSMC','VSMC-NOS'),
                  glom_extra=ctypes <- c('POD','MC','EC-GC','EC-NOS', 'I-EC', 'I-VSMC','VSMC-NOS',
                                         'MYOF','cycMYOF','FIB','M-FIB','dM-FIB','aFIB','dFIB',
                                         'B','PL','T','NKT','MAST','MAC-M2','cycMNP','MDC','cDC','pDC','ncMON','N'))


for (cond in conditions){
  for (model in names(modellist)){
    folder = paste0('sn_cellchat/',cond,'_',model)
    dir.create(folder)
    
    data.input = atlas@assays$RNA@data # normalized data matrix
    meta = atlas@meta.data # a dataframe with rownames containing cell mata data
    cell.use = rownames(meta)[meta$condition.l2 == cond & 
                                meta$new_types %in% modellist[[model]]] # extract the cell names from disease data
    
    # Prepare input data for CelChat analysis
    data.input = data.input[, cell.use]
    meta = meta[cell.use, ]
    # meta = data.frame(labels = meta$labels[cell.use], row.names = colnames(data.input)) # manually create a dataframe consisting of the cell labels
    unique(meta$new_types) # check the cell labels
    
    cellchat <- createCellChat(object = data.input, meta = meta, group.by = "new_types")
    
    cellchat <- addMeta(cellchat, meta = meta)
    cellchat <- setIdent(cellchat, ident.use = "new_types") # set "labels" as default cell identity
    levels(cellchat@idents) # show factor levels of the cell labels
    groupSize <- as.numeric(table(cellchat@idents)) # number of cells in each cell group
    
    CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data
    showDatabaseCategory(CellChatDB)
    
    # Show the structure of the database
    dplyr::glimpse(CellChatDB$interaction)
    
    # use a subset of CellChatDB for cell-cell communication analysis
    # CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling
    # use all CellChatDB for cell-cell communication analysis
    CellChatDB.use <- CellChatDB # simply use the default CellChatDB
    
    # set the used database in the object
    cellchat@DB <- CellChatDB.use
    
    # subset the expression data of signaling genes for saving computation cost
    cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
    future::plan("multicore", workers = 4) # do parallel
    
    cellchat <- identifyOverExpressedGenes(cellchat)
    cellchat <- identifyOverExpressedInteractions(cellchat)
    
    # project gene expression data onto PPI (Optional: when running it, USER should set `raw.use = FALSE` in the function `computeCommunProb()` in order to use the projected data)
    # cellchat <- projectData(cellchat, PPI.human)
    
    cellchat <- computeCommunProb(cellchat)
    cellchat <- filterCommunication(cellchat, min.cells = 10)
    
    df.net <- subsetCommunication(cellchat)
    write.csv(df.net,paste0(folder,'/dfnet.csv'),quote = F,row.names = F)
    df.path <- subsetCommunication(cellchat,slot.name = "netP")
    write.csv(df.net,paste0(folder,'/dfnetP.csv'),quote = F,row.names = F)
    
    cellchat <- computeCommunProbPathway(cellchat)
    df.net <- subsetCommunication(cellchat)
    write.csv(df.net,paste0(folder,'/dfpathnet.csv'),quote = F,row.names = F)
    df.path <- subsetCommunication(cellchat,slot.name = "netP")
    write.csv(df.net,paste0(folder,'/dfpathP.csv'),quote = F,row.names = F)
    
    cellchat <- aggregateNet(cellchat)
    groupSize <- as.numeric(table(cellchat@idents))
    par(mfrow = c(1,2), xpd=TRUE)
    pdf(paste0(folder,'/aggregated.pdf'))
    netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
    netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
    dev.off()
    
    mat <- cellchat@net$weight
    par(mfrow = c(1,7), xpd=TRUE)
    pdf(paste0(folder,'/per_type.pdf'),width = 20)
    for (i in 1:nrow(mat)) {
      mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
      mat2[i, ] <- mat[i, ]
      netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
    }
    dev.off()
    
    save(cellchat,file=paste0(folder,'/cellchat.RData'))
    
  }
}
