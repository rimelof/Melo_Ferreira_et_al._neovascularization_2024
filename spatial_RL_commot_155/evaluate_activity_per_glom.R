library(Seurat)
library(heatmap3)
library(PCAtools)
library(Boruta)

gloms <- readRDS('../neighborhoods/gloms_neighborhood_res2_020524.RDS')

glom_ind=0
gloms@meta.data$glom_index <- NA
all_spots <- rownames(gloms@meta.data)
while (length(all_spots) > 0){
  spot <- all_spots[1]
  print(spot)
  smp <- gloms@meta.data[spot,'orig.ident']
  img <- stringr::str_replace_all(smp,'-','.')
  smpspots <- gloms@images[[img]]@coordinates
  glom <- spot
  neigh <- spot
  while(length(neigh) > 0){
    i <- smpspots[neigh[1],"row"]
    j <- smpspots[neigh[1],"col"]
    neigh <- c(neigh,
               rownames(smpspots[(smpspots$row == i-1 & smpspots$col == j-1) |
                                   (smpspots$row == i-1 & smpspots$col == j+1) |
                                   (smpspots$row == i & smpspots$col == j-2) |
                                   (smpspots$row == i & smpspots$col == j+2) |
                                   (smpspots$row == i+1 & smpspots$col == j-1) |
                                   (smpspots$row == i+1 & smpspots$col == j+1),]))
    neigh <- neigh[!neigh %in% glom]
    glom <- c(glom,neigh)
  }
  glom_ind <- glom_ind+1
  gloms@meta.data[glom,"glom_index"] <- paste0('G',glom_ind)
  all_spots <-  all_spots[!all_spots %in% glom]
}

unique(gloms@meta.data$glom_index)




activity <- c()
for (gl in unique(gloms@meta.data$glom_index)){
  sample <- unique(gloms@meta.data[gloms@meta.data$glom_index == gl,"orig.ident"])
  print(c(length(sample),nrow(gloms@meta.data[gloms$glom_index == gl,])))
  sum_receiver <- read.csv(paste0('simple_path_pairs_per_glom/',sample,'/obsm/commot-cellchat-sum-receiver.csv'),
                           row.names = 1)
  rnames <- rownames(gloms@meta.data[gloms$glom_index == gl,])
  spots <- stringr::str_split(rnames,'_',simplify = T)[,2]
  sum_receiver <- sum_receiver[spots,]
  rownames(sum_receiver) <- rnames
  sum_receiver$glom_index <- gl
  activity <- rbind(activity,sum_receiver)
}

# write.csv(activity,'simple_path_pairs_per_glom/activity_per_glom.csv',row.names = T,col.names = T,quote = F)
# #activity <- read.csv('simple_path_pairs_per_glom/activity_per_glom.csv',row.names = 1)
# 
# act_mat <- activity
# act_mat <- aggregate(.~glom_index,act_mat,sum)
# row.names(act_mat) <- act_mat$glom_index
# act_mat$glom_index <- NULL
# act_mat <- as.matrix(act_mat)
# log_act <- log(act_mat+1)
# act_mat <- act_mat[,apply(act_mat,2,sd)>0]
# log_act <- log_act[,apply(log_act,2,sd)>0]
# 
# glom_cond <- unique(gloms@meta.data[,c("glom_index","condition")])
# row.names(glom_cond) <- glom_cond$glom_index
# glom_cond$glom_index <- NULL
# glom_cond <- glom_cond[rownames(act_mat),]
# 
# pdf('glom_heatmap.pdf')
# hm <- heatmap3(act_mat,RowSideLabs = glom_cond)
# dev.off()
# 
# pdf('glom_heatmap_log.pdf')
# hm <- heatmap3(log_act,RowSideLabs = as.factor(glom_cond))
# dev.off()
# 
# 
# act_mat <- activity
# act_mat <- aggregate(.~glom_index,act_mat,sum)
# row.names(act_mat) <- act_mat$glom_index
# act_mat$glom_index <- NULL
# act_mat <- as.matrix(act_mat)
# 
# glom_cond <- unique(gloms@meta.data[,c("glom_index","condition")])
# row.names(glom_cond) <- glom_cond$glom_index
# #glom_cond$glom_index <- NULL
# glom_cond <- glom_cond[rownames(act_mat),]
# 
# p <- pca(t(act_mat), metadata = glom_cond, removeVar = 0.1)
# 
# biplot(p,
#        #lab = paste0(p$metadata$Age, ' años'),
#        colby = 'condition',
#        hline = 0, vline = 0,
#        legendPosition = 'right')
# 
# 
# act_mat <- activity
# act_mat$glom_index <- NULL
# act_mat <- as.matrix(act_mat)
# 
# glom_cond <- gloms@meta.data[,c("glom_index","condition","seurat_clusters")]
# glom_cond <- glom_cond[rownames(act_mat),]
# 
# p <- pca(t(act_mat), metadata = glom_cond, removeVar = 0.1)
# 
# biplot(p,
#        #lab = paste0(p$metadata$Age, ' años'),
#        colby = 'seurat_clusters',
#        hline = 0, vline = 0,
#        legendPosition = 'right')


act_mat <- activity
act_mat <- aggregate(.~glom_index,act_mat,sum)
row.names(act_mat) <- act_mat$glom_index
act_mat$glom_index <- NULL

glom_cond <- unique(gloms@meta.data[,c("glom_index","condition")])
row.names(glom_cond) <- glom_cond$glom_index
act_mat$condition <- factor(glom_cond[rownames(act_mat),"condition"])

boruta_output <- Boruta(condition ~ ., data = act_mat, doTrace = 2)
boruta_signif <- names(boruta_output$finalDecision[boruta_output$finalDecision %in% c("Confirmed", "Tentative")])  # collect Confirmed and Tentative variables
boruta_signif
boruta_signif <- names(boruta_output$finalDecision[boruta_output$finalDecision %in% c("Confirmed")])  # collect Confirmed and Tentative variables
boruta_signif
write.csv(boruta_signif,'boruta_signif_02192024.csv',quote = F)

rl_pairs <- boruta_signif[1:38]
rl_pairs <- stringr::str_split(stringr::str_remove(rl_pairs,'r\\.'),'\\.',simplify = T)
write.csv(rl_pairs,'boruta_signif_rlpairs_02202024.csv',quote = F)


sort(boruta_signif[95:128])

act_mat <- activity[,c(boruta_signif,'glom_index')]
act_mat <- aggregate(.~glom_index,act_mat,sum)
row.names(act_mat) <- act_mat$glom_index
act_mat$glom_index <- NULL
act_mat <- reshape2::melt(as.matrix(act_mat))
act_mat$condition <- gloms@meta.data[act_mat$Var1,'condition']
act_mat <- reshape2::dcast(act_mat,Var2 ~ condition,fun.aggregate = mean)
act_mat$FC <- act_mat$DKD / act_mat$Ref

act_mat <- activity[,c(boruta_signif,'glom_index')]
act_mat <- aggregate(.~glom_index,act_mat,sum)
row.names(act_mat) <- act_mat$glom_index
act_mat$glom_index <- NULL
act_mat <- as.matrix(act_mat)
#act_mat <- act_mat[,apply(act_mat,2,sd)>0]

glom_cond <- unique(gloms@meta.data[,c("glom_index","condition")])
row.names(glom_cond) <- glom_cond$glom_index
glom_cond$glom_index <- NULL
glom_cond <- glom_cond[rownames(act_mat),]

pdf('glom_heatmap_boruta.pdf')
hm <- heatmap3(act_mat,RowSideLabs = glom_cond)
dev.off()


act_mat <- activity[,c(boruta_signif,'glom_index')]
act_mat <- aggregate(.~glom_index,act_mat,sum)
rownames(act_mat) <- act_mat$glom_index
act_mat$glom_index <- NULL
act_mat <- as.matrix(act_mat)

glom_cond <- unique(gloms@meta.data[,c("glom_index","condition")])
row.names(glom_cond) <- glom_cond$glom_index
#glom_cond$glom_index <- NULL
glom_cond <- glom_cond[rownames(act_mat),]

p <- pca(t(act_mat), metadata = glom_cond, removeVar = 0.1)

biplot(p,
       #lab = paste0(p$metadata$Age, ' años'),
       colby = 'condition',
       hline = 0, vline = 0,
       legendPosition = 'right')
