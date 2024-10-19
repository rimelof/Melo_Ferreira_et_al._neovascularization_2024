library(SCPA)
library(msigdbr)
library(Seurat)
library(tidyverse)
library(rio)
library(ggplot2)
library(ggrepel)
library(stringr)
#library(rentrez)

glom <- readRDS('../identify_glom_cell_types/glom_reclustered.RDS')
DimPlot(glom)

path_gobp <- msigdbr('Homo sapiens','C5','GO:BP') %>% format_pathways()
path_kegg <- msigdbr('Homo sapiens','C2','CP:KEGG') %>% format_pathways()

path_list <- list()

scpa_glom_gobp <- compare_seurat(glom,
                            group1 = "condition.l2", 
                            group1_population = c("DKD", "Ref"),
                            pathways = path_gobp,
                            downsample = 500)
path_list[['glom_gobp']] <- scpa_glom_gobp

scpa_glom_kegg <- compare_seurat(glom,
                                 group1 = "condition.l2", 
                                 group1_population = c("DKD", "Ref"),
                                 pathways = path_kegg,
                                 downsample = 500)
path_list[['glom_kegg']] <- scpa_glom_kegg

comparisons <- list(ec =  c("I-EC", "EC-GC",'EC-NOS'),
                    mc = c("I-VSMC","MC","VSMC-NOS"),
                    pod = c("I-POD","POD"))

for (comp in names(comparisons)){
  scpa <- compare_seurat(glom,
                         group1 = 'new_types',
                         group1_population = comparisons[[comp]],
                         #group2 = "condition.l2", 
                         #group2_population = c("DKD", "Ref"),
                         pathways = path_gobp,
                         downsample = 500)
  path_list[[paste0('glom_',comp,'_gobp')]] <- scpa
  
  scpa <- compare_seurat(glom,
                         group1 = 'new_types',
                         group1_population = comparisons[[comp]],
                         #group2 = "condition.l2", 
                         #group2_population = c("DKD", "Ref"),
                         pathways = path_kegg,
                         downsample = 500)
  path_list[[paste0('glom_',comp,'_kegg')]] <- scpa
}

saveRDS(path_list,'path_list.RDS')
#path_list <- readRDS('path_list.RDS')
export(path_list,'path_list.xlsx')
path_list <- import_list('path_list_annot.xlsx')

paths <- path_list[[3]][!is.na(path_list[[3]]$Labels),1]
path_list[[3]] <- path_list[[3]][order(path_list[[3]]$qval,decreasing = T),]
path_list[[3]]$Rank <- as.numeric(rownames(path_list[[3]]))
#path_list[[3]] <- path_list[[3]][rev(order(path_list[[3]]$Rank)),]
colnames(path_list[[3]])[5] <- 'color'
path_list[[3]]$Labels <- ifelse(is.na(path_list[[3]]$color),NA,
                                str_to_title(str_replace_all(str_remove(path_list[[3]]$Pathway,'GOBP_'),'_',' ')))
#path_list[[3]]$color <- ifelse(!is.na(path_list[[3]]$Labels),"#FF8C00","#333333")
pdf('pathway_dkd_ec.pdf',width = 4,height = 5.5)
ggplot(path_list[[3]],aes(x=Rank,y=qval,color=color,label=Labels))+
  geom_point()+
  geom_text_repel(max.overlaps = 50,max.time = .5,color='black',size=4,min.segment.length = 0,nudge_x = 2)+
  scale_color_manual(values=c("#FF8C00","#333333"))+
  theme_classic()+
  theme(panel.border = element_rect(color = "black", fill=NA, linewidth = 1.5),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        legend.position = 'none')+
  ylab('Q value')
dev.off()

paths <- path_list[[5]][!is.na(path_list[[5]]$Labels),1]
path_list[[5]] <- path_list[[5]][order(path_list[[5]]$qval,decreasing = T),]
path_list[[5]]$Rank <- as.numeric(rownames(path_list[[5]]))
#path_list[[5]] <- path_list[[5]][rev(order(path_list[[5]]$Rank)),]
colnames(path_list[[5]])[5] <- 'color'
path_list[[5]]$Labels <- ifelse(is.na(path_list[[5]]$color),NA,
                                str_to_title(str_replace_all(str_remove(path_list[[5]]$Pathway,'GOBP_'),'_',' ')))
#path_list[[5]]$color <- ifelse(!is.na(path_list[[5]]$Labels),"#FF8C00","#333333")
pdf('pathway_dkd_mc.pdf',width = 4,height = 5.5)
ggplot(path_list[[5]],aes(x=Rank,y=qval,color=color,label=Labels))+
  geom_point()+
  geom_text_repel(max.overlaps = 50,max.time = .5,color='black',size=4,min.segment.length = 0,nudge_x = 2)+
  scale_color_manual(values=c("#483D8B","#333333"))+
  theme_classic()+
  theme(panel.border = element_rect(color = "black", fill=NA, linewidth = 1.5),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        legend.position = 'none')+
  ylab('Q value')
dev.off()

paths <- path_list[[7]][!is.na(path_list[[7]]$Labels),1]
path_list[[7]] <- path_list[[7]][order(path_list[[7]]$qval,decreasing = T),]
path_list[[7]]$Rank <- as.numeric(rownames(path_list[[7]]))
#path_list[[7]] <- path_list[[7]][rev(order(path_list[[7]]$Rank)),]
colnames(path_list[[7]])[6] <- 'color'
path_list[[7]]$Labels <- ifelse(is.na(path_list[[7]]$color),NA,
                                str_to_title(str_replace_all(str_remove(path_list[[7]]$Pathway,'GOBP_'),'_',' ')))
#path_list[[7]]$color <- ifelse(!is.na(path_list[[7]]$Labels),"#FF8C00","#333333")
pdf('pathway_dkd_pod.pdf',width = 4,height = 5.5)
ggplot(path_list[[7]],aes(x=Rank,y=qval,color=color,label=Labels))+
  geom_point()+
  geom_text_repel(max.overlaps = 50,max.time = .5,color='black',size=4,min.segment.length = 0,nudge_x = 2)+
  scale_color_manual(values=c("#DB7295","#333333"))+
  theme_classic()+
  theme(panel.border = element_rect(color = "black", fill=NA, linewidth = 1.5),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        legend.position = 'none')+
  ylab('Q value')
dev.off()

stglom <- readRDS('../create_glom_ST/gloms_manuscript_012224.RDS')
pathst <- list()
pathst[['spatial_gobp']] <- compare_seurat(stglom,
                                           assay = "SCT",
                                           group1 = "condition", 
                                           group1_population = c("DKD", "Ref"),
                                           pathways = path_gobp,
                                           downsample = 1824)
pathst[['spatial_kegg']] <- compare_seurat(stglom,
                                           assay = "SCT",
                                           group1 = "condition", 
                                           group1_population = c("DKD", "Ref"),
                                           pathways = path_kegg,
                                           downsample = 1824)
saveRDS(pathst,'path_st.RDS')
#pathst <- readRDS('path_st.RDS')
export(pathst,'path_st.xlsx')

path_l1_dkd <- list()
for (comp in unique(glom$subclass.l1)){
  scpa <- compare_seurat(glom,
                         group1 = "condition.l2", 
                         group1_population = c("DKD", "Ref"),
                         group2 = "subclass.l1", 
                         group2_population = comp,
                         pathways = path_gobp,
                         downsample = 500)
  path_l1_dkd[[paste0('glom_',comp,'_gobp')]] <- scpa
  
  scpa <- compare_seurat(glom,
                         group1 = "condition.l2", 
                         group1_population = c("DKD", "Ref"),
                         group2 = "subclass.l1", 
                         group2_population = comp,
                         pathways = path_kegg,
                         downsample = 500)
  path_l1_dkd[[paste0('glom_',comp,'_kegg')]] <- scpa
}

saveRDS(path_l1_dkd,'path_l1_DKD.RDS')
#path_list <- readRDS('path_list.RDS')
names(path_l1_dkd)[5:6] <- c("glom_VSM-P_gobp","glom_VSM-P_kegg")
export(path_l1_dkd,'path_l1_DKD.xlsx')
path_l1_dkd <- import_list('path_l1_DKD_annot.xlsx')



paths <- path_l1_dkd[[1]][!is.na(path_l1_dkd[[1]]$Labels),1]
path_l1_dkd[[1]] <- path_l1_dkd[[1]][order(path_l1_dkd[[1]]$qval,decreasing = T),]
path_l1_dkd[[1]]$Rank <- as.numeric(rownames(path_l1_dkd[[1]]))
#path_l1_dkd[[1]] <- path_l1_dkd[[1]][rev(order(path_l1_dkd[[1]]$Rank)),]
colnames(path_l1_dkd[[1]])[6] <- 'color'
path_l1_dkd[[1]]$Labels <- ifelse(is.na(path_l1_dkd[[1]]$color),NA,
                                str_to_title(str_replace_all(str_remove(path_l1_dkd[[1]]$Pathway,'GOBP_'),'_',' ')))
#path_l1_dkd[[1]]$color <- ifelse(!is.na(path_l1_dkd[[1]]$Labels),"#FF8C00","#333333")
pdf('pathway_all_ec.pdf',width = 4,height = 5.5)
ggplot(path_l1_dkd[[1]],aes(x=Rank,y=qval,color=color,label=Labels))+
  geom_point()+
  geom_text_repel(max.overlaps = 50,max.time = .5,color='black',size=4,min.segment.length = 0,nudge_x = 2)+
  scale_color_manual(values=c("#FF8C00","#333333"))+
  theme_classic()+
  theme(panel.border = element_rect(color = "black", fill=NA, linewidth = 1.5),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        legend.position = 'none')+
  ylab('Q value')
dev.off()

paths <- path_l1_dkd[[5]][!is.na(path_l1_dkd[[5]]$Labels),1]
path_l1_dkd[[5]] <- path_l1_dkd[[5]][order(path_l1_dkd[[5]]$qval,decreasing = T),]
path_l1_dkd[[5]]$Rank <- as.numeric(rownames(path_l1_dkd[[5]]))
#path_l1_dkd[[5]] <- path_l1_dkd[[5]][rev(order(path_l1_dkd[[5]]$Rank)),]
colnames(path_l1_dkd[[5]])[6] <- 'color'
path_l1_dkd[[5]]$Labels <- ifelse(is.na(path_l1_dkd[[5]]$color),NA,
                                str_to_title(str_replace_all(str_remove(path_l1_dkd[[5]]$Pathway,'GOBP_'),'_',' ')))
#path_l1_dkd[[5]]$color <- ifelse(!is.na(path_l1_dkd[[5]]$Labels),"#FF8C00","#333333")
pdf('pathway_all_mc.pdf',width = 4,height = 5.5)
ggplot(path_l1_dkd[[5]],aes(x=Rank,y=qval,color=color,label=Labels))+
  geom_point()+
  geom_text_repel(max.overlaps = 50,max.time = .5,color='black',size=4,min.segment.length = 0,nudge_x = 2)+
  scale_color_manual(values=c("#483D8B","#333333"))+
  theme_classic()+
  theme(panel.border = element_rect(color = "black", fill=NA, linewidth = 1.5),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        legend.position = 'none')+
  ylab('Q value')
dev.off()

paths <- path_l1_dkd[[3]][!is.na(path_l1_dkd[[3]]$Labels),1]
path_l1_dkd[[3]] <- path_l1_dkd[[3]][order(path_l1_dkd[[3]]$qval,decreasing = T),]
path_l1_dkd[[3]]$Rank <- as.numeric(rownames(path_l1_dkd[[3]]))
#path_l1_dkd[[3]] <- path_l1_dkd[[3]][rev(order(path_l1_dkd[[3]]$Rank)),]
colnames(path_l1_dkd[[3]])[6] <- 'color'
path_l1_dkd[[3]]$Labels <- ifelse(is.na(path_l1_dkd[[3]]$color),NA,
                                str_to_title(str_replace_all(str_remove(path_l1_dkd[[3]]$Pathway,'GOBP_'),'_',' ')))
#path_l1_dkd[[3]]$color <- ifelse(!is.na(path_l1_dkd[[3]]$Labels),"#FF8C00","#333333")
pdf('pathway_all_pod.pdf',width = 4,height = 5.5)
ggplot(path_l1_dkd[[3]],aes(x=Rank,y=qval,color=color,label=Labels))+
  geom_point()+
  geom_text_repel(max.overlaps = 50,max.time = .5,color='black',size=4,min.segment.length = 0,nudge_x = 2)+
  scale_color_manual(values=c("#DB7295","#333333"))+
  theme_classic()+
  theme(panel.border = element_rect(color = "black", fill=NA, linewidth = 1.5),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        legend.position = 'none')+
  ylab('Q value')
dev.off()









# search_pubmed <- function(keyword1, keyword2) {
#   query <- paste("(",keyword1, ") AND (", keyword2,")")
#   search_results <- entrez_search(db="pubmed", term=query)
#   return(search_results)
# }
# 
# 
# for (i in 3:length(path_list)){
#   path_list[[i]]$num.results <- 0
#   path_list[[i]]$pmid.results <- ''
#   for (p in 1:nrow(path_list[[i]])){
#     pw <- stringr::str_split(path_list[[i]][p,1],'_',simplify = T)
#     pw <- paste(pw[1,2:ncol(pw)],collapse = ' ')
#     
#     result_a<- search_pubmed("Diabetic Kidney Disease",
#                              pw)
#     path_list[[i]][p,'num.results'] <- result_a$count
#     path_list[[i]][p,'pmid.results'] <- paste(result_a$ids,collapse = ' ') 
#   }
# }
# 
