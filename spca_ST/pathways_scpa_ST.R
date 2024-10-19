library(SCPA)
library(msigdbr)
library(Seurat)
library(tidyverse)
library(rio)
library(ggplot2)
library(ggrepel)
library(stringr)

st_neigh <- readRDS('../neighborhoods/gloms_neighborhood_res2_020524.RDS')

path_gobp <- msigdbr('Homo sapiens','C5','GO:BP') %>% format_pathways()
path_kegg <- msigdbr('Homo sapiens','C2','CP:KEGG') %>% format_pathways()



comparisons <- list(ec =  c("1", "20", "18"),
                    pod = c("15", "11" , "3", "14", "19"))

path_list <- list()
for (comp in names(comparisons)){
  scpa <- compare_seurat(st_neigh,
                         group1 = 'seurat_clusters',
                         group1_population = comparisons[[comp]],
                         #group2 = "condition.l2", 
                         #group2_population = c("DKD", "Ref"),
                         pathways = path_gobp,
                         assay = 'SCT',
                         downsample = 500)
  path_list[[paste0('glom_',comp,'_gobp')]] <- scpa
  
  scpa <- compare_seurat(st_neigh,
                         group1 = 'seurat_clusters',
                         group1_population = comparisons[[comp]],
                         #group2 = "condition.l2", 
                         #group2_population = c("DKD", "Ref"),
                         pathways = path_kegg,
                         assay = 'SCT',
                         downsample = 500)
  path_list[[paste0('glom_',comp,'_kegg')]] <- scpa
}
