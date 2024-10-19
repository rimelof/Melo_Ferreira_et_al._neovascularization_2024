library(pathfindR)
library(ggplot2)
library(Seurat)
library(EnhancedVolcano)



degglom <- read.csv('markers_iVSMC_vs_MC.csv')

pdf('volcano_iVSMC_vs_MC.pdf')
plot(  EnhancedVolcano(degglom,x = "avg_log2FC",y = "p_val_adj",lab = degglom$genes,
                       title = NULL,subtitle = NULL,legendPosition = 'none',
                       pCutoff = 0.05,FCcutoff = .5))
dev.off()
degglom <- degglom[c(6,2,5)]

gs <- 'mc'

setwd('pathways/')

for(subset in c('all','top200')){ 
  # degglom <- degglom[order(abs(degglom[,2])),]
  # degglom <- degglom[1:200,]
  for (pw in c('GO-BP','KEGG')){
    print(c(gs,subset,pw))
    #run_pathfindR(degglom,gene_sets = pw,max_gset_size = 500,pin_name_path = "STRING",)
    
    namepath = 'STRING'
    geneset = pw
    output_df <- run_pathfindR(degglom, output_dir = paste0(gs,'_',subset,'_',
                                                            namepath,'_',geneset,'/'),
                               gene_sets = geneset,
                               min_gset_size = 1,
                               max_gset_size = 500,
                               p_val_threshold=0.05,
                               pin_name_path = namepath,
                               plot_enrichment_chart=T,
                               n_processes = 1
    )
    
    enrichment_chart(output_df)
    RA_clustered <- cluster_enriched_terms(output_df, plot_dend = FALSE, plot_clusters_graph = FALSE,use_description=TRUE)
    
    RA_selected <- subset(RA_clustered,Cluster %in% 1:15)
    
    pdf(paste0(gs,'_',subset,'_',
               namepath,'_',geneset,'_path.pdf'),height = 20, width = 10)
    plot(enrichment_chart(RA_clustered, plot_by_cluster = TRUE))
    dev.off()
    
    saveRDS(output_df,paste0(gs,'_',subset,'_',
                             namepath,'_',geneset,'_outputdf.RDS'))
  }
  degglom <- degglom[order(abs(degglom[,2]),decreasing = T),]
  degglom <- degglom[1:200,]
}



setwd('../')


gs <- 'ecgc'

setwd('pathways/')

for(subset in c('all','top200')){ 
  # degglom <- degglom[order(abs(degglom[,2])),]
  # degglom <- degglom[1:200,]
  for (pw in c('GO-BP','KEGG')){
    namepath = 'STRING'
    geneset = pw
    output_df <- readRDS(paste0(gs,'_',subset,'_',
                                namepath,'_',geneset,'_outputdf.RDS'))
    
    pdf(paste0(gs,'_',subset,'_',
               namepath,'_',geneset,'simple_path.pdf'),height = 20, width = 10)
    plot(enrichment_chart(output_df))
    dev.off()
    
  }
}

setwd('../')


library(clusterProfiler)
#library(org.Hs.eg.db)

degglom <- read.csv('markers_iVSMC_vs_MC.csv')

# #genelist <- -log10(degglom$p_val_adj)
# #genelist <- abs(degglom$avg_log2FC)
# genelist <- degglom$avg_log2FC
# genelist <- degglom$signal2noise
# names(genelist) <- degglom$genes
# genelist <- sort(genelist,decreasing = T)
# gosetrich <- gseGO(
#   gene=genelist,
#   OrgDb = 'org.Hs.eg.db',
#   keyType = 'SYMBOL',
#   ont = 'BP',
#   pvalueCutoff = .05,
#   minGSSize = 10
# )
# 
# names(genelist) <- degglom$entrez
# keggsetrich <- gseKEGG(
#   gene=genelist,
#   organism = 'hsa',
#   keyType = 'ncbi-geneid',
#   pvalueCutoff = .05,
#   minGSSize = 10
# )
# 
# keggmodrich <- gseMKEGG(
#   gene=genelist,
#   organism = 'hsa',
#   keyType = 'ncbi-geneid',
#   pvalueCutoff = .05,
#   minGSSize = 10
# )

degglom <- degglom[degglom$p_val_adj < .05 &
                     abs(degglom$avg_log2FC) > .5, ]

# write.table(degglom[,c(6,2,5)],'mc_deg_4_ctpathwat.txt',sep = "\t",
#             col.names = F,row.names = F, quote = F)

goover <- enrichGO(
  gene=degglom$genes,
  OrgDb = 'org.Hs.eg.db',
  keyType = 'SYMBOL',
  ont = 'BP',
  pvalueCutoff = .05,
  minGSSize = 30
)


degglom$entrez <- bitr(degglom$genes,fromType = 'SYMBOL',toType = 'ENTREZID',
                       OrgDb = 'org.Hs.eg.db',drop = F)$ENTREZID

keggover <- enrichKEGG(
  gene=degglom$entrez,
  organism = 'hsa',
  keyType = 'ncbi-geneid',
  pvalueCutoff = .05,
  minGSSize = 30
)


write.csv(goover@result,'pathways_vsmc_go.csv',quote = F,row.names = F)
pdf('pathways_vsmc_go.pdf',width = 6,height = 8)
dotplot(goover,showCategory = 30,font.size=10)
dev.off()
write.csv(keggover@result,'pathways_vsmc_kegg.csv',quote = F,row.names = F)
pdf('pathways_vsmc_kegg.pdf',width = 6,height = 8)
dotplot(keggover,showCategory = 30,font.size=10)
dev.off()


# library(aPEAR)
# 
# p <- enrichmentNetwork(goover@result[1:30,], drawEllipses = TRUE, fontSize = 2.5)
# p




# output_df <- readRDS('degglom_0R.csv_all_STRING_GO-BP_outputdf.RDS')
# RA_clustered <- cluster_enriched_terms(output_df, plot_dend = FALSE, plot_clusters_graph = FALSE,use_description=TRUE)
# RA_selected <- subset(RA_clustered,Cluster %in% 1:15)
# pdf('degglom_0R.csv_all_STRING_GO-BP_path_final.pdf',height = 20, width = 10)
# plot(enrichment_chart(RA_clustered, plot_by_cluster = TRUE))+
#   scale_color_gradient(low = '#EEEEEE',high='#00930B')
# dev.off()
# 
# output_df <- readRDS('degglom_4R.csv_all_STRING_GO-BP_outputdf.RDS')
# RA_clustered <- cluster_enriched_terms(output_df, plot_dend = FALSE, plot_clusters_graph = FALSE,use_description=TRUE)
# RA_selected <- subset(RA_clustered,Cluster %in% 1:15)
# pdf('degglom_4R.csv_all_STRING_GO-BP_path_final.pdf',height = 20, width = 10)
# plot(enrichment_chart(RA_clustered, plot_by_cluster = TRUE))+
#   scale_color_gradient(low = '#EEEEEE',high='#FC0DBF')
# dev.off()
# 
# output_df <- readRDS('degglom_7R.csv_all_STRING_GO-BP_outputdf.RDS')
# RA_clustered <- cluster_enriched_terms(output_df, plot_dend = FALSE, plot_clusters_graph = FALSE,use_description=TRUE)
# RA_selected <- subset(RA_clustered,Cluster %in% 1:15)
# pdf('degglom_7R.csv_all_STRING_GO-BP_path_final.pdf',height = 10, width = 7)
# plot(enrichment_chart(RA_clustered, plot_by_cluster = TRUE))+
#   scale_color_gradient(low = '#EEEEEE',high='#8C0BD4')
# dev.off()
