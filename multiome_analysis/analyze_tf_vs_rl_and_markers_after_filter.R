library(rio)
library(qgraph)
library(ggplot2)
library(ggrepel)

gc_pec <- import_list('may2_table_filtered_peaks/Trajectory_EC-GC.pEC_May_02_2024.xlsx')
gc_iec <- import_list('may2_table_filtered_peaks/Trajectory_EC-GC.I-EC_May_02_2024.xlsx')

iec_tfs <- gc_iec[["SCmega"]][gc_iec[["SCmega"]]$GENE %in% c('SEMA6A','PLXNA2')&
                                abs(gc_iec[["SCmega"]]$correlation_TFxG) > .7,"TF"]

pec_tfs <- gc_pec[["SCmega"]][(gc_pec[["SCmega"]]$GENE %in% c('WNTB2')|
                                 gc_pec[["SCmega"]]$GENE %in% c('FN1')|
                                 #gc_pec[["SCmega"]]$GENE %in% c('JAG2')|
                                 #gc_pec[["SCmega"]]$GENE %in% c('JAG1')|
                                 gc_pec[["SCmega"]]$GENE %in% c('DLL4'))&
                                abs(gc_pec[["SCmega"]]$correlation_TFxG) > .7,]

iec_graph <- qgraph(gc_iec[["SCmega"]][,c("TF","GENE","correlation_TFxG")])
iec_cent <- centrality_auto(iec_graph)
pec_graph <- qgraph(gc_pec[["SCmega"]][,c("TF","GENE","correlation_TFxG")])
pec_cent <- centrality_auto(pec_graph)

iec_cent$node.centrality$Rank <- rank(-iec_cent$node.centrality$OutExpectedInfluence,ties.method = 'random')
iec_cent$node.centrality$label <- row.names(iec_cent$node.centrality)
iec_cent$node.centrality[!iec_cent$node.centrality$label %in% c('MEF2A','MEF2C','TRPS1'),'label'] <- NA
iec_cent$color <- 'gray'
iec_cent$node.centrality[iec_cent$node.centrality$label %in% c('MEF2A','MEF2C','TRPS1'),'color'] <- '#ff8c00'
pdf('centrality_iec.pdf',width = 4.5,height = 2)
ggplot(iec_cent$node.centrality,aes(y=OutExpectedInfluence,x=Rank,color=color,label=label))+
  geom_point(size=2.3)+
  geom_text_repel(min.segment.length = 30)+
  scale_color_manual(values=c('#ff8c00','gray'))+
  theme_classic()+
  theme(panel.border = element_rect(color = "black", fill=NA, linewidth = 1.5),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        legend.position="none")
dev.off()

pec_cent$node.centrality$Rank <- rank(-pec_cent$node.centrality$OutExpectedInfluence,ties.method = 'random')
pec_cent$node.centrality$label <- row.names(pec_cent$node.centrality)
pec_cent$node.centrality[!pec_cent$node.centrality$label %in% c('MEF2A','MEF2C','TRPS1'),'label'] <- NA
pec_cent$color <- 'gray'
pec_cent$node.centrality[pec_cent$node.centrality$label %in% c('MEF2A','MEF2C','TRPS1'),'color'] <- '#ff8c00'
pdf('centrality_pec.pdf',width = 4.5,height = 2)
ggplot(pec_cent$node.centrality,aes(y=OutExpectedInfluence,x=Rank,color=color,label=label))+
  geom_point(size=2.3)+
  geom_text_repel(min.segment.length = 30)+
  scale_color_manual(values=c('#ff8c00','gray'))+
  theme_classic()+
  theme(panel.border = element_rect(color = "black", fill=NA, linewidth = 1.5),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        legend.position="none")
dev.off()

mat <- gc_iec[["SCmega"]][,c("TF","GENE","correlation_TFxG")]
nrow(mat[mat$TF == 'GATA5',])


gc_vsmc <- import_list('may2_table_filtered_peaks/Trajectory_MC.VSMC-NS.I-VSMC_May_02_2024.xlsx')
vsmc_graph <- qgraph(gc_vsmc[["SCmega"]][,c("TF","GENE","correlation_TFxG")])
vsmc_cent <- centrality_auto(vsmc_graph)

write.table(gc_pec[["SCmega"]][gc_pec[["SCmega"]]$TF == 'MEF2C' & 
                                 gc_pec[["SCmega"]]$correlation_TFxG >.5 ,"GENE"],
            'mef2c_targets_ec-gc_pEC_trajectory.txt',quote = F,row.names = F,col.names = F)
write.table(gc_iec[["SCmega"]][gc_iec[["SCmega"]]$TF == 'MEF2C' & 
                                 gc_iec[["SCmega"]]$correlation_TFxG >.5 ,"GENE"],
            'mef2c_targets_ec-gc_iec_trajectory.txt',quote = F,row.names = F,col.names = F)

combined <- sort(unique(c(gc_pec[["SCmega"]][gc_pec[["SCmega"]]$TF == 'MEF2C' & 
                                               gc_pec[["SCmega"]]$correlation_TFxG >.5 ,"GENE"],
                          gc_iec[["SCmega"]][gc_iec[["SCmega"]]$TF == 'MEF2C' & 
                                               gc_iec[["SCmega"]]$correlation_TFxG >.5 ,"GENE"],
                          'MEF2C','MEF2D','TRPS1','SEMA6A','PLXNA2')))
write.table(combined,
            'mef2c_targets_combined_trajectory.txt',quote = F,row.names = F,col.names = F)

supptable <- rbind(gc_pec[["SCmega"]][gc_pec[["SCmega"]]$TF %in% c('MEF2A','MEF2C','TRPS1') & 
                                        gc_pec[["SCmega"]]$correlation_TFxG >.5 ,],
                   gc_iec[["SCmega"]][gc_iec[["SCmega"]]$TF %in% c('MEF2A','MEF2C','TRPS1') & 
                                        gc_iec[["SCmega"]]$correlation_TFxG >.5 ,])

write.csv(supptable,'supplemental_TF.csv',quote = F,row.names = F)
