library(Seurat)
library(rio)

sn <- readRDS('mesangial_snRNAApr_01_2024.RDS')
Idents(sn) <- sn$predsnRNA0.5

markers <- FindMarkers(sn,ident.1='pEC',ident.2=c('EC-GC','I-EC','EC-NOS'))
markers$gene <- rownames(markers)
markers <- markers[abs(markers$avg_log2FC) > 2,]
markers <- markers[order(markers$avg_log2FC),]

rl <- read.csv('../sn_RL/sn_cellchat_pEC/comparison_glom_all/selected_rl.csv',header = F)

gc_nos_pec <- import_list('Trajectory_EC-GC.EC-NOS.pEC_Mar_27_2024.xlsx')
gc_nos_pec_rl <- rbind(gc_nos_pec[[1]][gc_nos_pec[[1]]$GENE %in% c(rl$V1) |
                                      gc_nos_pec[[1]]$TF %in% c(rl$V1),],
                    gc_nos_pec[[2]][gc_nos_pec[[2]]$GENE %in% c(rl$V1) |
                                      gc_nos_pec[[2]]$TF %in% c(rl$V1),],
                    gc_nos_pec[[3]][gc_nos_pec[[3]]$GENE %in% c(rl$V1) |
                                      gc_nos_pec[[3]]$TF %in% c(rl$V1),],
                    gc_nos_pec[[4]][gc_nos_pec[[4]]$GENE %in% c(rl$V1) |
                                      gc_nos_pec[[4]]$TF %in% c(rl$V1),])
gc_nos_pec_rl_simple <- unique(gc_nos_pec_rl[,c("TIME","TF","GENE","correlation_TFxG")])
gc_nos_pec_mk <- rbind(gc_nos_pec[[1]][gc_nos_pec[[1]]$GENE %in% c(markers$gene) |
                                      gc_nos_pec[[1]]$TF %in% c(markers$gene),],
                    gc_nos_pec[[2]][gc_nos_pec[[2]]$GENE %in% c(markers$gene) |
                                      gc_nos_pec[[2]]$TF %in% c(markers$gene),],
                    gc_nos_pec[[3]][gc_nos_pec[[3]]$GENE %in% c(markers$gene) |
                                      gc_nos_pec[[3]]$TF %in% c(markers$gene),],
                    gc_nos_pec[[4]][gc_nos_pec[[4]]$GENE %in% c(markers$gene) |
                                      gc_nos_pec[[4]]$TF %in% c(markers$gene),])
gc_nos_pec_mk_simple <- unique(gc_nos_pec_mk[,c("TIME","TF","GENE","correlation_TFxG")])
# gc_nos_pec <- rbind(gc_nos_pec[[1]][gc_nos_pec[[1]]$GENE %in% c(markers$gene,rl$V1) |
#                                          gc_nos_pec[[1]]$TF %in% c(markers$gene,rl$V1),],
#                        gc_nos_pec[[2]][gc_nos_pec[[2]]$GENE %in% c(markers$gene,rl$V1) |
#                                          gc_nos_pec[[2]]$TF %in% c(markers$gene,rl$V1),],
#                        gc_nos_pec[[3]][gc_nos_pec[[3]]$GENE %in% c(markers$gene,rl$V1) |
#                                          gc_nos_pec[[3]]$TF %in% c(markers$gene,rl$V1),],
#                        gc_nos_pec[[4]][gc_nos_pec[[4]]$GENE %in% c(markers$gene,rl$V1) |
#                                          gc_nos_pec[[4]]$TF %in% c(markers$gene,rl$V1),])


gc_iec <- import_list('Trajectory_EC-GC.I-EC_Apr_01_2024.xlsx')
gc_iec_rl <- rbind(gc_iec[[1]][gc_iec[[1]]$GENE %in% c(rl$V1) |
                                         gc_iec[[1]]$TF %in% c(rl$V1),],
                       gc_iec[[2]][gc_iec[[2]]$GENE %in% c(rl$V1) |
                                         gc_iec[[2]]$TF %in% c(rl$V1),],
                       gc_iec[[3]][gc_iec[[3]]$GENE %in% c(rl$V1) |
                                         gc_iec[[3]]$TF %in% c(rl$V1),],
                       gc_iec[[4]][gc_iec[[4]]$GENE %in% c(rl$V1) |
                                         gc_iec[[4]]$TF %in% c(rl$V1),])
gc_iec_rl_simple <- unique(gc_iec_rl[,c("TIME","TF","GENE","correlation_TFxG")])
gc_iec_mk <- rbind(gc_iec[[1]][gc_iec[[1]]$GENE %in% c(markers$gene) |
                                         gc_iec[[1]]$TF %in% c(markers$gene),],
                       gc_iec[[2]][gc_iec[[2]]$GENE %in% c(markers$gene) |
                                         gc_iec[[2]]$TF %in% c(markers$gene),],
                       gc_iec[[3]][gc_iec[[3]]$GENE %in% c(markers$gene) |
                                         gc_iec[[3]]$TF %in% c(markers$gene),],
                       gc_iec[[4]][gc_iec[[4]]$GENE %in% c(markers$gene) |
                                         gc_iec[[4]]$TF %in% c(markers$gene),])
gc_iec_mk_simple <- unique(gc_iec_mk[,c("TIME","TF","GENE","correlation_TFxG")])



gc_pec <- import_list('Trajectory_EC-GC.pEC_Apr_01_2024.xlsx')
gc_pec_rl <- rbind(gc_pec[[1]][gc_pec[[1]]$GENE %in% c(rl$V1) |
                                         gc_pec[[1]]$TF %in% c(rl$V1),],
                       gc_pec[[2]][gc_pec[[2]]$GENE %in% c(rl$V1) |
                                         gc_pec[[2]]$TF %in% c(rl$V1),],
                       gc_pec[[3]][gc_pec[[3]]$GENE %in% c(rl$V1) |
                                         gc_pec[[3]]$TF %in% c(rl$V1),],
                       gc_pec[[4]][gc_pec[[4]]$GENE %in% c(rl$V1) |
                                         gc_pec[[4]]$TF %in% c(rl$V1),])
gc_pec_rl_simple <- unique(gc_pec_rl[,c("TIME","TF","GENE","correlation_TFxG")])
gc_pec_mk <- rbind(gc_pec[[1]][gc_pec[[1]]$GENE %in% c(markers$gene) |
                                         gc_pec[[1]]$TF %in% c(markers$gene),],
                       gc_pec[[2]][gc_pec[[2]]$GENE %in% c(markers$gene) |
                                         gc_pec[[2]]$TF %in% c(markers$gene),],
                       gc_pec[[3]][gc_pec[[3]]$GENE %in% c(markers$gene) |
                                         gc_pec[[3]]$TF %in% c(markers$gene),],
                       gc_pec[[4]][gc_pec[[4]]$GENE %in% c(markers$gene) |
                                         gc_pec[[4]]$TF %in% c(markers$gene),])
gc_pec_mk_simple <- unique(gc_pec_mk[,c("TIME","TF","GENE","correlation_TFxG")])

mk_vsmc <- FindMarkers(sn,ident.1='MC',ident.2=c('I-VSMC','VSMC-NOS'))
mk_vsmc$gene <- rownames(mk_vsmc)
mk_vsmc <- mk_vsmc[abs(mk_vsmc$avg_log2FC) > 2,]
mk_vsmc <- mk_vsmc[order(mk_vsmc$avg_log2FC),]

mc_vsmcnos_ivsmc <- import_list('Trajectory_MC.VSMC-NOS.I-VSMC_Mar_27_2024.xlsx')
mc_vsmcnos_ivsmc_rl <- rbind(mc_vsmcnos_ivsmc[[1]][mc_vsmcnos_ivsmc[[1]]$GENE %in% c(rl$V1) |
                                 mc_vsmcnos_ivsmc[[1]]$TF %in% c(rl$V1),],
                   mc_vsmcnos_ivsmc[[2]][mc_vsmcnos_ivsmc[[2]]$GENE %in% c(rl$V1) |
                                 mc_vsmcnos_ivsmc[[2]]$TF %in% c(rl$V1),],
                   mc_vsmcnos_ivsmc[[3]][mc_vsmcnos_ivsmc[[3]]$GENE %in% c(rl$V1) |
                                 mc_vsmcnos_ivsmc[[3]]$TF %in% c(rl$V1),],
                   mc_vsmcnos_ivsmc[[4]][mc_vsmcnos_ivsmc[[4]]$GENE %in% c(rl$V1) |
                                 mc_vsmcnos_ivsmc[[4]]$TF %in% c(rl$V1),])
mc_vsmcnos_ivsmc_rl_simple <- unique(mc_vsmcnos_ivsmc_rl[,c("TIME","TF","GENE","correlation_TFxG")])
mc_vsmcnos_ivsmc_mk <- rbind(mc_vsmcnos_ivsmc[[1]][mc_vsmcnos_ivsmc[[1]]$GENE %in% c(mk_vsmc$gene) |
                                 mc_vsmcnos_ivsmc[[1]]$TF %in% c(mk_vsmc$gene),],
                   mc_vsmcnos_ivsmc[[2]][mc_vsmcnos_ivsmc[[2]]$GENE %in% c(mk_vsmc$gene) |
                                 mc_vsmcnos_ivsmc[[2]]$TF %in% c(mk_vsmc$gene),],
                   mc_vsmcnos_ivsmc[[3]][mc_vsmcnos_ivsmc[[3]]$GENE %in% c(mk_vsmc$gene) |
                                 mc_vsmcnos_ivsmc[[3]]$TF %in% c(mk_vsmc$gene),],
                   mc_vsmcnos_ivsmc[[4]][mc_vsmcnos_ivsmc[[4]]$GENE %in% c(mk_vsmc$gene) |
                                 mc_vsmcnos_ivsmc[[4]]$TF %in% c(mk_vsmc$gene),])
mc_vsmcnos_ivsmc_mk_simple <- unique(mc_vsmcnos_ivsmc_mk[,c("TIME","TF","GENE","correlation_TFxG")])
