library(rio)
library(ggplot2)

data <- import('Eadon_Ferreira_mef2c_targets_ec-gc_pEC_trajectory.xlsx')

ec_ptc <- data[2:108,1:5]
colnames(ec_ptc) <- data[1,1:5]

ec_gc <- data[2:108,7:11]
colnames(ec_gc) <- data[1,7:11]

ec_aea <- data[2:108,13:17]
colnames(ec_aea) <- data[1,13:17]

ec_lym <- data[2:108,19:23]
colnames(ec_lym) <- data[1,19:23]


ec_gc$logFC_T2Di_T2D_Ec.GC <- as.numeric(ec_gc$logFC_T2Di_T2D_Ec.GC)
ec_gc$logFC_T2D_HC_Ec.GC <- as.numeric(ec_gc$logFC_T2D_HC_Ec.GC)
ec_gc$logFC_T2Di_HC_Ec.GC <- ec_gc$logFC_T2Di_T2D_Ec.GC + ec_gc$logFC_T2D_HC_Ec.GC
ec_gc$P.Value_T2D_HC <- as.numeric(ec_gc$P.Value_T2D_HC)
ec_gc$P.Value_T2Di_T2D <- as.numeric(ec_gc$P.Value_T2Di_T2D)

genes <- rownames(ec_gc[!is.na(ec_gc$logFC_T2D_HC_Ec.GC),])

dfplot <- ec_gc[genes,c("gene","logFC_T2D_HC_Ec.GC","logFC_T2Di_T2D_Ec.GC","logFC_T2Di_HC_Ec.GC")]
colnames(dfplot) <- c("gene","logFC_T2D_HC","logFC_T2Di_T2D","logFC_T2Di_HC")
geneorder <- dfplot$gene[order(dfplot$logFC_T2D_HC,decreasing = T)]
dfplot <- reshape2::melt(dfplot)
dfplot$gene <- factor(dfplot$gene,levels = geneorder)
dfplot$variable <- factor(dfplot$variable,levels=rev(c("logFC_T2D_HC","logFC_T2Di_T2D","logFC_T2Di_HC")))
pdf('heatmap_fold_change.pdf',height = 2,width = 8)
ggplot(dfplot,aes(x=gene,y=variable,fill = value))+
  geom_tile()+
  scale_fill_gradient2()+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90,vjust = .5))
dev.off()

dfcount <- ec_gc
dfcount$restore <- dfcount$logFC_T2D_HC_Ec.GC * dfcount$logFC_T2Di_HC_Ec.GC
#dfcount[is.na(dfcount$P.Value_T2D_HC),"P.Value_T2D_HC"] <- 1
dfcount <- dfcount[!is.na(dfcount$P.Value_T2D_HC),]
fishmat <- matrix(0,ncol=2,nrow = 2) 
fishmat[1,1] <- nrow(dfcount[dfcount$P.Value_T2D_HC < 0.01 & dfcount$restore > 0,])
fishmat[2,1] <- nrow(dfcount[dfcount$P.Value_T2D_HC > 0.01 & dfcount$restore > 0,])
fishmat[1,2] <- nrow(dfcount[dfcount$P.Value_T2D_HC < 0.01 & dfcount$restore < 0,])
fishmat[2,2] <- nrow(dfcount[dfcount$P.Value_T2D_HC > 0.01 & dfcount$restore < 0,])

test <- fisher.test(fishmat)

forest <- as.data.frame(matrix(0,ncol=3,nrow = 1,
                               dimnames = list(1,c('low','est','high'))))
forest[1,'est'] <- test$estimate#log(test$estimate)
forest[1,'low'] <- exp(log(test$estimate) - 1.96*sqrt(1/fishmat[1,1] + 1/fishmat[2,1] + 
                                                        1/fishmat[1,2] + 1/fishmat[2,2]))
forest[1,'high'] <- exp(log(test$estimate) + 1.96*sqrt(1/fishmat[1,1] + 1/fishmat[2,1] + 
                                                         1/fishmat[1,2] + 1/fishmat[2,2]))

forest$cell <- rownames(forest)
pdf('forest_sglt2_recover.pdf',width = 3,height = 3)
ggplot(data=forest, aes(x=cell, y=est, ymin=low, ymax=high)) +
  geom_pointrange(size = .7,linewidth = 1,shape = 2) + 
  geom_hline(yintercept=1, lty=2) +  # add a dotted line at x=1 after flip
  coord_flip() +  # flip coordinates (puts labels on y axis)
  #scale_x_discrete(limits=rownames(fisher_results))+
  xlab("Cell type") + ylab("OR (95% CI)")+
  # scale_color_manual(values = c(rep('#8DCFFF',10),rep('#dddddd',6),rep('#FFA829',5)),
  #                    name = element_blank(),
  #                    guide = guide_legend(reverse = TRUE))+
  theme_classic()+
  theme(panel.border = element_rect(color = "black", fill=NA, linewidth = 1.5),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10))
dev.off()




ec_gc <- import('mef2c_targets_combined_trajectory.xlsx',skip=1)

rownames(ec_gc) <- ec_gc$gene
genes <- rownames(ec_gc[!is.na(ec_gc$logFC_T2D_HC_EC.GC),])

# iec <- read.table('../multiome_analysis/mef2c_targets_ec-gc_iec_trajectory.txt')
# genes <- genes[genes %in% iec$V1]

ec_gc <- ec_gc[genes,]
genes <- rownames(ec_gc[ec_gc$adj.P.Value_T2D_HC < .05 & ec_gc$P.Value_T2Di_T2D...10 < .05,])


dfplot <- ec_gc[genes,c("gene","logFC_T2D_HC_EC.GC","logFC_T2Di_T2D_EC.GC","logFC_T2Di_HC_EC.GC")]
colnames(dfplot) <- c("gene","logFC_T2D_HC","logFC_T2Di_T2D","logFC_T2Di_HC")
geneorder <- dfplot$gene[order(dfplot$logFC_T2D_HC,decreasing = T)]
dfplot <- reshape2::melt(dfplot)
dfplot$gene <- factor(dfplot$gene,levels = geneorder)
dfplot$variable <- factor(dfplot$variable,levels=rev(c("logFC_T2D_HC","logFC_T2Di_T2D","logFC_T2Di_HC")))
pdf('heatmap_fold_change_new.pdf',height = 2,width = 8)
ggplot(dfplot,aes(x=gene,y=variable,fill = value))+
  geom_tile()+
  scale_fill_gradient2()+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90,vjust = .5))
dev.off()


#genes <- rownames(ec_gc[!is.na(ec_gc$logFC_T2D_HC_Ec.GC),])
#ec_gc[genes,'logFC_T2Di_HC_Ec.GC'] <- ec_gc[genes,'logFC_T2Di_T2D_Ec.GC'] + 
#  ec_gc[genes,'logFC_T2D_HC_Ec.GC']