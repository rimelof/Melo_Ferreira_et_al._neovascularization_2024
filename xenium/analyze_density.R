library(ggplot)
library(stringr)

files <- list.files('../glom_selections_xenium/','coordinates..csv')
files <- files[!grepl('Selection',files)]

samples <- str_split(files,'_',simplify = T)[,1]
names(files) <- samples

reads <- read.csv('../glom_selections_xenium/3723_Selection_1_transcript_stats..csv',skip = 3)
counts <- data.frame(row.names = reads$Gene)
density <- data.frame(row.names = reads$Gene)
for (smp in names(files)){
  #smp <- names(files)[1]
  coords <- read.csv(paste0('../glom_selections_xenium/',files[smp]),skip = 2)
  for (sel in unique(coords$Selection)){
    #sel <- unique(coords$Selection)[1]
    reads <- read.csv(paste0('../glom_selections_xenium/',smp,'_',str_replace(sel,' ','_'),'_transcript_stats..csv'),
                      skip = 3)
    rownames(reads) <- reads$Gene
    counts <- cbind(counts, reads[,'Count',drop=F])
    colnames(counts)[ncol(counts)] <- paste0(smp,'_G',str_split(sel,' ',simplify = T)[,2])
    density <- cbind(density, reads[,"Density..µm..2.",drop=F])
    colnames(density)[ncol(density)] <- paste0(smp,'_G',str_split(sel,' ',simplify = T)[,2])
    
  }
  
}
#3775_Selection_102_transcript_stats..csv
#40775_Selection_21_transcript_stats..csv


density2 <- reshape2::melt(as.matrix(density))
density2$sample <- str_split(density2$Var2,'_',simplify = T)[,1]
density2$condition <- ifelse(density2$sample %in% c("3723","3775","f59"),
                            'Ref','DKD')
density2$condition <- factor(density2$condition,levels = c('Ref','DKD'))

densiplot <- density2[density2$Var1 %in% c('NPHS2','PECAM1','TAGLN','SEMA6A','PLXNA2'),]
densiplot$Var1 <- factor(densiplot$Var1,levels=c('NPHS2','PECAM1','TAGLN','SEMA6A','PLXNA2'))
pdf('violin_densities.pdf',width = 4.5,height = 3)
ggplot(densiplot,aes(x=Var1,y=value,fill=condition))+
  geom_violin(scale = 'width')+
  scale_fill_manual(values = c('#8DCFFF','#FFA829'))+
  theme_classic()
dev.off()

densiratio <- as.data.frame(t(density[c('PECAM1','SEMA6A'),]))
densiratio$ratio <- densiratio$SEMA6A / densiratio$PECAM1
densiratio$sample <- str_split(rownames(densiratio),'_',simplify = T)[,1]
densiratio$condition <- ifelse(densiratio$sample %in% c("3723","3775","f59"),
                             'Ref','DKD')
densiratio$condition <- factor(densiratio$condition,levels = c('Ref','DKD'))
densiratio$x <- 'SEMA6A / PECAM1 ' 

pdf('violin_ratio_sema6a.pdf',width = 2,height = 3)
ggplot(densiratio,aes(x=x,y=ratio,fill=condition))+
  geom_violin(scale = 'width')+
  scale_fill_manual(values = c('#8DCFFF','#FFA829'))+
  theme_classic()
dev.off()

densiratio <- as.data.frame(t(density[c('PECAM1','TAGLN','PLXNA2'),]))
densiratio$ratio <- densiratio$PLXNA2 / (densiratio$PECAM1 +densiratio$TAGLN)
densiratio$sample <- str_split(rownames(densiratio),'_',simplify = T)[,1]
densiratio$condition <- ifelse(densiratio$sample %in% c("3723","3775","f59"),
                               'Ref','DKD')
densiratio$condition <- factor(densiratio$condition,levels = c('Ref','DKD'))
densiratio$x <- 'PLXNA2 / (PECAM1 + TAGLN)' 

pdf('violin_ratio_plxna2.pdf',width = 2,height = 3)
ggplot(densiratio,aes(x=x,y=ratio,fill=condition))+
  geom_violin(scale = 'width')+
  scale_fill_manual(values = c('#8DCFFF','#FFA829'))+
  theme_classic()
dev.off()


densiplot <- density2[density2$Var1 %in% c('PDGFR','PDGFRB'),]
pdf('violin_densities_PDGFR.pdf',width = 2,height = 3)
ggplot(densiplot,aes(x=Var1,y=value,fill=condition))+
  geom_violin(scale = 'width')+
  scale_fill_manual(values = c('#8DCFFF','#FFA829'))+
  theme_classic()
dev.off()
