library(Seurat)

samples <- read.csv('../create_glom_ST/list_samples.csv')

samples$iu_id <- stringr::str_split(samples$x,'_',simplify = T)[,3]

kpmp <- rio::import('KPMP Biopsy LOG SHEET.xlsx')
kpmp <- kpmp[!is.na(kpmp$`IU ID`),]
rownames(kpmp) <- kpmp$`IU ID`

samples$study_id <- kpmp[samples$iu_id,'KPMP Study        ID']
samples$central_hub_id <- kpmp[samples$iu_id,'Central Hub        ID']

write.csv(samples[,3:4],'List_of_samples_for_Rosas.csv')
