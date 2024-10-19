#Script to generate prediction scores feature plots and piecharts plots

library(Seurat)
library(Matrix)
library(data.table)
library(dplyr)
library(gt)
library(igraph)
library(RColorBrewer)
library(cowplot)
library(rio)
library(stringr)
library(ggplot2)
library(reshape2)


#Table with celltypes labels and colors

coltable <- read.csv('../identify_glom_cell_types/colors_new_types.csv')
coltable <- rbind(coltable,
                  list(New.types = 'Other',Colors = '#DDDDDD'))
 
stmerged <- readRDS('gloms_manuscript_020824.RDS')

DefaultAssay(stmerged) <- 'pred.new.types'
Idents(stmerged) <- stmerged@meta.data$orig.ident



#Size factor varies for each sample
sizes=c(3,4,.6,
        3,1,2,
        3,1,1,
        1,1,.5,
        .4,.4,.4,
        .5,4,4,
        4,8,2,
        4,.5,.8,
        4,5,.4,
        .6,1,8,
        1,1.2,1,
        1,1,2)
unique(stmerged@meta.data$orig.ident)

for (i in 1:36){ #36
  sample <-  unique(stmerged@meta.data$orig.ident)[i]
  sz <- sizes[i]

  smp_split <- unlist(str_split(sample,'_'))
  img_folder <- smp_split[1]
  img_file <- list.files(paste0('../../spatial_samples/share_',sample),paste0(sample,'.tif'))
  
  spatial <- subset(stmerged,idents = sample)
  
  getPalette <- colorRampPalette(brewer.pal(9, "YlOrRd"))
  cell_types_all <- coltable$New.types
  
  img <- stringr::str_replace_all(sample,'-','.')
  if (ncol(spatial) > 1){
    # pdf(paste0('piecharts/',sample,'_features.pdf'),width=8,height = 4)
    # for (ctype in cell_types_all){
    #   p1 <- SpatialFeaturePlot(spatial,ctype,crop = F,pt.size.factor = sz,images = img)+ 
    #     scale_fill_gradientn(colors=rev(brewer.pal(n = 11, name = "Spectral")),
    #                          limits=c(0,0.25),
    #                          oob = scales::squish)
    #   plot(p1)
    # }
    # dev.off()
    
    
    
    #spript to generate an plot piecharts adapted from SPOTlight source
    metadata_ds <- as.data.frame(t(spatial@assays$pred.new.types@data))
    metadata_ds <- metadata_ds[,colnames(metadata_ds) != 'max']
    metadata_ds$Other <- rowSums(metadata_ds[,!colnames(metadata_ds) %in% coltable$New.types])
    
    rownames(coltable) <- coltable[[1]]
    coltable <- coltable[rownames(coltable) %in% colnames(metadata_ds),]
    metadata_ds <- metadata_ds[,coltable[[1]]]
    
    rownames(coltable) <- coltable[[1]]
    cols <- coltable[colnames(metadata_ds)[colSums(metadata_ds) > 0],][[2]]
    
    
    spatial_coord <- data.frame(spatial@images[[img]]@coordinates) %>%
      tibble::rownames_to_column("barcodeID") %>%
      dplyr::inner_join(metadata_ds %>% tibble::rownames_to_column("barcodeID"),
                        by = "barcodeID")
    
    
    
    img <- tiff::readTIFF(paste0('../../spatial_samples/share_',sample,'/',img_file))
    
    # Convert image to grob object
    img_grob <- grid::rasterGrob(img,
                                 interpolate = FALSE,
                                 width = grid::unit(1, "npc"),
                                 height = grid::unit(1, "npc"))
    
    ## Plot spatial scatterpie plot
    scatterpie_plt <- suppressMessages(
      ggplot2::ggplot() +
        ggplot2::annotation_custom(
          grob = img_grob,
          xmin = 0,
          xmax = ncol(img),
          ymin = 0,
          ymax = -nrow(img)) +
        scatterpie::geom_scatterpie(
          data = spatial_coord,
          ggplot2::aes(x = imagecol,
                       y = imagerow),
          cols = coltable[[1]],
          color = NA,
          alpha = 1,
          pie_scale = sz) +
        ggplot2::scale_y_reverse() +
        ggplot2::ylim(nrow(img), 0) +
        ggplot2::xlim(0, ncol(img)) +
        cowplot::theme_half_open(11, rel_small = 1) +
        ggplot2::theme_void() +
        ggplot2::coord_fixed(ratio = 1,
                             xlim = NULL,
                             ylim = NULL,
                             expand = TRUE,
                             clip = "on")+
        ggplot2::scale_fill_manual(values = cols))
    print(c(i,paste0('piecharts/',sample,'_pie.pdf')))
    pdf(paste0('piecharts/',sample,'_pie.pdf'),width=12,height = 8)
    plot(scatterpie_plt)
    dev.off()
  }  
  
}

##########
# level1 #
##########


### ORGANIZE COLTABLE

coltable[3,1] <- 'VSM/P'
coltable[5,1] <- 'EC'
coltable <- coltable[c(1,3,5,8),]
rownames(coltable) <- coltable[,1]

for (i in 1:36){ #33
  sample <-  unique(stmerged@meta.data$orig.ident)[i]
  sz <- sizes[i]
  
  smp_split <- unlist(str_split(sample,'_'))
  img_folder <- smp_split[1]
  img_file <- list.files(paste0('../../spatial_samples/share_',sample),paste0(sample,'.tif'))
  
  spatial <- subset(stmerged,idents = sample)
  
  getPalette <- colorRampPalette(brewer.pal(9, "YlOrRd"))
  cell_types_all <- coltable$New.types
  
  img <- stringr::str_replace_all(sample,'-','.')
  if (ncol(spatial) > 1){
     
    #spript to generate an plot piecharts adapted from SPOTlight source
    metadata_ds <- as.data.frame(t(spatial@assays$pred.new.types@data))
    metadata_ds <- metadata_ds[,colnames(metadata_ds) != 'max']
    
    ### AGGREGATE COLUMNS AND REMOVE LEVEL 2
    #metadata_ds$POD <- rowSums(metadata_ds[,c('POD','I-POD')])
    #metadata_ds$`I-POD` <- NULL
    metadata_ds$EC <- rowSums(metadata_ds[,c("EC-GC","EC-NOS","I-EC")])
    metadata_ds[,c("EC-GC","EC-NOS","I-EC")] <- NULL
    metadata_ds$`VSM/P` <- rowSums(metadata_ds[,c("MC","VSMC-NOS","I-VSMC")])
    metadata_ds[,c("MC","VSMC-NOS","I-VSMC")] <- NULL
    
    metadata_ds$Other <- rowSums(metadata_ds[,!colnames(metadata_ds) %in% coltable$New.types])
    
    rownames(coltable) <- coltable[[1]]
    coltable <- coltable[rownames(coltable) %in% colnames(metadata_ds),]
    metadata_ds <- metadata_ds[,coltable[[1]]]
    
    rownames(coltable) <- coltable[[1]]
    cols <- coltable[colnames(metadata_ds)[colSums(metadata_ds) > 0],][[2]]
    
    
    spatial_coord <- data.frame(spatial@images[[img]]@coordinates) %>%
      tibble::rownames_to_column("barcodeID") %>%
      dplyr::inner_join(metadata_ds %>% tibble::rownames_to_column("barcodeID"),
                        by = "barcodeID")
    
    
    
    img <- tiff::readTIFF(paste0('../../spatial_samples/share_',sample,'/',img_file))
    
    # Convert image to grob object
    img_grob <- grid::rasterGrob(img,
                                 interpolate = FALSE,
                                 width = grid::unit(1, "npc"),
                                 height = grid::unit(1, "npc"))
    
    ## Plot spatial scatterpie plot
    scatterpie_plt <- suppressMessages(
      ggplot2::ggplot() +
        ggplot2::annotation_custom(
          grob = img_grob,
          xmin = 0,
          xmax = ncol(img),
          ymin = 0,
          ymax = -nrow(img)) +
        scatterpie::geom_scatterpie(
          data = spatial_coord,
          ggplot2::aes(x = imagecol,
                       y = imagerow),
          cols = coltable[[1]],
          color = NA,
          alpha = 1,
          pie_scale = sz) +
        ggplot2::scale_y_reverse() +
        ggplot2::ylim(nrow(img), 0) +
        ggplot2::xlim(0, ncol(img)) +
        cowplot::theme_half_open(11, rel_small = 1) +
        ggplot2::theme_void() +
        ggplot2::coord_fixed(ratio = 1,
                             xlim = NULL,
                             ylim = NULL,
                             expand = TRUE,
                             clip = "on")+
        ggplot2::scale_fill_manual(values = cols))
    print(c(i,paste0('piecharts/',sample,'_level1_pie.pdf')))
    pdf(paste0('piecharts/',sample,'_level1_pie.pdf'),width=12,height = 8)
    plot(scatterpie_plt)
    dev.off()
  }  
  
}






