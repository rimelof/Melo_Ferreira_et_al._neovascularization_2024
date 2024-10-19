library(Seurat)
library(sp)

xenium <- LoadXenium('../../xenium/Xenium_1_02122024/output-XETG00126__0010200__f59__20240214__210015/')
xenium <- subset(xenium,subset = nCount_Xenium > 0)

selections <- read.csv('../glom_selections_xenium/f59_coordinates..csv',
                       row.names = NULL,skip = 2)

coords <- as.data.frame(xenium@images$fov@boundaries$centroids@coords)
rownames(coords) <- xenium@images$fov@boundaries$centroids@cells

in_selec1 <- point.in.polygon(coords$x,coords$y,
                              selections[selections$Selection == 'Selection 1','X'],
                              selections[selections$Selection == 'Selection 1','Y'])

cells <- rownames(coords)[in_selec1==1]

glom <- subset(xenium,cells = cells)
ImageDimPlot(glom)
