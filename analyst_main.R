#author: David Lenci

library("tidyverse")
library("Seurat")
library("patchwork")
library("DESeq2")

read_data <- function(path){
  
  data <- readRDS(file=path)
  
  return(data)
  
}
data <- read_data("/usr4/bf527/dlenci/Documents/BF528-proj4-nogit/data/GSM2230760_seurat.rda")

data

get_markers <- function(seurat_data, markers) {
  
  data_markers <- FindAllMarkers(seurat_data, min.pct = 0.2) %>%
    group_by(cluster) %>%
    filter(gene%in%markers)

  
  return(data_markers)
  
}
markers <- c("SST", "INS", "COL1A1", "COL1A2", "FN1", "REG1B", "SPINK1", "CTRB2",
             "GCG", "KRT19", "ACER3", "RPS6KA5", "GC", "KRT19", "MMP7", "KRT7",
             "FCER1G", "IFI30", "TTR", "CHGB", "CRYBA2", "COL6A3", "ACP5")
data_markers <- get_markers(data, markers)

new.cluster.ids <- c("Delta", "Beta", "Stellate", "Acinar", "Alpha",
                     "Ductal", "Beta", "Oligodendrocyte", "Hepatocyte", 
                     "Stellate", "Ductal", "Acinar", "Macrophage")
names(new.cluster.ids) <- levels(data)
data <- RenameIdents(data, new.cluster.ids)


data_t <- RunUMAP(data, dims = 1:9)
png("umap.png")
DimPlot(data_t, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
dev.off()

png("hm.png")
DoHeatmap(data, features = data_markers$gene, label = F)
dev.off()