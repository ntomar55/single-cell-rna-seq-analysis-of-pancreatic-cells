#---------------------------------
#Title: BF528 | Saxophone | Project 4 | Programmer
#Author: Nikita Tomar
#Date: 04/27/2022
#---------------------------------

library(Seurat)
library(tximport)
library(fishpond)
library(plyr)
library(dplyr)
library(Matrix)
library(RColorBrewer)
library('EnsDb.Hsapiens.v79')


#Part 1
#Load the UMI counts matrix - Salmon Quant Alevin file
files <- file.path("/usr4/bf527/ntomar/saxophone/data_p4/alevin_output/alevin/quants_mat.gz")
file.exists(files)

txi <- tximport(files, type="alevin")
txicts <- txi$counts

#convert the Ensembl ids to Gene Symbols


#creating Seurat Object
panc_cells <- CreateSeuratObject(counts = txi$counts, project = "panc_cells", min.cells = 3, min.features = 200)
panc_cells
head(panc_cells)
#rownames(txicts)
#mapping gene symbols
ens_ids <- rownames(txicts)
ens_ids <- sub("[.][0-9]*", "", ens_ids)
dim(ens_ids)
length(ens_ids)
rownames(txicts) <- ens_ids
geneIDs <- ensembldb::select(EnsDb.Hsapiens.v79, keys= ens_ids, keytype = "GENEID", columns = c("SYMBOL","GENEID"))
length(geneIDs$SYMBOL)
txicts <- txicts[rownames(txicts) %in% geneIDs$GENEID,]
rownames(txicts) <- geneIDs$SYMBOL

rownames(txicts)
geneIDs$SYMBOL

panc_cells[["percent.mt"]] <- PercentageFeatureSet(panc_cells, pattern = "^MT-") #calculates percentages of counts and does quality check
head(panc_cells@meta.data, 5)
x <- panc_cells[["percent.mt"]]

sum(panc_cells@meta.data$percent.mt)
#visualize QC metrics as violin plots
VlnPlot(panc_cells, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) #visualises as violin plot

plot1 <- FeatureScatter(panc_cells, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(panc_cells, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot3 <- FeatureScatter(panc_cells, feature1 = "nFeature_RNA", feature2 = "percent.mt")
plot1 + plot2 + plot3

#filtering low quality cells 
panc_cells <- subset(panc_cells, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

#determine number of cells after filtering out low-quality cells
ncol(x=panc_cells)


#Check plots after filtering out low-quality cells
VlnPlot(panc_cells, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(panc_cells, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(panc_cells, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot3 <- FeatureScatter(panc_cells, feature1 = "nFeature_RNA", feature2 = "percent.mt")
plot1 + plot2 + plot3

######## PART 2 ######################

#normalize the data
panc_cells <- NormalizeData(panc_cells)

#identify highly variable features 
panc_cells <- FindVariableFeatures(panc_cells, selection.method = "vst", nfeatures = 2000)

#Top 10 most highly variable genes
top10 <- head(VariableFeatures(panc_cells), 10)
top10

plot4<- VariableFeaturePlot(panc_cells)
plot5 <- LabelPoints(plot = plot4, points = top10, repel = TRUE, xnudge = 0, ynudge = 0)
plot4 + plot5

#scaling the data
all.genes <- rownames(panc_cells)
panc_cells <- ScaleData(panc_cells, features = all.genes)

#perform linear dimensional reduction
panc_cells <- RunPCA(panc_cells, features = VariableFeatures(object = panc_cells))

#Method 1 - Look at PCA results
print(panc_cells[["pca"]], dims = 1:5, nfeatures = 5)
DimPlot(panc_cells, reduction = "pca")

#Method 2 - Look at PCA results
VizDimLoadings(panc_cells, dims = 1:2, reduction = "pca")

#Method 3 - Look at PCA results
DimPlot(panc_cells, reduction = "pca")


#Method 4 - DimHeatmap
DimHeatmap(panc_cells, dims = 1, cells = 500, balanced = TRUE)


#testing for association between observed and latent variables
panc_cells <- JackStraw(panc_cells, num.replicate = 100)
panc_cells <- ScoreJackStraw(panc_cells, dims = 1:20)
JackStrawPlot(panc_cells, dims = 1:20)
ElbowPlot(panc_cells)

#Part 3
#Determine clusters
panc_cells <- FindNeighbors(panc_cells, dims=1:10)
panc_cells <- FindClusters(panc_cells, resolution = 0.5)

#Cluster IDs for first 5 cells
head(Idents(panc_cells), 5)

#Non-Linear Dimensional Reduction
panc_cells <- RunUMAP(panc_cells, dims=1:10)
DimPlot(panc_cells, reduction = "umap")

#Save to RDS
saveRDS(panc_cells, "/usr4/bf527/ntomar/saxophone/data_p4/programmer/panc_cells.rds")

#Create pie chart to display relative proportions of count in each cluster
#stores the cluster number and count of cell number
cluster_counts <- as.data.frame(table(Idents(panc_cells)))
cluster_counts <- paste(cluster_counts$Var1, cluster_counts$Freq, sep=" - ")

#stores the cluster number and the relative proportions
cluster_proportions <- as.data.frame(prop.table(table(Idents(panc_cells))))
cluster_proportions$Freq <- cluster_proportions$Freq %>% `*`(100) %>% round(3)
cluster_proportions$Percents <- cluster_proportions$Freq
cluster_proportions$Percents <- paste(cluster_proportions$Percents, "%", sep="")

#plot the piechart
cols = colorRampPalette(brewer.pal(8, "Dark2"))(14)
pie(cluster_proportions$Freq, labels=cluster_proportions$Percents, col=cols, main = "Relative Proportions of Cell Numbers\nFor the Identified Clusters")
legend(1.45, 1, legend=cluster_counts, cex = 0.8, fill = cols, bty="n", title="Clusters")
