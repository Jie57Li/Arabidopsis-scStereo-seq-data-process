##### Sub cluster of Mesophyll (R Script for Figure 3)

# 1. library
#-------------------------------------------------------------------------------
# library(reshape2)
# library(ggplot2)
# library(factoextra)
library(Seurat)
library(dplyr)
library(patchwork)
library(openxlsx)
setwd("Projects/Fig3")

# 2. load data and save the subtype
#-------------------------------------------------------------------------------
obj <- readRDS("../RDS/Run.rds")
Mesophyll <- subset(obj,idents = "Mesophyll")
# saveRDS(Mesophyll,"../RDS/Mesophyll_subtype_pre.rds")

# 3. Cluster-pre 
#-------------------------------------------------------------------------------
## pre_prosses (var-genes)
Mesophyll <- FindVariableFeatures(Mesophyll, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(Mesophyll), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(Mesophyll)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1
plot2

all.genes <- rownames(Mesophyll)
Mesophyll <- ScaleData(Mesophyll, features = all.genes)
Mesophyll <- RunPCA(Mesophyll, features = VariableFeatures(object = Mesophyll))
# Examine and visualize PCA results a few different ways
print(Mesophyll[["pca"]], dims = 1:5, nfeatures = 5)
DimPlot(Mesophyll, reduction = "pca")
DimHeatmap(Mesophyll, dims = 1, cells = 500, balanced = TRUE)
ElbowPlot(Mesophyll)
Mesophyll <- FindNeighbors(Mesophyll, dims = 1:10)
Mesophyll <- FindClusters(Mesophyll, resolution = 0.5)
# If you haven't installed UMAP, you can do so via reticulate::py_install(packages =
# 'umap-learn')
Mesophyll <- RunUMAP(Mesophyll, dims = 1:10)
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
pre1 <- DimPlot(Mesophyll, reduction = "umap",pt.size = 2)
pre2 <- DimPlot(Mesophyll, reduction = "umap",pt.size = 2)+ NoLegend()
pre3 <- DimPlot(Mesophyll, reduction = "umap",group.by = "cell_type",pt.size = 2)
pre4 <- DimPlot(Mesophyll, reduction = "umap",group.by = "cell_type",pt.size = 2) + NoLegend()


# 4. DEGs-PCA/UMAP 
#-------------------------------------------------------------------------------
Mesophyll@active.ident <- factor(Mesophyll@meta.data$cell_type)
names(Mesophyll@active.ident) <- rownames(Mesophyll@meta.data)
DEGs_all <- FindMarkers(Mesophyll,ident.1 = "palisad",min.pct = 0,logfc.threshold = 0)
DEGs_all$gene <- rownames(DEGs_all)

genes1 <- DEGs_all[DEGs_all$p_val<0.05 & DEGs_all$avg_log2FC < -0.35,"gene"]
genes2 <- DEGs_all[DEGs_all$p_val<0.05 & DEGs_all$avg_log2FC > 0.35,"gene"]
Mesophyll2 <- Mesophyll
Mesophyll2 <- RunPCA(Mesophyll2, features = c(genes1,genes2))
# Examine and visualize PCA results a few different ways
print(Mesophyll2[["pca"]], dims = 1:4, nfeatures = 5)
DimPlot(Mesophyll2, reduction = "pca")
DimHeatmap(Mesophyll2, dims = 1:10, cells = 300, balanced = TRUE)
ElbowPlot(Mesophyll2)
Mesophyll2 <- FindNeighbors(Mesophyll2, dims =1:3)
Mesophyll2 <- FindClusters(Mesophyll2, resolution = 1.9)
table(Mesophyll2@meta.data[,c("seurat_clusters","cell_type")])


Mesophyll2 <- RunUMAP(Mesophyll2, dims = 1:3)
# Mesophyll2 <- RunTSNE(Mesophyll2, dims = 1:4)
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
p1 <- DimPlot(Mesophyll2, reduction = "umap",pt.size = 2) + NoLegend()
p2 <- DimPlot(Mesophyll2, reduction = "umap",pt.size = 2)
DimPlot(Mesophyll2, reduction = "umap",group.by = "cell_type",pt.size = 2)

# 6. Heatmap
#-------------------------------------------------------------------------------
library(reshape2)
library(ggplot2)
library(factoextra)
library(Seurat)
library(dplyr)
library(patchwork)
library(ggplot2)
library(cowplot)
library(reshape2)
library(factoextra)
library(ComplexHeatmap)
library(grid)
library(circlize)
library("GetoptLong")

mat <- as.data.frame(Mesophyll2@assays$RNA@data)
# cell type
types <- obj@meta.data[,c("orig.ident","cell_type","pca_type")]
celltype <- c("palisad","sponge")

names <- "Palisade_Sponge"
gene1 <- as.data.frame(genes1)
gene2 <- as.data.frame(genes2)
gene1$type <- "Palisade"
gene2$type <- "Sponge"
colnames(gene1) <- c("gene","type")
colnames(gene2) <- c("gene","type")
degs <- rbind(gene1,gene2)

## 3. the data of heatmap
reslut <- "Mesophyll/00.DEG-heatmap/"
if ( !dir.exists(reslut)) {
  print("0. Create specified directory...")
  dir.create(reslut)
}
celltype_S <- celltype
cells <- c(rownames(types[types$cell_type == celltype[1],]),
           rownames(types[types$cell_type == celltype[2],]))
data_fc <- expm1(mat[degs$gene,cells])
data_exp <- data_fc + 1

cell1 <- rownames(types[types$cell_type == celltype[1],])
cell2 <- rownames(types[types$cell_type == celltype[2],])
data_fc[,cell1] <- data_exp[,cell1]/rowMeans(data_exp[,cell2])
data_fc[,cell2] <- data_exp[,cell2]/rowMeans(data_exp[,cell1])

Groups<- c()
for (z in colnames(data_fc)) {
  if (z %in% rownames(types[types$cell_type == "palisad",])) Groups <- c(Groups, 'Palisade')
  if (z %in% rownames(types[types$cell_type == "sponge",])) Groups <- c(Groups, 'Sponge')
}

collist <- list(Groups=c("Palisade"="#DD322E","Sponge"="#40845C"))
ha <-HeatmapAnnotation(df = data.frame(Type = Groups),col = collist)
f2 <- colorRamp2(c(-2,0,2), c("#265297","#EEEEEE","#DD322E"))
df_p <- log2(data_fc)
df_p$gene <- rownames(df_p)
data_ed <- merge(degs,df_p,by = "gene")
data_ed$type <- factor(data_ed$type,levels = c("Palisade","Sponge"))
rownames(data_ed) <- data_ed$gene
h <- Heatmap(as.matrix(data_ed[,c(-1,-2)]), name ="-log2 (FC)", column_title = " ",row_title = " ",
             row_names_side = "left",row_dend_side = "left",col = f2,
             cluster_columns = F,cluster_rows = F,split = data_ed$type,
             show_column_names = F,show_row_names = F,
             row_names_gp = gpar(fontsize = 10,font=3),
             top_annotation = ha)
h      
