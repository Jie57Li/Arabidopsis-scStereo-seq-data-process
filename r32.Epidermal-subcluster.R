##### Sub cluster of Epidermal (R Script for Figure 3)

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

## 2. load data and save the subtype
obj <- readRDS("../RDS/Run-PCA-871cells.rds")
Epidermal <- subset(obj,idents = "Epidermal")
# saveRDS(sub_obj,"../a0.RDS/Epidermal_subtype_pre.rds")

# 4. Cluster-pre 
#-------------------------------------------------------------------------------
## pre_prosses (var-genes)
Epidermal <- FindVariableFeatures(Epidermal, selection.method = "vst", nfeatures = 2000)

## Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(Epidermal), 10)

## plot variable features with and without labels
plot1 <- VariableFeaturePlot(Epidermal)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1
plot2

all.genes <- rownames(Epidermal)
Epidermal <- ScaleData(Epidermal, features = all.genes)
Epidermal <- RunPCA(Epidermal, features = VariableFeatures(object = Epidermal))
# Examine and visualize PCA results a few different ways
print(Epidermal[["pca"]], dims = 1:5, nfeatures = 5)
DimPlot(Epidermal, reduction = "pca")
DimHeatmap(Epidermal, dims = 1, cells = 500, balanced = TRUE)
ElbowPlot(Epidermal)
Epidermal <- FindNeighbors(Epidermal, dims = 1:10)
Epidermal <- FindClusters(Epidermal, resolution = 0.5)
# If you haven't installed UMAP, you can do so via reticulate::py_install(packages =
# 'umap-learn')
Epidermal <- RunUMAP(Epidermal, dims = 1:10)
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(Epidermal, reduction = "umap")
pre1 <- DimPlot(Epidermal, reduction = "umap",group.by = "cell_type",pt.size = 2)
# pdf("05.Subcluster/Epidermal/01.Cluster/01.Pre-Cluster.pdf",width = 4.5,height = 4.1)
# pre1
# dev.off()

# 5. DEGs-PCA/UMAP 
#-------------------------------------------------------------------------------
Epidermal@active.ident <- factor(Epidermal@meta.data$cell_type)
names(Epidermal@active.ident) <- rownames(Epidermal@meta.data)
## DEGs between lowerepidermal and uperepidermal
DEGs_all <- FindMarkers(Epidermal,ident.1 = "lowerepidermal",min.pct = 0,logfc.threshold = 0)
DEGs_all$gene <- rownames(DEGs_all)
# write.xlsx(DEGs_all,"05.Subcluster/Epidermal/Epidermal-DEGs-all.xlsx")

# load the degs
#genes1 <- read.table("05.Subcluster/Epidermal/Up_Low_Epidermal_genes.txt",header = T)
genes1 <- DEGs_all[DEGs_all$p_val<0.05 & DEGs_all$avg_log2FC < -0.25,"gene"]
genes2 <- DEGs_all[DEGs_all$p_val<0.05 & DEGs_all$avg_log2FC > 0.25,"gene"]
# write.table(as.data.frame(genes1),"../Fig3/05.Subcluster/Epidermal/Upperepidermal.txt",row.names = F,quote = F,sep = "\t")
# write.table(as.data.frame(genes2),"../Fig3/05.Subcluster/Epidermal/Lowerepidermal.txt",row.names = F,quote = F,sep = "\t")
Epidermal2 <- Epidermal
Epidermal2 <- RunPCA(Epidermal2, features = c(genes1,genes2))
# Examine and visualize PCA results a few different ways
print(Epidermal2[["pca"]], dims = 1:5, nfeatures = 5)
DimPlot(Epidermal2, reduction = "pca")
DimHeatmap(Epidermal2, dims = 1, cells = 500, balanced = TRUE)
ElbowPlot(Epidermal2)
Epidermal2 <- FindNeighbors(Epidermal2, dims = 1:5)
Epidermal2 <- FindClusters(Epidermal2, resolution = 2)
table(Epidermal2@meta.data[,c("seurat_clusters","cell_type")])
Epidermal2 <- RunUMAP(Epidermal2, dims = 1:5)
# p1 <- DimPlot(Epidermal2, reduction = "umap",pt.size = 2) + NoLegend()
# p2 <- DimPlot(Epidermal2, reduction = "umap",pt.size = 2)
col <- c("#F39700","#739AFF")
p3 <- DimPlot(Epidermal2, reduction = "umap",group.by = "cell_type",pt.size = 2,cols = col) + NoLegend()
p4 <- DimPlot(Epidermal2, reduction = "umap",group.by = "cell_type",pt.size = 2,cols = col) 
save.image("Epidermal/Epidermal-subcluster.RData")
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

mat <- as.data.frame(Epidermal2@assays$RNA@data)
# cell type
types <- Epidermal2@meta.data[,c("orig.ident","cell_type","pca_type")]
celltype <- c("upperepidermal","lowerepidermal")
names <- "Upper_Lowerepidermal"
gene1 <- as.data.frame(genes1)
gene2 <- as.data.frame(genes2)
gene1$type <- "Upper"
gene2$type <- "Lower"
colnames(gene1) <- c("gene","type")
colnames(gene2) <- c("gene","type")
degs <- rbind(gene1,gene2)

# 3. the data of heatmap
reslut <- "0.DEG-heatmap/"
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
  if (z %in% rownames(types[types$cell_type == "upperepidermal",])) Groups <- c(Groups, 'Upper')
  if (z %in% rownames(types[types$cell_type == "lowerepidermal",])) Groups <- c(Groups, 'Lower')
}

collist <- list(Groups=c("Upper"="#DD322E","Lower"="#40845C"))
ha <-HeatmapAnnotation(df = data.frame(Type = Groups),col = collist)
f2 <- colorRamp2(c(-2,0,2), c("#265297","#EEEEEE","#DD322E"))

df_p <- log2(data_fc)
df_p$gene <- rownames(df_p)
data_ed <- merge(degs,df_p,by = "gene")
data_ed$type <- factor(data_ed$type,levels = c("Upper","Lower"))
rownames(data_ed) <- data_ed$gene
mark_gene <- c("AT5G10310","AT2G17840")
gene_pos <- which(rownames(data_ed) %in% mark_gene)
labs <- rownames(data_ed[gene_pos,])
row_anno <-  rowAnnotation(mark_gene = anno_mark(at = gene_pos, 
                                                 labels = labs))

h <- Heatmap(as.matrix(data_ed[,c(-1,-2)]), name ="-log2 (FC)", column_title = " ",
             row_title = " ",
             row_names_side = "left",row_dend_side = "left",col = f2,
             cluster_columns = F,cluster_rows = F,split = data_ed$type,
             show_column_names = F,show_row_names = F,
             row_names_gp = gpar(fontsize = 10,font=3),right_annotation = row_anno,
             top_annotation = ha)
h 

