###### Trajectory analysis of vascular cells

# 1. Install and Library Packages 
#-------------------------------------------------------------------------------
# install.packages("BiocManager")
# BiocManager::install("monocle")
library(monocle)
library(RColorBrewer)
library(ggplot2)
library(Seurat)
library(ggpubr)
setwd("Projects/Fig5")

# 2. Find marker gene 
#-------------------------------------------------------------------------------
rds1 <- readRDS("../RDS/Run.rds")
## Vascular marker
sub <- rds1[,grep(pattern = "vascu", rownames(rds1@meta.data))]
vas <- FindMarkers(sub, ident.1 = 0, ident.2 = 3, min.pct = 0.14)
vas <- rownames(vas[which(vas$p_val<0.05),])
data <- data.frame(sub@assays$RNA@counts)
names(data) <- gsub("\\.", "-", names(data))
type <- data.frame(sub@meta.data)
type <- type[c("cell_type", "reassign_7cut")]
type$reassign_7cut <-paste("Position", type$reassign_7cut, sep="")

# 3. Data Sets
#-------------------------------------------------------------------------------
expr_matrix <- data[,rownames(type)]   # expression matrix
gene_annotation <- data.frame(gene_id=rownames(expr_matrix),
                              gene_short_name=rownames(expr_matrix))
rownames(gene_annotation) <- gene_annotation[,1]   # gene annotation
sample_sheet <- type   # corresponding cell annotation information

# 4. Preprocessing
#-------------------------------------------------------------------------------
pos <- which(rownames(gene_annotation) %in% rownames(expr_matrix))
gene_annotation <- gene_annotation[pos,]
expr_matrix <- expr_matrix[rownames(gene_annotation),rownames(sample_sheet)]
pd <- new("AnnotatedDataFrame", data = sample_sheet)   # cell annotation object
fd <- new("AnnotatedDataFrame", data = gene_annotation)   # gene annotation object

## Choose appropriate distribution according to expression data type
## UMI (default input): for negative binomal distribution data
cd <- newCellDataSet(as(as.matrix(expr_matrix), "sparseMatrix"), phenoData = pd,
                     featureData = fd, expressionFamily = negbinomial.size(),
                     lowerDetectionLimit = 0.5)
cd <- estimateSizeFactors(cd)   # help eliminate differences in mRNA capture between cells
cd <- estimateDispersions(cd)   # for subsequent differential expression analysis
cd <- detectGenes(cd, min_expr = 0.5)   # counts the number of genes expressed in each cell
expressed_genes <- row.names(subset(fData(cd), num_cells_expressed > nrow(sample_sheet) * 0.01))

# 5. Construction of Single Cell Pseudotime 
#-------------------------------------------------------------------------------
ordering_genes <- vas
ordering_genes <- intersect(ordering_genes, expressed_genes)
cd <- setOrderingFilter(cd, ordering_genes)
## plot_ordering_genes(cd)

## Order Cells by Progress
cd <- reduceDimension(cd, max_components = 2, method = 'DDRTree')
cd <- orderCells(cd, reverse = F)   # reverse = F by default

## generate color palette
getPalette <- colorRampPalette(brewer.pal(6, "Set1"))

# 6. Visualize in reduced dimensional Space 
#-------------------------------------------------------------------------------
# Visualize Cell Trajectory
color1 <- c("#CC503B","#0072B5","#D68518","#55A368","#7976B2","#779AAE","#F9DB8D","#F54E91")

# 7. Visualize in reduced dimensional Space 
#-------------------------------------------------------------------------------
# Visualize Cell Trajectory
pc1 <- plot(plot_cell_trajectory(cd, show_cell_names = F, color_by = "cell_type", cell_size = 2.5)) +
  scale_color_manual(values = "#B9C8A0")
pc2 <- plot(plot_cell_trajectory(cd, show_cell_names = F, color_by = "reassign_7cut", cell_size = 2.5))
# scale_color_manual(values = color1))
pc3 <- plot(plot_cell_trajectory(cd, show_cell_names = F,
                                 color_by = "Pseudotime", cell_size = 2.5) +
              scale_color_viridis_c())
# pc4 <- plot(plot_cell_trajectory(cd, show_cell_names = F, color_by = "State",cell_size = 3))

p <- ggarrange(pc1, pc2, pc3, nrow = 1)


# 8. Do heatmap 
#-------------------------------------------------------------------------------
# Whole
diff_test_res <- differentialGeneTest(cd[ordering_genes,],
                                      fullModelFormulaStr = "~sm.ns(Pseudotime)")
sig_gene_names <- row.names(subset(diff_test_res, qval < 0.4, pval < 0.05))

tst <- plot_pseudotime_heatmap(cd[sig_gene_names,],
                               num_clusters = 3,
                               cores = 1,
                               show_rownames = F)
plot_pseudotime_heatmap(cd[sig_gene_names,],
                        num_clusters = 3,
                        cores = 1,
                        show_rownames = F)

# 9. Plot canonical genes
#-------------------------------------------------------------------------------
library(ggpubr)
## Main figures
## Start
cds_subset <- cd[c("AT1G78370"),]
pg1 <- plot_genes_in_pseudotime(cds_subset, color_by = "reassign_7cut")
## Middle
cds_subset <- cd[c("AT5G44580"),]
pg2 <- plot_genes_in_pseudotime(cds_subset, color_by = "reassign_7cut")
## End
cds_subset <- cd[c("AT1G09932"),]
pg3 <- plot_genes_in_pseudotime(cds_subset, color_by = "reassign_7cut")
p1 <- ggarrange(pg1, pg2, pg3, nrow = 1)
## Supplementary figures
cds_subset <- cd[c("AT1G26630"),]
pg1 <- plot_genes_in_pseudotime(cds_subset, color_by = "reassign_7cut")
## Middle
cds_subset <- cd[c("AT5G42530"),]
pg2 <- plot_genes_in_pseudotime(cds_subset, color_by = "reassign_7cut")
## End
cds_subset <- cd[c("AT5G40690"),]
pg3 <- plot_genes_in_pseudotime(cds_subset, color_by = "reassign_7cut")
cds_subset <- cd[c("AT4G33980"),]
pg4 <- plot_genes_in_pseudotime(cds_subset, color_by = "reassign_7cut")
p2 <- ggarrange(pg1, pg2, pg3, pg4, nrow = 1)
