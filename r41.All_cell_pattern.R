##### R version 4.0.4 (2021-02-15)
##### Jie Li
##### 2021/03/12
##### export the expressed genes

# 1. library
# library(reshape2)
# library(ggplot2)
# library(factoextra)
library(Seurat)
library(dplyr)
library(patchwork)
library(openxlsx)
library(ggrepel)
library(ggplot2)

setwd("F:/BGI/Projects/2021/2.植物空间组/2.Analyses/02.AT-leaves/Fig4/")

# 2. load data
obj <- readRDS("../a0.RDS/reassigned_merged_umap_res0.4_0311v1.rds")
meta <- read.table("Pattern-markers/00.Area_data.txt",header = T)


# 2. ave (leaf + cell type)
obj@active.ident <- factor(meta$area_new)
names(obj@active.ident) <- rownames(meta)

# 4. Dotplot
#sub_obj <- subset(obj,celltype_ed = "Mesophyll")
obj@active.ident <- factor(obj@active.ident,levels = c("Area0","Area1","Area2","Area3"))
d_g <- c("AT1G14150", "AT1G19150","AT2G20890","AT2G39470", "AT5G66190")
d <- DotPlot(object = obj, features = d_g,cols = c("#C0C0C0", "#E41C12"),scale = T) #+ RotatedAxis() 
d

pdf("../../../4.Figs/Figs_0423/Fig4/Fig4d.pdf",width = 8.3,height = 3)
d
dev.off()
