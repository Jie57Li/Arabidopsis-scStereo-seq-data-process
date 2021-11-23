###### R Script for Figure 2
# 1. library
#-------------------------------------------------------------------------------
# library(reshape2)
# library(ggplot2)
# library(factoextra)
library(Seurat)
library(dplyr)
library(patchwork)
library(openxlsx)
setwd("Projects/Fig2")

# 2. load data
#-------------------------------------------------------------------------------
obj <- readRDS("../Rds/Run_pre.rds")
meta <- obj@meta.data
meta$celltype_ed <- "NA"
meta[meta$cell_type %in% c("lowerepidermal","upperepidermal"),"celltype_ed"] <- "Epidermal"
meta[meta$cell_type %in% c("palisad","sponge"),"celltype_ed"] <- "Mesophyll"
meta[meta$cell_type %in% c("vascular"),"celltype_ed"] <- "Vascular"
meta[meta$cell_type %in% c("guardcells"),"celltype_ed"] <- "Guardcells"
meta$pca_type <- paste(meta$orig.ident,meta$celltype_ed,sep = "_")
obj@meta.data <- meta

# 3. Annotation
#-------------------------------------------------------------------------------
obj@active.ident <- factor(obj@meta.data$pca_type)
names(obj@active.ident) <- rownames(obj@meta.data)
ave <- AverageExpression(obj,slot = "counts")
data <- as.data.frame(ave$RNA)
type <- meta[,c("orig.ident","celltype_ed","pca_type")]
colnames(type) <- c("leaf","cell_type","sample")
type <- type[!duplicated(type$sample),]
#write.table(types,"heatmap-type.txt",row.names = T,sep = "\t",quote = F)
obj@active.ident <- factor(obj@meta.data$celltype_ed)
names(obj@active.ident) <- rownames(obj@meta.data)

library(ggrepel)
library(ggplot2)
celltype <- c("Vascular","Epidermal","Mesophyll","Guardcells")

# 4. DEGs and PCA
#-------------------------------------------------------------------------------
degs <- list()
n <-0
for (i in 1:3) {
  for (j in (i+1):4) {
    n <- n + 1
    var_genes <- FindMarkers(obj,ident.1 = celltype[i],ident.2 = celltype[j],logfc.threshold = 0.2,min.pct = 0.25,only.pos = F)
    var_genes <- var_genes[var_genes$p_val < 0.05, ]
    names <- paste(celltype[i],celltype[j],sep = "_")
    degs[[names]] <- var_genes
    reslut <- paste(n,names,sep = ".")
    if ( !dir.exists(reslut)) {
      print("0. Create specified directory...")
      dir.create(reslut)
    }
    
    # data
    data_t <- t(log10(data[rownames(var_genes),type[type$cell_type %in% celltype[c(i,j)],"sample"]]+1))
    
    # pca
    pca <- prcomp(data_t, scale=T,center = T)
    pca.table <- as.data.frame(pca$x)
    
    ### percent of PCs
    percent <- round((summary(pca)$importance[2,]*100)[1:4],2)
    print(percent)
    write.table(percent,paste(reslut,"/percent.txt",sep = ""),row.names = T,quote = F,sep = "\t")
    Plot
    pca.table$sample <- rownames(pca.table)
    dp <- merge(type,pca.table,by="sample")

    p <- ggplot(dp, aes(PC1, PC2, color = cell_type ,shape=cell_type)) +
      geom_point(size=4) + theme_bw() +
      labs(x = "PC1", y = "PC2") +
      theme(panel.background = element_rect(fill='transparent'),
            panel.grid=element_line(color='grey'),panel.border=element_rect(fill='transparent',color='black'),
            axis.title=element_text(size=15)) +
      theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank()) +
      theme(axis.text.x = element_text(size = 10,  color = "black", face = "bold", vjust = 0.5, hjust = 0.5)) +
      theme(axis.text.y = element_text(size = 10,  color = "black", face = "bold", vjust = 0.5, hjust = 0.5)) +
      geom_text_repel(aes(PC1, PC2, label = leaf))
    p
    # pdf(paste(reslut,"/01-PCA-Plot.pdf",sep = ""),width = 6,height = 5)
    # print(p)
    # dev.off()
    rownames(dp) <- dp$sample
    pca_res <- dp[,c("PC1","PC2","PC3")]
    pca_type <- as.data.frame(dp[,"cell_type"])
    rownames(pca_type) <- rownames(dp)
    write.table(pca_res,paste(reslut,"/PCA.txt",sep = ""),row.names = T,quote = F,sep = "\t")
    write.table(pca_type,paste(reslut,"/type.txt",sep = ""),row.names = T,quote = F,sep = "\t")
  }
}

# 5. the data of the heatmap
#-------------------------------------------------------------------------------
degs[[1]]$genes <- rownames(degs[[1]])
degs[[2]]$genes <- rownames(degs[[2]])
degs[[3]]$genes <- rownames(degs[[3]])
degs[[4]]$genes <- rownames(degs[[4]])
degs[[5]]$genes <- rownames(degs[[5]])
degs[[6]]$genes <- rownames(degs[[6]])
write.xlsx(degs,"DEGs-heatmap.xlsx")
saveRDS(obj,"../Rds/Run-PCA-871cells.rds")

# 6. find all markers for each cell type
#-------------------------------------------------------------------------------
markers <- FindAllMarkers(obj)
write.xlsx(markers,"DEGs-All-celltype.txt")

# 7. Dotplot 
#-------------------------------------------------------------------------------
d_g <- c("AT1G01170","AT5G59690","AT5G67600","AT4G04330","AT4G35490",
         "AT5G42270","AT4G16146","AT3G62530","AT2G04030","AT2G45860",
         "AT1G67865","AT4G00430","AT3G12260","AT1G19910","AT2G40880",
         "ATCG00020","AT1G20693","AT1G68660","AT4G28740","AT5G56030")
obj@active.ident <- factor(obj@active.ident,levels = c("Guardcells","Mesophyll","Epidermal","Vascular"))
d <- DotPlot(object = obj, features = d_g,cols = c("#C0C0C0", "#E41C12"),scale = T) + RotatedAxis() 
d
VlnPlot(obj, features = d_g,ncol = 5,pt.size = 0)
#--------------------------------End--------------------------------------------
