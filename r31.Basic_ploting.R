##### R Script for Figure 3 and Fig S3

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

# 2. load data
#-------------------------------------------------------------------------------
obj <- readRDS("../a0.RDS/Run-PCA-871cells.rds")
obj@active.ident <- factor(obj@active.ident,levels = c("Vascular","Epidermal","Mesophyll","Guardcells"))

# 2. Vlnplot (Fig. 3a and Fig. S3a)
#-------------------------------------------------------------------------------
marker_ed <- c("AT1G78370","AT2G45860","AT5G53940","AT5G56030","AT5G67600","AT2G47710","AT1G16890")
col <- c("#A1B198","#C82C30","#364829","#83A4A0")
VlnPlot(obj, features = marker_ed,ncol = 4,pt.size = 0,cols = col)

# 3. Enrichment
#-------------------------------------------------------------------------------
library(DOSE)
library(topGO)
library(clusterProfiler)
library(pathview)
library(ggplot2)
# BiocManager::install('org.At.tair.db')
# BiocManager::install('pathview')
library(org.At.tair.db)   
keytypes(org.At.tair.db)
library(openxlsx)

file <- list.files("./") # the file of markers for each cell type
go <- list()
## Run GO   
for (j in 1:4) {
  i <- file[j]
  name <- gsub(pattern = "\\.txt$",replacement = "",x = i)
  result <- name
  if ( !dir.exists(result)) {
    print("0. Create specified directory...")
    dir.create(result)
  }
  data <- read.table(i,header = T)
  gene <- as.vector(data$genes)
  genelits = bitr(gene, 
                  fromType="TAIR", #输入为SYMBOL格式
                  toType="ENTREZID",  # 转为ENTERZID格式
                  OrgDb="org.At.tair.db")
  write.table(genelits,paste(result,"/01.",name,"_genelist.txt",sep = ""),sep = "\t",row.names = FALSE,quote = FALSE)
  id <- as.vector(genelits[,2])
  
  ### GO
  ego.all <- enrichGO(
    gene          = id,
    keyType = "ENTREZID",
    OrgDb         = org.At.tair.db,
    ont           = "ALL",    # optional : "CC" "BP" "MF" "ALL"
    pAdjustMethod = "BH",
    pvalueCutoff  = 1,
    qvalueCutoff  = 1,
    readable      = FALSE)
  go[[name]] <- as.data.frame(ego.all)
  barplot(ego.all, showCategory=20,title="EnrichmentGO_ALL")
  # dotplot(ego.all,showCategory=20,title="EnrichmentGO_ALL_dot")
  write.table(ego.all,paste(result,"/02.",name,"_GO.reaults.txt",sep = ""),row.names = FALSE,sep = "\t")
}
write.xlsx(go,"AT_GO_result.xlsx")

# 4. heatmap of Enrichment (selected)
#-------------------------------------------------------------------------------
library(ggplot2)
library(cowplot)
library(reshape2)
library(factoextra)
library(ComplexHeatmap)
library(grid)
library(circlize)
library("GetoptLong")

f2 <- colorRamp2(c(0,3), c("#C0C0C0", "#E41C12"))
df_p <- read.table("GO.txt",header = T,sep = "\t")
l <- dim(df_p)[2]
df_p[,3:l] <- -log10(df_u[,3:l])
p <- Heatmap(as.matrix(df_p[,c(-1,-2)]), name ="-log10 (p)", 
             column_title = " ",row_title = " ",
             row_names_side = "left",row_dend_side = "left",col = f2,
             cluster_columns = F,cluster_rows = F,split = df_p$Tissue,
             show_column_names = T,show_row_names = T,
             row_names_gp = gpar(fontsize = 10,font=3),
             row_dend_width = unit(10, "mm"))
#-------------------------------------------------------------------------------
