###### Expression pattern analysis of all cells and Enrichment

# 1. Library
#-------------------------------------------------------------------------------
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("Mfuzz")                   
# browseVignettes("Mfuzz")                          
library(Mfuzz)
library(Seurat)
library(DOSE)
library(org.At.tair.db)
library(topGO)
library(clusterProfiler)
library(pathview)
library(ggplot2)
setwd("Projects/Fig4")

# 2. Load data
#-------------------------------------------------------------------------------
s_merge <- readRDS("../RDS/Run.rds")
# s_merge <- s_merge[,grep(pattern="epi", rownames(s_merge@meta.data))]

all <- FindAllMarkers(s_merge, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
all <- unique(all$gene)
AverageExp <- as.data.frame(AverageExpression(s_merge)$RNA)
AverageExp <- AverageExp[all,]

# 3. mfuzz for pattern
#-------------------------------------------------------------------------------
data <- AverageExp
gene_tpm <- data.matrix(data)
eset <- new("ExpressionSet",exprs = gene_tpm)
gene.r <- filter.NA(eset, thres = 0.25)
## mean fill NA values
gene.f <- fill.NA(gene.r, mode="mean")
## knn/wknn方法表现更好，但是计算起来比较复杂
# gene.f <- fill.NA(gene.r, mode="knn")
gene.f <- fill.NA(gene.r, mode="wknn")
## 过滤标准差为0的基因
tmp <- filter.std(gene.f, min.std = 0)
## Standardise
gene.s <- standardise(tmp)
## Cluster numbers
c <- 3
## Estimate best m
m <- mestimate(gene.s)
## Cluster
cl <- mfuzz(gene.s, c = c, m = m)
## Check genes numbers in each cluster
cl$size
## Export gene cluster
write.table(cl$cluster, "All_cells/All_cell_clusters.txt", 
            quote = F, row.names = T, col.names = F, sep = "\t")

## Draw line chart
pdf("All_cells/All_cell_patterns.pdf", width = 6, height = 7)
#p <- 
  mfuzz.plot(gene.s, cl, mfrow=c(1,3), new.window = F, 
                time.labels = c("0", "1", "2", "3"), min.mem = 0.41)
dev.off()

# 4. Enrichment
#-------------------------------------------------------------------------------
genes <- as.data.frame(cl$cluster)
colnames(genes) <- "cluster"
genes$id <- rownames(genes)
file <- c("Cluster1","Cluster2","Cluster3")
#cluster <- c("1","2","3","4")
for (i in 1:length(file)) {
  dir.create(paste("All_cells/Enrichment/",i,sep = ""))
  id2 <- as.vector(genes[which(genes$cluster == i),2])
  # Transform gene names
    genelits = bitr(id2, 
                    fromType="TAIR",
                    toType="ENTREZID", 
                    OrgDb="org.At.tair.db")
    id <- as.vector(genelits[,2])
    
  ## 3. GO
  ego.all <- enrichGO(
    gene          = id,
    keyType = "ENTREZID",
    OrgDb         = org.At.tair.db,
    ont           = "BP",    # options : "CC" "BP" "MF" "ALL"
    pAdjustMethod = "BH",
    pvalueCutoff  = 1,
    qvalueCutoff  = 1,
    readable      = FALSE) # whether mapping gene ID to gene Name

  barplot(ego.all, showCategory=20,title="EnrichmentGO_ALL")
  # dotplot(ego.all,showCategory=20,title="EnrichmentGO_ALL_dot")
  write.table(ego.all,paste("All_cells/Enrichment/",i,"/Cluster",i,"_GO.results.txt",sep = ""),row.names = FALSE,sep = "\t")

  ## 4. KEGG
  ekk.all <- enrichKEGG(gene         = id2,
                        organism     = "ath",
                        pvalueCutoff = 1)
  barplot(ekk.all, showCategory=20,title="EnrichmentGO_ALL")
  write.table(ekk.all, paste("All_cells/Enrichment/",i,"/Cluster",i,"_KEGG.results.txt",sep = ""),row.names = FALSE,sep="\t")
#   
}

