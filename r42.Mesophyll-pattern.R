#------------------------------------------------#
#               Spatial AT Analysis              #
#     Editor: Ruiying Chen   TEL:18995624500     #
#              Edited on: 2021.3.7               #
#          Expression pattern analysis           #
#------------------------------------------------#

# Library
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("Mfuzz")  Mfuzz: "2.50.0"                 
# browseVignettes("Mfuzz")                          
library(Mfuzz)
library(Seurat)
library(DOSE)
library(org.At.tair.db)
library(topGO)
library(clusterProfiler)
library(pathview)
library(ggplot2)
library(openxlsx)

# Load data
# mfuzz for pattern
setwd("G:/BGI/Projects/2021/2.植物空间组/2.Analyses/02.AT-leaves/Fig4/")
data <- read.table("Mesophyll_DEG_average_0312.txt", header = TRUE, row.names = 1)[,c(4,3,2,1)]
outdir <- "Mesophyll/"
if (!dir.exists(outdir)){
  dir.create(outdir,recursive = T) # mkdir -p
}

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
c <- 2
## Estimate best m
m <- mestimate(gene.s)
## Cluster
cl <- mfuzz(gene.s, c = c, m = m)
## Check genes numbers in each cluster
cl$size
## Export gene ID
# write.table(cl$cluster, "results/pattern/Mild_mfuzz.txt", quote = F, row.names = T, col.names = F, sep = "\t")

## Draw line chart
pdf(paste(outdir,"Mesophyll-pattern.pdf",sep = ""), width = 5, height = 5)
mfuzz.plot(gene.s, cl, mfrow=c(2,2), new.window = FALSE, 
           time.labels = c("0", "1", "2", "3"))
dev.off()

# # 3. Enrichment
genes <- as.data.frame(cl$cluster)
colnames(genes) <- "cluster"
genes$id <- rownames(genes)
go <- list()
kegg <- list()

file <- c("Cluster1","Cluster2")
#cluster <- c("1","2","3","4")
for (i in 1:length(file)) {
  result <- paste(outdir,"/Cluster-",i,sep = "")
  if (!dir.exists(result)){
    dir.create(result,recursive = T) # mkdir -p
  }
  names <- paste("Clustre",i,sep = "")
  id2 <- as.vector(genes[which(genes$cluster == i),2])
  # Transform gene names
  genelits = bitr(id2, 
                  fromType="TAIR",
                  toType="ENTREZID", 
                  OrgDb="org.At.tair.db")
  write.table(genelits,paste(result,"/01.genelist.txt",sep = ""),quote = F,sep = "\t",row.names = F)
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
  write.table(ego.all,paste(result,"/02.Cluster",i,"_GO.results.txt",sep = ""),row.names = FALSE,sep = "\t")
  go[[names]] <- ego.all
  
  ## 4. KEGG
  ekk.all <- enrichKEGG(gene         = id2,
                        organism     = "ath",
                        pvalueCutoff = 1)
  barplot(ekk.all, showCategory=20,title="EnrichmentGO_ALL")
  write.table(ekk.all, paste(result,"/03.Cluster",i,"_KEGG.results.txt",sep = ""),row.names = FALSE,sep="\t")
  kegg[[names]] <- ekk.all
}

write.xlsx(go,paste(outdir,"/GO-result.xlsx",sep = ""))
write.xlsx(kegg,paste(outdir,"/KEGG-result.xlsx",sep = ""))
