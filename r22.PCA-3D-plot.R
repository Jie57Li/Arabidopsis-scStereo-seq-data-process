##### 3D Plot of PCA (R Script for Figure 2 and Fig S2)

# 1. library
#-------------------------------------------------------------------------------
library("scatterplot3d")
setwd("Projects/Fig2")

# 2. load data
#-------------------------------------------------------------------------------
# Vascular and Epidermal
data <- read.table("1.Vascular_Epidermal/Vascular-Epidermal-PCA.txt",header = T)
iris <- data
iris$color <- "NA"
iris[iris$Type == "Epidermal","color"] <- "#9FD0F0"
iris[iris$Type == "Vascular","color"] <- "#844F9C"
scatterplot3d(iris$PC1,iris$PC2,iris$PC3,main="Basic 3D Scatter Plot",
              color = iris$color,pch = 16,angle = 20)

# Vascular and Mesophyll
pca <- read.table("2.Vascular_Mesophyll/PCA.txt",header = T)
pca$sample <- rownames(pca)
type <- read.table("2.Vascular_Mesophyll/type.txt",header = T,sep = "\t")
type$sample <- rownames(type)
data <- merge(pca,type,by="sample")
colnames(data) <- c("sample","PC1","PC2","PC3","type")

iris <- data
iris$color <- "NA"
iris[iris$type == "Mesophyll","color"] <- "#2F4827"
iris[iris$type == "Vascular","color"] <- "#844F9C"
scatterplot3d(iris$PC1,iris$PC2,iris$PC3,main="Basic 3D Scatter Plot",
              color = iris$color,pch = 16,angle = 8,cex.axis =1)

# Vascular and Guardcells
pca <- read.table("3.Vascular_Guardcells/PCA.txt",header = T)
pca$sample <- rownames(pca)
type <- read.table("3.Vascular_Guardcells/type.txt",header = T,sep = "\t")
type$sample <- rownames(type)
data <- merge(pca,type,by="sample")
colnames(data) <- c("sample","PC1","PC2","PC3","type")

iris <- data
iris$color <- "NA"
iris[iris$type == "Guardcells","color"] <- "#83A7A2"
iris[iris$type == "Vascular","color"] <- "#844F9C"
scatterplot3d(iris$PC1,iris$PC2,iris$PC3,main="Basic 3D Scatter Plot",
              color = iris$color,pch = 16,angle = -40,cex.axis =1)

# Epidermal and Mesophyll
pca <- read.table("4.Epidermal_Mesophyll/PCA.txt",header = T)
pca$sample <- rownames(pca)
type <- read.table("4.Epidermal_Mesophyll/type.txt",header = T,sep = "\t")
type$sample <- rownames(type)
data <- merge(pca,type,by="sample")
colnames(data) <- c("sample","PC1","PC2","PC3","type")

iris <- data
iris$color <- "NA"
iris[iris$type == "Epidermal","color"] <- "#9FD0F0"
iris[iris$type == "Mesophyll","color"] <- "#2F4827"
scatterplot3d(iris$PC1,iris$PC2,iris$PC3,main="Basic 3D Scatter Plot",
              color = iris$color,pch = 16,angle =-146,cex.axis =1)

# Epidermal and Guardcells
pca <- read.table("5.Epidermal_Guardcells/PCA.txt",header = T)
pca$sample <- rownames(pca)
type <- read.table("5.Epidermal_Guardcells/type.txt",header = T,sep = "\t")
type$sample <- rownames(type)
data <- merge(pca,type,by="sample")
colnames(data) <- c("sample","PC1","PC2","PC3","type")

iris <- data
iris$color <- "NA"
iris[iris$type == "Epidermal","color"] <- "#9FD0F0"
iris[iris$type == "Guardcells","color"] <- "#83A7A2"
scatterplot3d(iris$PC2,iris$PC1,iris$PC3,main="Basic 3D Scatter Plot",
              color = iris$color,pch = 16,angle =78,cex.axis =1)

# Mesophyll and Guardcells
pca <- read.table("6.Mesophyll_Guardcells/PCA.txt",header = T)
pca$sample <- rownames(pca)
type <- read.table("6.Mesophyll_Guardcells/type.txt",header = T,sep = "\t")
type$sample <- rownames(type)
data <- merge(pca,type,by="sample")
colnames(data) <- c("sample","PC1","PC2","PC3","type")

iris <- data
iris$color <- "NA"
iris[iris$type == "Mesophyll","color"] <- "#2F4827"
iris[iris$type == "Guardcells","color"] <- "#83A7A2"
scatterplot3d(iris$PC2,iris$PC1,iris$PC3,main="Basic 3D Scatter Plot",
              color = iris$color,pch = 16,angle =-145,cex.axis =1)
#-------------------------------------------------------------------------------
