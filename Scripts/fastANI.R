library(pheatmap)
library(reshape2)
library(readxl)
library(dplyr)

data1 <- read.table("D:/fastANI")
annota_table <- read.csv("C:/Users/67477/OneDrive/ARGs_genes/new/vcontact/genome_by_genome_overview.csv")[,0:2] |> as.data.frame()
data2 <- dcast(data1,V1~V2,value.var = "V3",fill = 0) |> as.data.frame()
data3 <- merge(data2,annota_table,by.x = "V1",by.y = "Genome",all = F)
annota <- data3 |> select(1,ncol(data3))
rownames(data2) <- data2[,1]
data2 <- data2[,-1]
rownames(annota) <- annota[,1]
annota <- annota |> select(2)
pheatmap(data2,show_rownames=F,show_colnames=F,color = colorRampPalette(c("navy", "white", "firebrick3"))(50),annotation_row = annota
         ,cluster_rows = F,cluster_cols = F)
