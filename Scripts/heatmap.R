library(readxl)
library(dplyr)
library(pheatmap)
library(ggsci)

df <- read_xlsx("C:/Users/67477/Desktop/abundance_1/Abu_Transcriptome.xlsx",sheet = 3) |> as.data.frame()

rownames(df) <- df[,1]
df <- df[,-1]
# 创建一个因子并设置其水平

row_order <- c("Bioreactor", "Solid", 
               "Wastewater", "Other engineered", "Freshwater", "Marine", 
               "Thermal", "Soil", "Other enviroment", "Digestive", 
               "Plants", "Microbial", "Other Host-associated", "Unclassfied")
habitat_factor <- factor(df$habitat, levels = row_order)
# 对数据框进行排序
df_sorted <- df[order(habitat_factor), ]

annotation_row <- select(df_sorted,c(1:2))

df2 <- df_sorted[,-c(1:2)] 

annotation_col <- data.frame(
  Transcriptome = factor(ifelse(grepl("Waste", colnames(df2)), "Waste",
                        ifelse(grepl("TARAocean", colnames(df2)), "TARAocean",
                               ifelse(grepl("Freshwater", colnames(df2)), "Freshwater",
                                      ifelse(grepl("Soil", colnames(df2)), "Soil", "HumanGut")))))
)

rownames(annotation_col) <- colnames(df2)

ann_colors = list(
  Status = c(MIX = "red", Blank = "#FFFFFF", Prophage = "green"),
  habitat = c("Bioreactor" = "#E8D3C0",
              "Solid" = "#D89C7A",
              "Wastewater" = "#D6C38B",
              "Other engineered" = "#ECCED0",
              "Freshwater" = "#849B91",
              "Marine" = "#C2CEDC",
              "Thermal" = "#979771",
              "Soil" = "#686789",
              "Other enviroment" = "#B57C82",
              "Digestive" = "#B77F70",
              "Plants" = "#A79A89",
              "Microbial" = "#987465",
              "Other Host-associated" = "#D3D2D0",
              "Unclassfied" = "#808080")
)




pheatmap(df2, cluster_rows=F, cluster_cols = F, annotation_row = annotation_row,
         show_rownames = F, show_colnames = F,annotation_col = annotation_col,
         annotation_colors = ann_colors,
         color = c("#eeeeee","#808080","#00b8a9","#ffde7d","#f6416c","#9896f1","#3490de"))


