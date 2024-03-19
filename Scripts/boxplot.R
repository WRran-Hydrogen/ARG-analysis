library(ggplot2)
library(tidyr)
library(dplyr)

#读入表格
df <- read.table("C:/Users/67477/Downloads/box.txt",header = T, sep = "\t",na.strings = "NA")

#自定义变量参数
myorder <- c("Bioreactor", "Solid", "Wastewater", "Other engineered", 
             "Freshwater", "Marine", "Thermal", "Soil", 
             "Other enviroment", "Digestive", "Plants", 
             "Microbial", "Other Host-associated", "Unclassfied")

#数据透视表
df2 <- df %>%
  group_by(Habitat, PFAM) %>%
  summarise(Number = sum(Number)) %>%
  spread(key = PFAM, value = Number) %>%
  replace_na(list(PF00144 = 0, PF00905 = 0, PF13354 = 0))

df2 <- df |> select(2,4,5) |> spread(key = PFAM, value = Number) |> replace_na(list(PF00144 = 0, PF00905 = 0, PF13354 = 0))
#PCA分析
pca_result <- prcomp(df2[, c("PF00144", "PF00905", "PF13354")], scale. = TRUE)
#提取PCA坐标轴

pca_df <- data.frame(id = df$id, Habitat = df$Habitat, pca_result$x[, 1:2])
#提取PCA两坐标的贡献度
pca_summary <- summary(pca_result)
pc1_contribution <- pca_summary$importance[2, 1]
pc2_contribution <- pca_summary$importance[2, 2]
# 使用ggplot2绘制散点图
pca_df$Habitat <- factor(pca_df$Habitat, levels = myorder)
ggplot(pca_df, aes(x = PC1, y = PC2, color = Habitat)) +
  geom_point()+
  labs(x =  paste('PCA1:', pc1_contribution, '%'), y = paste('PCA2:', pc2_contribution, '%'), color = '')

biplot(pca_result)

ggplot(df,aes(x=Habitat,y=Hit,color=PFAM))+
  geom_point()+
  scale_y_continuous(breaks = c(0.0000001, 0.000001, 2))+
  theme_bw()

ggplot(df, aes(x = Habitat, y = Hit_number)) +
  geom_violin() +

  theme_bw()
