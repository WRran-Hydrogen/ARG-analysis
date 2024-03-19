library(treemapify)
library(ggplot2)

mycolors <- c("Bioreactor" = "#E8D3C0",
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


myorder <- c("Bioreactor", "Solid", "Wastewater", "Other engineered", 
             "Freshwater", "Marine", "Thermal", "Soil", 
             "Other enviroment", "Digestive", "Plants", 
             "Microbial", "Other Host-associated", "Unclassfied")

df <- read.table("C:/Users/67477/Downloads/abu.txt",header = T, sep = "\t")
ggplot(df, aes(area = Hit, fill = Habitat2, subgroup = Habitat1,label = Habitat2)) +
  geom_treemap(colour = "white") +
  geom_treemap_subgroup_border(colour = "black") +
  geom_treemap_text(fontface = "italic", colour = "black", place = "topleft",min.size = 0) +
  theme(legend.position = "bottom") +
  scale_fill_manual(values = mycolors,breaks = myorder) +
  theme(legend.position = "right")+
  labs(fill = "Habitat")
