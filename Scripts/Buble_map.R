library(ggplot2)
library(maps)
library(readxl)


data1 <- read_xlsx("c:/Users/67477/OneDrive/ARGs_genes/new/information.xlsx",sheet = 2) |> as.data.frame()
data1$Longitude <- as.numeric(data1$Longitude)
data1$Latitude <- as.numeric(data1$Latitude)
data1$Habitat <- factor(data1$Habitat,level = c("Bioreactor","Solid","Wastewater","Other engineered","Freshwater","Marine","Thermal","Soil","Other enviroment","Digestive","Plants","Human Respiratory","Microbial","Other Host-associated"))

# 加载世界地图数据
world_map <- map_data("world")

# 绘制世界地图
ggplot() +
  geom_polygon(data = world_map, aes(x = long, y = lat, group = group),fill = "lightgray") +
  coord_fixed(ratio = 1.3) +
  # 添加散点图
  geom_point(data = data1, aes(x = Longitude, y = Latitude, shape = pfam, color = Habitat, size = Hit)) +
  # 设置形状和颜色的映射关系
  scale_shape_manual(values = c("PF00905" = 15, "PF00144" = 16, "PF13354" = 17)) +
  scale_size(range = c(1,7)) +
  scale_color_manual(values = c("#E3DCA5","#E0E0A1","#F5F5DC","#D1D1B1","#9ACD32","#8FBC8F","#7CB68E","#6CA67C","#5D956B","#FFB6C1","#FFC0CB","#FFA7B3","#FF99A6","#FF8A98"))+
  labs(x = "longitude", y = "latitude", color = "habitat", shape = "pfam", size = "number") +
  # 自定义图例
  guides(color = guide_legend(title.position = "top"),
         shape = guide_legend(title.position = "top")) +
  # 自定义主题
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5),
        legend.background=element_blank(),
        panel.background=element_rect(fill='white'))

                                      