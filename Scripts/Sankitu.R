library(networkD3)
library(readxl)
library(dplyr)

data <- read_xlsx("C:/Users/67477/Desktop/aafile.xlsx",sheet=3)

links <- data |> group_by(habitat,predict,class,b_phylum) |> summarise(count = n()) |> ungroup()


# 将数据转换为节点和链接
nodes <- data.frame(name = unique(c(as.character(links$habitat),
                                    as.character(links$predict),
                                    as.character(links$class), 
                                    as.character(links$b_phylum))))

links1 <- data.frame(source = match(links$habitat, nodes$name) - 1,
                     target = match(links$class, nodes$name) - 1,
                     value = links$count,
                     group = links$habitat)
links2 <- data.frame(source = match(links$class, nodes$name) - 1,
                     target = match(links$b_phylum, nodes$name) - 1,
                     value = links$count,
                     group = links$class)
links3 <- data.frame(source = match(links$b_phylum, nodes$name) - 1,
                     target = match(links$predict, nodes$name) - 1,
                     value = links$count,
                     group = links$b_phylum)

links_all <- rbind(links1, links2, links3)

links_all <- links_all[complete.cases(links_all),]
# 定义颜色
colourScale <- 'd3.scaleOrdinal()
.domain([
  "Bioreactor",
  "Solid",
  "Wastewater",
  "Other engineered",
  "Freshwater",
  "Marine",
  "Thermal",
  "Soil",
  "Other enviroment",
  "Digestive",
  "Plants",
  "Microbial",
  "Other Host-associated",
  "Unclassfied"
])
.range([
  "#E8D3C0",
  "#D89C7A",
  "#D6C38B",
  "#ECCED0",
  "#849B91",
  "#C2CEDC",
  "#979771",
  "#686789",
  "#B57C82",
  "#B77F70",
  "#A79A89",
  "#987465",
  "#D3D2D0",
  "#808080"
]);'


# 绘制桑基图
sankeyNetwork(Links = links_all, Nodes = nodes, Source = "source",
              Target = "target", Value = "value", NodeID = "name",
              fontSize = 12, nodeWidth = 10,LinkGroup = 'group',colourScale=colourScale)

sankeyNetworkOutput(p, width = "100%", height = "500px")
