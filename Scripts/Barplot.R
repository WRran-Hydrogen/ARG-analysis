library(ggplot2)
library(dplyr)
library(readxl)
library(reshape2)
library(ggsci)

df <- read_xlsx("C:/Users/67477/Desktop/Abundance.xlsx",sheet = 3) |> as.data.frame()

# 使用mutate函数创建新列，这些新列是包含特定字符串的列的汇总
df_2 <- df %>%
  rowwise() %>%
  mutate(Freshwater = sum(c_across(starts_with("Fresh"))),
         Soil = sum(c_across(starts_with("Soil"))),
         HumanGut = sum(c_across(starts_with("HumanGut"))),
         Waste = sum(c_across(starts_with("Waste"))),
         TARAocean = sum(c_across(starts_with("TARAocean")))) %>%
  ungroup() %>%
  select(habitat, Status, Freshwater, Soil, HumanGut, Waste, TARAocean)


df_long <- melt(df_2, id.vars = c("habitat", "Status"), variable.name = "Environment", value.name = "Value")

# Environment&Status
df_sum <- df_long %>%
  group_by(Status, Environment) %>%
  summarise(Value = sum(Value))

df_percentage <- df_sum %>%
  group_by(Environment) %>%
  mutate(Percentage = Value / sum(Value) * 100)

# 使用ggplot函数创建堆叠柱状图
ggplot(df_percentage, aes(x = Environment, y = Percentage, fill = Status)) +
  geom_bar(stat = "identity", position = "fill",width = 0.3)+
  scale_fill_cosmic()+
  theme_classic()

# Habitat&Environment

df_sum2 <- df_long %>%
  group_by(habitat, Environment) %>%
  summarise(Value = sum(Value))

df_percentage2 <- df_sum2 %>%
  group_by(habitat) %>%
  mutate(Percentage = Value / sum(Value) * 100)

# 使用ggplot函数创建堆叠柱状图

colors <- c("Bioreactor" = "#E8D3C0",
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

ggplot(df_percentage2, aes(x = Environment, y = Percentage, 
                           fill = factor(habitat,levels = names(colors)))) +
  geom_bar(stat = "identity", position = "fill",width = 0.3)+
  scale_fill_manual(values = colors)+
  theme_classic()



