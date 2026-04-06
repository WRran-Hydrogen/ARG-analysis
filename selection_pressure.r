# HYPHY ABSREL

### New analysis for the selection pressure on the branches with significant positive selection
library(ggplot2)
library(dplyr)
library(tools)
library(readr)
library(aplot)

df <- read_tsv("c:\\Users\\WRran\\Desktop\\pf13354_info.tsv")
df2 <- read_tsv("c:\\Users\\WRran\\Desktop\\pf00905_info.tsv")
df3 <- read_tsv("c:\\Users\\WRran\\Desktop\\pf00144_info.tsv")
# 3. 数据清洗与节点分类打标签
df_13354 <- df %>%
  filter(LRT >= 0) %>%
  mutate(
    # 计算 -log10(p)
    neg_log_p = -log10(p_used + 1e-15),
    # 标记显著性 (控制颜色)
    Significance = case_when(
      p_used < 0.05 ~ internal_type,
      TRUE ~ "Not Significant"
    ),
  )

# 4. 绘制形状映射火山图
pf13354 <- ggplot(df_13354, aes(x = prop_omega_gt1, y = neg_log_p, 
                          color = Significance, shape = group)) +
  geom_point(size = 4, alpha = 0.8) +
  scale_color_manual(values = c("Not Significant" = "#d3d3d3", 
                                "virus" = "#3498db", 
                                "bacteria" = "#e74c3c", 
                                "Nested" = "#4daf4a")) +
  # 自定义形状 (16=实心圆/病毒, 15=实心方块/细菌, 17=实心三角/内部节点)
  scale_shape_manual(values = c("virus" = 16, 
                                "bacteria" = 15, 
                                "internal_node" = 17)) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", size = 0.6) +
  labs(
    title = "β-lactamase2 (PF13354)",
    x = expression(paste("Proportion of Sites with ", omega > 1)),
    y = expression(paste("-log"[10], " (p-value)")),
    color = "Selection Status",
    shape = "Lineage Type"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    plot.subtitle = element_text(
      hjust = 0.5, color = "#555555", margin = margin(b = 15)
    ),
    axis.title.x = element_text(margin = margin(t = 12), face = "bold"),
    axis.title.y = element_text(margin = margin(r = 12), face = "bold"),
    legend.position = "none",
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8)
  ) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 15))

## pf00905
df_00905 <- df2 %>%
  filter(LRT >= 0) %>%
  mutate(
    # 计算 -log10(p)
    neg_log_p = -log10(p_used + 1e-15),
    # 标记显著性 (控制颜色)
    Significance = case_when(
      p_used < 0.05 ~ internal_type,
      TRUE ~ "Not Significant"
    ),
  )

# 4. 绘制形状映射火山图
pf00905 <- ggplot(df_00905, aes(x = prop_omega_gt1, y = neg_log_p, 
                          color = Significance, shape = group)) +
  geom_point(size = 4, alpha = 0.8) +
  scale_color_manual(values = c("Not Significant" = "#d3d3d3", 
                                "virus" = "#3498db", 
                                "bacteria" = "#e74c3c", 
                                "Nested" = "#4daf4a")) +
  # 自定义形状 (16=实心圆/病毒, 15=实心方块/细菌, 17=实心三角/内部节点)
  scale_shape_manual(values = c("virus" = 16, 
                                "bacteria" = 15, 
                                "internal_node" = 17)) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", size = 0.6) +
  labs(
    title = "Transpeptidase (PF00905)",
    x = expression(paste("Proportion of Sites with ", omega > 1)),
    y = expression(paste("-log"[10], " (p-value)")),
    color = "Selection Status",
    shape = "Lineage Type"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    plot.subtitle = element_text(
      hjust = 0.5, color = "#555555", margin = margin(b = 15)
    ),
    axis.title.x = element_text(margin = margin(t = 12), face = "bold"),
    axis.title.y = element_text(margin = margin(r = 12), face = "bold"),
#    legend.position = "none",
    panel.grid.minor = element_blank()
  )

df_00144 <- df3 %>%
  filter(LRT >= 0) %>%
  mutate(
    # 计算 -log10(p)
    neg_log_p = -log10(p_used + 1e-15),
    # 标记显著性 (控制颜色)
    Significance = case_when(
      p_used < 0.05 ~ internal_type,
      TRUE ~ "Not Significant"
    ),
  )

# 4. 绘制形状映射火山图
pf00144 <- ggplot(df_00144, aes(x = prop_omega_gt1, y = neg_log_p, 
                          color = Significance, shape = group)) +
  geom_point(size = 4, alpha = 0.8) +
  scale_color_manual(values = c("Not Significant" = "#d3d3d3", 
                                "virus" = "#3498db", 
                                "bacteria" = "#e74c3c", 
                                "Nested" = "#4daf4a")) +
  # 自定义形状 (16=实心圆/病毒, 15=实心方块/细菌, 17=实心三角/内部节点)
  scale_shape_manual(values = c("virus" = 16, 
                                "bacteria" = 15, 
                                "internal_node" = 17)) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", size = 0.6) +
  labs(
    title = "β-lactamase (PF00144)",
    x = expression(paste("Proportion of Sites with ", omega > 1)),
    y = expression(paste("-log"[10], " (p-value)")),
    color = "Selection Status",
    shape = "Lineage Type"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    plot.subtitle = element_text(
      hjust = 0.5, color = "#555555", margin = margin(b = 15)
    ),
    axis.title.x = element_text(margin = margin(t = 12), face = "bold"),
    axis.title.y = element_text(margin = margin(r = 12), face = "bold"),
    legend.position = "none",
    panel.grid.minor = element_blank()
  )


pf13354 | pf00144 | pf00905


## HYPHY MG94 analysis
source("E:\\Onedrive-Yao\\OneDrive\\ARGs_genes\\new\\Scripts\\mg94_pairwise_plot_functions.r")

res_pf13354 <- make_mg94_family_plot(
  family_name = "PF13354",
  global_json = "c:\\Users\\WRran\\Desktop\\pf13354_MG94_global.json",
  full_json   = "c:\\Users\\WRran\\Desktop\\pf13354_MG94_part.json",
  bg_fg_json  = "c:\\Users\\WRran\\Desktop\\pf13354_MG94_BGvsFG.json",
  bg_nt_json  = "c:\\Users\\WRran\\Desktop\\pf13354_MG94_BGvsNT.json",
  fg_nt_json  = "c:\\Users\\WRran\\Desktop\\pf13354_MG94_FGvsNT.json"
)

res_pf13354$plot

res_pf00144 <- make_mg94_family_plot(
  family_name = "PF00144",
  global_json = "c:\\Users\\WRran\\Desktop\\pf00144_MG94_global.json",
  full_json   = "c:\\Users\\WRran\\Desktop\\pf00144_MG94_part.json",
  bg_fg_json  = "c:\\Users\\WRran\\Desktop\\pf00144_MG94_BGvsFG.json",
  bg_nt_json  = "c:\\Users\\WRran\\Desktop\\pf00144_MG94_BGvsNT.json",
  fg_nt_json  = "c:\\Users\\WRran\\Desktop\\pf00144_MG94_FGvsNT.json"
)

res_pf00144$plot

res_pf00905 <- make_mg94_family_plot(
  family_name = "PF00905",
  global_json = "c:\\Users\\WRran\\Desktop\\pf00905_MG94_global.json",
  full_json   = "c:\\Users\\WRran\\Desktop\\pf00905_MG94_part.json",
  bg_fg_json  = "c:\\Users\\WRran\\Desktop\\pf00905_MG94_BGvsFG.json",
  bg_nt_json  = "c:\\Users\\WRran\\Desktop\\pf00905_MG94_BGvsNT.json",
  fg_nt_json  = "c:\\Users\\WRran\\Desktop\\pf00905_MG94_FGvsNT.json"
)

res_pf00905$plot

res_pf13354$plot | res_pf00144$plot | res_pf00905$plot
