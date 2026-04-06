# ==============================================================================
# Script: Motif Usage Preference Analysis (Viral Habitats vs. Bacterial Control)
# ==============================================================================
library(tidyverse)
library(ggplot2)
# ------------------------------------------------------------------------------
# 1. Load Data
# ------------------------------------------------------------------------------
# Please ensure the file paths are correct for your working directory
edges <- read_csv("c:\\Users\\WRran\\Desktop\\motif_edge.csv")
nodes <- read_csv("c:\\Users\\WRran\\Desktop\\motif_node.csv")

# Define analytical categories
bacterial_habitat <- "BACTERIAL"
motif_categories <- c("M1", "M2", "M3")

# Extract Motif metadata (link each Motif ID to its Category M1/M2/M3)
motif_info <- nodes %>%
  filter(taxonomy %in% motif_categories) %>%
  select(Id, Motif_Category = taxonomy)

# Define Viral habitats dynamically (excluding Motifs and BACTERIAL)
viral_habitats <- nodes %>%
  filter(!taxonomy %in% c(motif_categories, bacterial_habitat)) %>%
  drop_na(taxonomy) %>%
  pull(taxonomy) %>%
  unique()

# ------------------------------------------------------------------------------
# 2. Viral Analysis (Background = Other Viral Habitats)
# ------------------------------------------------------------------------------
viral_habitat_totals <- nodes %>%
  filter(taxonomy %in% viral_habitats) %>%
  group_by(taxonomy) %>%
  summarise(Total_Sequences = n(), .groups = "drop")

viral_edge_stats <- edges %>%
  left_join(nodes %>% select(Id, taxonomy), by = c("Target" = "Id")) %>%
  rename(Sequence_Taxonomy = taxonomy) %>%
  filter(Sequence_Taxonomy %in% viral_habitats) %>%
  group_by(Source, Sequence_Taxonomy) %>%
  summarise(Observed_Count = n(), .groups = "drop")

viral_enrichment_df <- viral_edge_stats %>%
  left_join(viral_habitat_totals, by = c("Sequence_Taxonomy" = "taxonomy")) %>%
  mutate(Percentage = (Observed_Count / Total_Sequences) * 100)

viral_motif_totals <- viral_enrichment_df %>%
  group_by(Source) %>% summarise(Motif_Total = sum(Observed_Count), .groups = "drop")

viral_habitat_motif_totals <- viral_enrichment_df %>%
  group_by(Sequence_Taxonomy) %>% summarise(Habitat_Motif_Total = sum(Observed_Count), .groups = "drop")

viral_grand_total <- sum(viral_enrichment_df$Observed_Count)

viral_stats <- viral_enrichment_df %>%
  left_join(viral_motif_totals, by = "Source") %>%
  left_join(viral_habitat_motif_totals, by = "Sequence_Taxonomy") %>%
  mutate(
    A = Observed_Count, B = Motif_Total - A,
    C = Habitat_Motif_Total - A, D = viral_grand_total - A - B - C
  ) %>%
  rowwise() %>%
  mutate(
    p_value = tryCatch(fisher.test(matrix(c(A, B, C, D), nrow = 2))$p.value, error = function(e) 1)
  ) %>%
  ungroup() %>%
  mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
  left_join(motif_info, by = c("Source" = "Id"))

# ------------------------------------------------------------------------------
# 3. Bacterial Analysis (Background = All Viruses combined)
# ------------------------------------------------------------------------------
bact_total_seqs <- nodes %>% filter(taxonomy == bacterial_habitat) %>% nrow()

bact_edge_stats <- edges %>%
  left_join(nodes %>% select(Id, taxonomy), by = c("Target" = "Id")) %>%
  filter(taxonomy == bacterial_habitat) %>%
  group_by(Source) %>%
  summarise(Observed_Count = n(), .groups = "drop") %>%
  mutate(Sequence_Taxonomy = bacterial_habitat)

viral_pool_stats <- edges %>%
  left_join(nodes %>% select(Id, taxonomy), by = c("Target" = "Id")) %>%
  filter(taxonomy %in% viral_habitats) %>%
  group_by(Source) %>%
  summarise(Viral_Pool_Count = n(), .groups = "drop")

bact_total_motifs <- sum(bact_edge_stats$Observed_Count)
viral_total_motifs <- sum(viral_pool_stats$Viral_Pool_Count)

bact_stats <- bact_edge_stats %>%
  left_join(viral_pool_stats, by = "Source") %>%
  mutate(Viral_Pool_Count = replace_na(Viral_Pool_Count, 0)) %>%
  mutate(
    Percentage = (Observed_Count / bact_total_seqs) * 100,
    A = Observed_Count, B = Viral_Pool_Count,
    C = bact_total_motifs - A, D = viral_total_motifs - B
  ) %>%
  rowwise() %>%
  mutate(
    p_value = tryCatch(fisher.test(matrix(c(A, B, C, D), nrow = 2))$p.value, error = function(e) 1)
  ) %>%
  ungroup() %>%
  mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
  left_join(motif_info, by = c("Source" = "Id"))

# ------------------------------------------------------------------------------
# 4. Combine Results and Format for Integrated Plotting
# ------------------------------------------------------------------------------
col_needed <- c("Source", "Sequence_Taxonomy", "Percentage", "p_adj", "Motif_Category")

combined_stats <- bind_rows(
  viral_stats %>% select(all_of(col_needed)),
  bact_stats %>% select(all_of(col_needed))
) %>%
  drop_na(Motif_Category) %>%
  mutate(Sequence_Taxonomy = factor(Sequence_Taxonomy, levels = c(bacterial_habitat, viral_habitats)))

# Filter top motifs (e.g., reaching at least 2% in any group) to avoid overcrowded Y-axis.
# Adjust the '2.0' threshold up or down based on your dataset size.
top_motifs <- combined_stats %>%
  group_by(Source) %>%
  filter(max(Percentage) >= 2.0) %>% 
  pull(Source) %>%
  unique()

plot_data <- combined_stats %>% 
  filter(Source %in% top_motifs) %>%
  # Crucial Step: Create a new variable for Significance Coloring
  # If P.adj < 0.05, calculate -log10(P), otherwise set to NA (will be colored grey)
  mutate(LogP = ifelse(p_adj < 0.05, -log10(p_adj), NA))

# ------------------------------------------------------------------------------
# 5. Integrated Visualization (Usage % + Significance Color)
# ------------------------------------------------------------------------------
p_integrated <- ggplot(plot_data, aes(x = Sequence_Taxonomy, y = Source)) +
  
  # shape = 21 allows both a 'fill' (color inside) and a 'color' (border)
  geom_point(aes(size = Percentage, fill = LogP), shape = 21, color = "black", alpha = 0.85) +
  
  facet_grid(Motif_Category ~ ., scales = "free_y", space = "free_y") +
  
  # Color scale: NA values (non-significant) are set to grey
  # Significant values are colored from light yellow to deep red
  scale_fill_gradientn(
    colors = c("#ffffbf", "#fdae61", "#d7191c"), 
    na.value = "grey85",
    name = "Significance\n-log10(P.adj)\n[Grey = P ≥ 0.05]"
  ) +
  
  # Size scale based on usage percentage
  scale_size_continuous(range = c(2, 10), name = "Usage (%)") +
  
  theme_bw() +
  labs(
    title = "Motif Usage Preference and Significant Enrichment",
    subtitle = "Control: BACTERIAL vs Viruses | Viral Habitats: Specific Habitat vs Other Viruses",
    x = "Taxonomy / Habitat",
    y = "Motif Sequence"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 11, face = "bold", color = "black"),
    axis.text.y = element_text(size = 9, color = "black"),
    strip.background = element_rect(fill = "#404040"),
    strip.text = element_text(color = "white", size = 12, face = "bold"),
    panel.grid.major.x = element_line(color = "grey80", linetype = "dashed"),
    panel.grid.minor = element_blank(),
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 10)
  ) +
  
  # Dashed line strictly isolating the Bacterial control group
  geom_vline(xintercept = 1.5, linetype = "dashed", color = "black", size = 0.8)

print(p_integrated)