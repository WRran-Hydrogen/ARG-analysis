library(gggenomes)
library(readxl)
library(ggsci)

df_genes <- read_xlsx("C:/Users/67477/OneDrive/ARGs_genes/Inphage_gggenome.xlsx",sheet=3)
df_seqs <- read_xlsx("C:/Users/67477/OneDrive/ARGs_genes/Inphage_gggenome.xlsx",sheet=4)
gggenomes(df_genes) |>
  add_subfeats(df_seqs,.transform = "none")+
  geom_seq()+
  geom_bin_label()+
  geom_gene(color="#000000",fill="#FFFFFF",shape = 5,size = 5)+
  geom_feat(aes(color=type),linewidth=7,position="identity")+
  scale_color_d3(palette = "category20")+
  geom_feat_tag(aes(label=type))
