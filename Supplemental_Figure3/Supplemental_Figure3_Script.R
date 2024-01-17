#Supplemental Figure 3

require(toolboxH)
require(here)
require(patchwork)
require(ggplot2)

#S3
plotdata = fread("Data/h295_2_ergtabm_spearman.csv")

NogpaletteS3 <-  c("#A85C85","#542E42", "#00ABE7", "#0081AF", "#246A73", "#5CC1BC", "#368F8B", "#62C370")

plotdata[,.N, variable]
dput(plotdata$celltype %>% unique())
plotdata[,variable := factor(variable, levels = unique(variable))]
plotdata[,celltype:=factor(celltype, levels = c("Classical_Monocytes", "Non_Classical_Monocytes", "Immature_Neutrophils_type1", "Neutrophils", "NK_Cells", "CD4+_T_Cells", "CD8+_T_Cells", "B_Cells"))]
levels2 <- c("Classical_Monocytes", "Non_Classical_Monocytes", "Immature_Neutrophils_type1", "Neutrophils", "NK_Cells", "CD4+_T_Cells", "CD8+_T_Cells", "B_Cells") %>% str_replace_all("_", " ")
plotdata$celltype2 = plotdata$celltype %>% str_replace_all("_", " ")
plotdata[,celltype2:=factor(celltype2, levels = levels2)]
plotdata[,.N,.(celltype,celltype2)]

p_spearman = ggplot(plotdata, aes(celltype2, value^2, alpha = variable, fill = celltype)) +
  theme_minimal(base_size = 16) +
  labs(alpha = "")+
  geom_col(position = "dodge") +
  facet_grid(~species2) +
  # scale_alpha_manual(values = c(1, 0.4)) +
  theme(legend.position = "top")+
  scale_y_continuous(breaks = pretty_breaks(8))+
  scale_fill_manual(values = NogpaletteS3) +
  guides(fill = "none") + scale_alpha_manual(values  = c(0.5, 1), labels = c(bquote(R[before~humanization]^2), bquote(R[humanized]^2))) + ylab("Correlation gene expresion hamster vs. human") + xlab("") + scale_x_discrete(guide = guide_axis(angle = -45))

p_spearman
