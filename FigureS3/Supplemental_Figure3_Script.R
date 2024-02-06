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

sessionInfo()
# R version 4.3.0 (2023-04-21 ucrt)
# Platform: x86_64-w64-mingw32/x64 (64-bit)
# Running under: Windows 10 x64 (build 19045)
# 
# Matrix products: default
# 
# 
# locale:
#   [1] LC_COLLATE=German_Germany.utf8  LC_CTYPE=German_Germany.utf8    LC_MONETARY=German_Germany.utf8
# [4] LC_NUMERIC=C                    LC_TIME=German_Germany.utf8    
# 
# time zone: Europe/Berlin
# tzcode source: internal
# 
# attached base packages:
#   [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] ggplot2_3.4.2      patchwork_1.1.2    here_1.0.1         toolboxH_0.2.17    eulerr_7.0.0       testthat_3.1.10   
# [7] stringr_1.5.0      scales_1.2.1       readxl_1.4.3       RColorBrewer_1.1-3 png_0.1-8          fdrtool_1.2.17    
# [13] R.utils_2.12.2     R.oo_1.25.0        R.methodsS3_1.8.2  data.table_1.14.8 
# 
# loaded via a namespace (and not attached):
#   [1] gtable_0.3.3          dplyr_1.1.2           compiler_4.3.0        BiocManager_1.30.21.1 brio_1.1.3           
# [6] tidyselect_1.2.0      Rcpp_1.0.11           R6_2.5.1              generics_0.1.3        tibble_3.2.1         
# [11] munsell_0.5.0         rprojroot_2.0.3       pillar_1.9.0          rlang_1.1.1           utf8_1.2.3           
# [16] stringi_1.7.12        cli_3.6.1             withr_2.5.0           magrittr_2.0.3        grid_4.3.0           
# [21] rstudioapi_0.15.0     lifecycle_1.0.3       vctrs_0.6.3           glue_1.6.2            cellranger_1.1.0     
# [26] fansi_1.0.4           colorspace_2.1-0      tools_4.3.0           pkgconfig_2.0.3   
