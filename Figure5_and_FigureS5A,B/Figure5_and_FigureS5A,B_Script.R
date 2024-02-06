#Figure5 A-C and Supplemental Figure 5A-B


require(toolboxH) #devtools::install_github('holgerman/toolboxH'); https://github.com/holgerman/toolboxH
library(ggplot2)
library(dplyr)
library(Seurat)
library(cowplot)
library(patchwork)
library(ggrepel)
require(ggthemes) 
library(reshape2)
library(ggforce)
library(scales)


#Load Data
Blood_merged <- readRDS(".../h228_3_qc_integrated_rpca_cellanno_FILTEREDqc_FILTEREDpureCells.RDS")



#5A

eta1_tab2 = fread("Data/h245_7_eta1_tab.txt")

eta1_tab2[, contrast1v2 := factor(contrast1,
                                  levels = c("3", "4u5", "7", "d2", "d3", "d5", "e14","ld_D2", "ld_D3", "hd_D2", "hd_D3"),
                                  labels = c("WHO3", "WHO4&5", "WHO7", "2 dpi", "3 dpi", "5 dpi", "14 dpi","2 dpi ld", "3 dpi ld", "2 dpi hd", "3 dpi hd"))]




eta1_tab2[, celltype2 := factor(celltype, levels =  c("Classical_Monocytes", "Non_Classical_Monocytes", "Immature Neutrophils 1",  "Neutrophils", "NK_Cells", "CD4+_T_Cells", "CD8+_T_Cells", "B_Cells","DC", "Platelet",  "Plasmablasts","Megakaryocyte", "Immature Neutrophils 2", "Proliferating_Cells") %>% rev())]


p_eta1v2 = ggplot(eta1_tab2, aes(contrast1v2, celltype2,label = label, fill =eta1)) + geom_tile() +
  theme_minimal(base_size = 14) +
  facet_grid(.~species, scales = 'free', space = "free") + geom_text() + scale_fill_gradient_tableau(palette = "Red-Gold") + ylab("") + xlab("") + guides(fill = "none") + ggtitle("Percent transcriptome changed compared with control")
p_eta1v2




#5B

data_fig5b = fread("Data/h265_8_plotdata_FDR20_dotplot_effectsizes.txt.gz")

Cellorder5b = c("Classical\nMonocytes", "Non Classical\nMonocytes","DC", "Immature\nNeutrophils 1", "Neutrophils", "NK Cells", "CD4+ T Cells", "CD8+ T Cells", "B Cells")


data_fig5b$human_contrast1v2[data_fig5b$human_contrast1 =="3"] <- "WOS 3"
data_fig5b$human_contrast1v2[data_fig5b$human_contrast1 =="7"] <- "WOS 7"
data_fig5b$human_contrast1v2[data_fig5b$human_contrast1 =="4u5"] <- "WOS 4&5"

data_fig5b$hamster_contrast1v3[data_fig5b$hamster_contrast1v2 =="d2"] <- "c2 dpi"
data_fig5b$hamster_contrast1v3[data_fig5b$hamster_contrast1v2 =="d3"] <- "b3 dpi"
data_fig5b$hamster_contrast1v3[data_fig5b$hamster_contrast1v2 =="d5"] <- "a5 dpi"
data_fig5b$hamster_contrast1v3[data_fig5b$hamster_contrast1v2 =="e14"] <- "14 dpi"

data_fig5b$hamster_contrast1v3[data_fig5b$hamster_contrast1v2 =="ld_D2"] <- "d2 dpi ld"
data_fig5b$hamster_contrast1v3[data_fig5b$hamster_contrast1v2 =="ld_D3"] <- "c3 dpi ld"
data_fig5b$hamster_contrast1v3[data_fig5b$hamster_contrast1v2 =="hd_D2"] <- "b2 dpi hd"
data_fig5b$hamster_contrast1v3[data_fig5b$hamster_contrast1v2 =="hd_D3"] <- "a3 dpi hd"

data_fig5b$celltype3 <- factor(data_fig5b$celltype3, levels = c(Cellorder5b))


p2=ggplot(data_fig5b[human_contrast1v2 != "WOS 4&5" & n_intersectgenes>0], aes( human_contrast1v2, hamster_contrast1v3, col = -intersectgenes_sign_samedir_proz , size = n_signif_genes , alpha = p_binom<=0.05)) + geom_point() +
  facet_grid(hamster~celltype3, space = "free", scales = "free") +
  scale_size_continuous(range = c(2,9), breaks = c(5, 50, 100,250, 400, 1000))+
  theme_minimal(base_size = 14) +
  labs(size = "Number overlapping significant\ngenes (FDR20%)",
       color = "Percent human effect direction\n same in hamster",
       alpha = "Significant difference from 50% of\n same effect direction in hamster") +
  theme(legend.position = "top",legend.direction = "vertical") +
  scale_color_gradient2_tableau(guide = "legend", breaks = (c(0:5, 2.5)/-5)%>% sort(), labels = (c(0:5, 2.5)/5)%>% sort() %>% rev() *100,limits = c(-1,0)) +
  guides(col = guide_legend(nrow = 1,override.aes = list(size=5)),
         size = guide_legend(ncol = 2),
         alpha = guide_legend(override.aes = list(size=5))) +
  ylab("") +xlab("")

p2


#S5A
p3=ggplot(data_fig5b[human_contrast1v2 == "WOS 4&5" & n_intersectgenes>0], aes( human_contrast1v2, hamster_contrast1v3, col = -intersectgenes_sign_samedir_proz , size = n_signif_genes , alpha = p_binom<=0.05)) + geom_point() +
  facet_grid(hamster~celltype3, space = "free", scales = "free") +
  scale_size_continuous(range = c(2,9), breaks = c(5, 50, 100,250, 400, 1000))+
  theme_minimal(base_size = 14) +
  labs(size = "Number overlapping significant\ngenes (FDR20%)",
       color = "Percent human effect direction\n same in hamster",
       alpha = "Significant difference from 50% of\n same effect direction in hamster") +
  theme(legend.position = "top",legend.direction = "vertical") +
  scale_color_gradient2_tableau(guide = "legend", breaks = (c(0:5, 2.5)/-5)%>% sort(), labels = (c(0:5, 2.5)/5)%>% sort() %>% rev() *100,limits = c(-1,0)) +
  guides(col = guide_legend(nrow = 1,override.aes = list(size=5)),
         size = guide_legend(ncol = 2),
         alpha = guide_legend(override.aes = list(size=5))) +
  ylab("") +xlab("")

p3


#5C
data_fig5c = fread("Data/h245_7_inputdata_plot_h245_7_overlapping_pathways_no4u5_smaller.txt.gz")

data_fig5c$human2[data_fig5c$human =="3"] <- "WOS 3"
data_fig5c$human2[data_fig5c$human =="7"] <- "WOS 7"
data_fig5c$human2[data_fig5c$human =="4u5"] <- "WOS 4&5"

data_fig5c$hamster2[data_fig5c$hamster =="d2"] <- "2 dpi"
data_fig5c$hamster2[data_fig5c$hamster =="d3"] <- "3 dpi"
data_fig5c$hamster2[data_fig5c$hamster =="d5"] <- "5 dpi"
data_fig5c$hamster2[data_fig5c$hamster =="e14"] <- "14 dpi"

data_fig5c$hamster2[data_fig5c$hamster =="ld_D2"] <- "2 dpi ld"
data_fig5c$hamster2[data_fig5c$hamster =="ld_D3"] <- "3 dpi ld"
data_fig5c$hamster2[data_fig5c$hamster =="hd_D2"] <- "2 dpi hd"
data_fig5c$hamster2[data_fig5c$hamster =="hd_D3"] <- "3 dpi hd"

data_fig5c[,human2:=factor(human2)]
dput(data_fig5c$hamster2 %>% unique())
data_fig5c[,hamster2:=factor(hamster2, levels = c("2 dpi", "3 dpi", "5 dpi", "2 dpi ld", "3 dpi ld", "2 dpi hd", "3 dpi hd"))]


p_duplipathways_no4u5 = ggplot(data_fig5c[human2 !="WOS 4&5"], aes(hamster2, human2,  col = variable, size= value)) +
  geom_point( alpha = 0.1) +
  geom_point( alpha = 1, pch = 1) +
  facet_grid(celltype +term_name~., scales= "free_x",switch  = "x") + theme_minimal(base_size = 14) +
  labs(size = "Pathway enrichment factor",
       color = "Species") +
  theme(strip.text.y.left =  element_text(angle = 0, hjust = 0),
        strip.text.y.right =   element_text(angle = 0, hjust = 0),
        legend.position = "top",
        legend.text.align = 0,
        legend.box = "vertical") +
  guides(color = guide_legend(override.aes = list(size=5))) +
  scale_size_continuous(range =  c(1.5,11)) +
  scale_color_manual(values = c("#FB3640", "#0D3B66"))+
  scale_fill_manual(values = c("#FB3640", "#0D3B66"))

p_duplipathways_no4u5


#S5B
p_duplipathways_only4u5 = ggplot(data_fig5c[human2 =="WOS 4&5"], aes(hamster2, human2,  col = variable, size= value)) +
  geom_point( alpha = 0.1) +
  geom_point( alpha = 1, pch = 1) +
  facet_grid(celltype +term_name~., scales= "free_x",switch  = "x") + theme_minimal(base_size = 14) +
  labs(size = "Pathway enrichment factor",
       color = "Species") +
  theme(strip.text.y.left =  element_text(angle = 0, hjust = 0),
        strip.text.y.right =   element_text(angle = 0, hjust = 0),
        legend.position = "top",
        legend.text.align = 0,
        legend.box = "vertical") +
  guides(color = guide_legend(override.aes = list(size=5))) +
  scale_size_continuous(range =  c(1.5,11)) +
  scale_color_manual(values = c("#FB3640", "#0D3B66"))+
  scale_fill_manual(values = c("#FB3640", "#0D3B66"))

p_duplipathways_only4u5

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
#   [1] ggforce_0.4.1           reshape2_1.4.4          ggthemes_4.2.4          ggrepel_0.9.3           patchwork_1.1.2        
# [6] cowplot_1.1.1           Seurat_4.9.9.9058       SeuratObject_4.9.9.9091 sp_2.0-0                dplyr_1.1.2            
# [11] ggplot2_3.4.2           toolboxH_0.2.17         eulerr_7.0.0            testthat_3.1.10         stringr_1.5.0          
# [16] scales_1.2.1            readxl_1.4.3            RColorBrewer_1.1-3      png_0.1-8               fdrtool_1.2.17         
# [21] R.utils_2.12.2          R.oo_1.25.0             R.methodsS3_1.8.2       data.table_1.14.8      
# 
# loaded via a namespace (and not attached):
#   [1] deldir_1.0-9           pbapply_1.7-2          gridExtra_2.3          rlang_1.1.1            magrittr_2.0.3        
# [6] RcppAnnoy_0.0.21       spatstat.geom_3.2-4    matrixStats_1.0.0      ggridges_0.5.4         compiler_4.3.0        
# [11] vctrs_0.6.3            pkgconfig_2.0.3        fastmap_1.1.1          ellipsis_0.3.2         utf8_1.2.3            
# [16] promises_1.2.0.1       purrr_1.0.1            jsonlite_1.8.7         goftest_1.2-3          later_1.3.1           
# [21] tweenr_2.0.2           spatstat.utils_3.0-3   irlba_2.3.5.1          parallel_4.3.0         cluster_2.1.4         
# [26] R6_2.5.1               ica_1.0-3              spatstat.data_3.0-1    stringi_1.7.12         reticulate_1.30       
# [31] parallelly_1.36.0      brio_1.1.3             lmtest_0.9-40          scattermore_1.2        cellranger_1.1.0      
# [36] Rcpp_1.0.11            tensor_1.5             future.apply_1.11.0    zoo_1.8-12             sctransform_0.3.5     
# [41] httpuv_1.6.11          Matrix_1.5-4           splines_4.3.0          igraph_1.5.0.1         tidyselect_1.2.0      
# [46] abind_1.4-5            rstudioapi_0.15.0      spatstat.random_3.1-5  spatstat.explore_3.2-1 codetools_0.2-19      
# [51] miniUI_0.1.1.1         listenv_0.9.0          lattice_0.21-8         tibble_3.2.1           plyr_1.8.8            
# [56] shiny_1.7.4.1          withr_2.5.0            ROCR_1.0-11            Rtsne_0.16             future_1.33.0         
# [61] fastDummies_1.7.3      survival_3.5-5         polyclip_1.10-4        fitdistrplus_1.1-11    pillar_1.9.0          
# [66] BiocManager_1.30.21.1  KernSmooth_2.23-20     plotly_4.10.2          generics_0.1.3         RcppHNSW_0.4.1        
# [71] munsell_0.5.0          globals_0.16.2         xtable_1.8-4           glue_1.6.2             lazyeval_0.2.2        
# [76] tools_4.3.0            RSpectra_0.16-1        RANN_2.6.1             leiden_0.4.3           dotCall64_1.0-2       
# [81] grid_4.3.0             tidyr_1.3.0            colorspace_2.1-0       nlme_3.1-162           cli_3.6.1             
# [86] spatstat.sparse_3.0-2  spam_2.9-1             fansi_1.0.4            viridisLite_0.4.2      uwot_0.1.16           
# [91] gtable_0.3.3           digest_0.6.33          progressr_0.13.0       farver_2.1.1           htmlwidgets_1.6.2     
# [96] htmltools_0.5.5        lifecycle_1.0.3        httr_1.4.6             mime_0.12              MASS_7.3-58.4   
