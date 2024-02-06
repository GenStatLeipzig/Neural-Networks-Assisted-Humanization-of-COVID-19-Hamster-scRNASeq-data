#Figure 6 and Supplemental Figure 5C

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

#Load Data
Blood_merged <- readRDS(".../h228_3_qc_integrated_rpca_cellanno_FILTEREDqc_FILTEREDpureCells.RDS")

#Figure6
data_fig6 = fread("Data/h265_8_plotdata_top5_genes_same_direction_no4u5.txt.gz")
allcelltypes = unique(data_fig6$celltype)

data_fig6[,hamster_contrast1:=factor(hamster_contrast1, levels = c("d2", "d3", "d5", "e14", "ld_D2", "ld_D3", "hd_D2", "hd_D3"))]
data_fig6[,human_contrast1:=factor(human_contrast1)]
dput(data_fig6$hamster_contrast1 %>% unique())

p_similgeneslist_no4u5 = lapply(allcelltypes, function(mycelltype) {
  # mycelltype = allcelltypes[2]
  myplotdat = data_fig6[celltype ==mycelltype & human_contrast1 != "4u5"  & topgen_samedir_nr_no4u5  <=10]
  
  if(nrow(myplotdat)==0) {message("No overlapping genes found for ",mycelltype );return(NULL)}
  myplot = ggplot(myplotdat, aes(hamster_contrast1 , human_contrast1 ,  col = Species, pch = direction , size= abslogFC)) +
    geom_point( alpha = 0.6) +
    facet_grid(features~., switch  = "x") + theme_minimal(base_size = 18) +
    labs(size = "abs. LogFoldChange",
         color = "Species") +
    theme(strip.text.y.left =  element_text(angle = 0, hjust = 0),
          strip.text.y.right =   element_text(angle = 0, hjust = 0),
          legend.position = "top",
          legend.text.align = 0,
          legend.box = "horizontal") +
    scale_shape_manual(values = c(6,2),drop=FALSE)+
    scale_size_continuous(range =  c(1.5,11), limits = c(0.5,11)) +
    scale_color_manual(values = c("#FB3640", "#0D3B66"))+
    scale_fill_manual(values = c("#FB3640", "#0D3B66")) +
    guides(color = guide_legend(override.aes = list(size=5), ncol = 1),
           pch = guide_legend(override.aes = list(size=5), ncol = 2),
           size = guide_legend(ncol = 1, override.aes = list(pch = 2))
    ) +
    scale_x_discrete(guide = guide_axis(angle = -45), drop=FALSE) +
    ggtitle(mycelltype)
  myplot
}
)
p_similgeneslist_no4u5 = p_similgeneslist_no4u5[sapply(p_similgeneslist_no4u5, is.null)==F]
plot_fig6 = patchwork::wrap_plots(p_similgeneslist_no4u5) + plot_layout(guides = "collect",ncol = 3)
plot_fig6


#Figure S5C
p_similgeneslist_just4u5 = lapply(allcelltypes, function(mycelltype) {
  # mycelltype = allcelltypes[6]
  myplotdat = data_fig6[celltype ==mycelltype & human_contrast1 == "4u5"  &topgen_samedir_nr   <=10]
  levels(myplotdat$human_contrast)
  
  if(nrow(myplotdat)==0) {message("No overlapping genes found for ",mycelltype );return(NULL)}
  myplot = ggplot(myplotdat, aes(hamster_contrast1 , human_contrast1,  col = Species, pch = direction , size= abslogFC)) +
    geom_point( alpha = 0.6) +
    facet_grid(features~., scales= "free_y",switch  = "x") + theme_minimal(base_size = 18) +
    labs(size = "abs. LogFoldChange",
         color = "Species") +
    theme(strip.text.y.left =  element_text(angle = 0, hjust = 0),
          strip.text.y.right =   element_text(angle = 0, hjust = 0),
          legend.position = "top",
          legend.text.align = 0,
          legend.box = "horizontal") +
    scale_shape_manual(values = c(6,2),drop=FALSE)+
    scale_size_continuous(range =  c(1.5,11), limits = c(0.5,11)) +
    scale_color_manual(values = c("#FB3640", "#0D3B66"))+
    scale_fill_manual(values = c("#FB3640", "#0D3B66")) +
    guides(color = guide_legend(override.aes = list(size=5), ncol = 1),
           pch = guide_legend(override.aes = list(size=5), ncol = 2),
           size = guide_legend(ncol = 1, override.aes = list(pch = 2))
    ) +
    scale_x_discrete(guide = guide_axis(angle = -45), drop=FALSE)+
    scale_y_discrete(drop=FALSE) +
    ylab("WOS")+
    xlab("")+
    
    ggtitle(mycelltype)
  myplot
}
)
p_similgeneslist_just4u5 = p_similgeneslist_just4u5[sapply(p_similgeneslist_just4u5, is.null)==F]
plot_supplfig5C = patchwork::wrap_plots(p_similgeneslist_just4u5) + plot_layout(guides = "collect",ncol = 3)
plot_supplfig5C

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
# [51] miniUI_0.1.1.1         listenv_0.9.0          plyr_1.8.8             lattice_0.21-8         tibble_3.2.1          
# [56] shiny_1.7.4.1          withr_2.5.0            ROCR_1.0-11            Rtsne_0.16             future_1.33.0         
# [61] fastDummies_1.7.3      survival_3.5-5         polyclip_1.10-4        fitdistrplus_1.1-11    pillar_1.9.0          
# [66] BiocManager_1.30.21.1  KernSmooth_2.23-20     plotly_4.10.2          generics_0.1.3         RcppHNSW_0.4.1        
# [71] munsell_0.5.0          globals_0.16.2         xtable_1.8-4           glue_1.6.2             lazyeval_0.2.2        
# [76] tools_4.3.0            RSpectra_0.16-1        RANN_2.6.1             leiden_0.4.3           dotCall64_1.0-2       
# [81] grid_4.3.0             tidyr_1.3.0            colorspace_2.1-0       nlme_3.1-162           cli_3.6.1             
# [86] spatstat.sparse_3.0-2  spam_2.9-1             fansi_1.0.4            viridisLite_0.4.2      uwot_0.1.16           
# [91] gtable_0.3.3           digest_0.6.33          progressr_0.13.0       farver_2.1.1           htmlwidgets_1.6.2     
# [96] htmltools_0.5.5        lifecycle_1.0.3        httr_1.4.6             mime_0.12              MASS_7.3-58.4 
