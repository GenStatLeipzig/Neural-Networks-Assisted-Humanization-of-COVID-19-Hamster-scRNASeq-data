#Supplemental Figure 1

require(toolboxH)
require(here)
require(patchwork)
require(ggplot2)
require(UpSetR)

#S1
plotdata = fread("Data/h224_3_overlap_genes_all_datasets_plotdat.txt.gz")
upset(plotdata,nsets = ncol(plotdata)-1,order.by =c("freq" ,"degree"), decreasing = c(TRUE, TRUE), set_size.scale_max =35000, text.scale = 1.3, set_size.show = TRUE,mainbar.y.label = "Overlap Genes", nintersects = ncol(plotdata))


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
#   [1] UpSetR_1.4.0       ggplot2_3.4.2      patchwork_1.1.2    here_1.0.1         toolboxH_0.2.17    eulerr_7.0.0      
# [7] testthat_3.1.10    stringr_1.5.0      scales_1.2.1       readxl_1.4.3       RColorBrewer_1.1-3 png_0.1-8         
# [13] fdrtool_1.2.17     R.utils_2.12.2     R.oo_1.25.0        R.methodsS3_1.8.2  data.table_1.14.8 
# 
# loaded via a namespace (and not attached):
#   [1] gtable_0.3.3          dplyr_1.1.2           compiler_4.3.0        BiocManager_1.30.21.1 brio_1.1.3           
# [6] Rcpp_1.0.11           tidyselect_1.2.0      gridExtra_2.3         plyr_1.8.8            R6_2.5.1             
# [11] generics_0.1.3        tibble_3.2.1          rprojroot_2.0.3       munsell_0.5.0         pillar_1.9.0         
# [16] rlang_1.1.1           utf8_1.2.3            stringi_1.7.12        cli_3.6.1             withr_2.5.0          
# [21] magrittr_2.0.3        grid_4.3.0            rstudioapi_0.15.0     lifecycle_1.0.3       vctrs_0.6.3          
# [26] glue_1.6.2            cellranger_1.1.0      fansi_1.0.4           colorspace_2.1-0      tools_4.3.0          
# [31] pkgconfig_2.0.3  
