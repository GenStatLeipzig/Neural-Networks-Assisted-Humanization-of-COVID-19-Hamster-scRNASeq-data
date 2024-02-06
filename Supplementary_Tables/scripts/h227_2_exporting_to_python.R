require(toolboxH)
require(here)
require(Seurat)


# load hamster MA ----
hamsterMA = readRDS(here("R/results/h226_1_seurat_hamsterMA_humanNames.rds"))
hamsterMA_cellanno = hamsterMA@meta.data %>% data.table(keep.rownames = T)
hamsterMA_cellanno[,.N, .(lowQC_or_doublet_cell,lowQC_or_doublet_cluster)]


hamsterPR = readRDS(here("R/results/h226_1_seurat_hamsterPR_humanNames.rds"))
hamsterPR
hamsterPR_cellanno = hamsterPR@meta.data %>% data.table(keep.rownames = T)
hamsterPR_cellanno[,.N, .(lowQC_or_doublet_cell,lowQC_or_doublet_cluster)]


human = readRDS(here("R/results/h220_5_qc_human_singlestudy.rds"))
human

DefaultAssay(human) = "RNA"
human= DietSeurat(human, layers =  c("counts", "data") , assays = "RNA")
human$species = "human"
human_cellanno = human@meta.data %>% data.table(keep.rownames = T)
human_cellanno[,.N, .(lowQC_or_doublet_cell,lowQC_or_doublet_cluster)]


## Check overlaps ----
qlist1 = venn3(rownames(hamsterMA),
               rownames(hamsterPR),
               rownames(human))


qlist_humanMA = venn2(rownames(hamsterMA),

               rownames(human))

qlist_humanPR = venn2(rownames(hamsterPR),

               rownames(human))



goodcells = fread(here("R/results/h225_5_qc_integrated_rpca_cellanno_FILTERED.txt.gz"))


goodcells[, .N, .(lowQC_or_doublet_cell,lowQC_or_doublet_cluster,lowQC_or_doublet_clusterIntegrated,QC_bad_2_exclude)]

qlist_goodcells = venn4(colnames(hamsterMA),
               colnames(hamsterPR),
               colnames(human),
               goodcells$rn)


qlist4b = venn3(colnames(hamsterMA),
               colnames(hamsterPR),
               colnames(human))


## merge human MA cells ----

humanMA = merge(human[qlist_humanMA$q1, colnames(human) %in% goodcells$rn],
                hamsterMA[qlist_humanMA$q1, colnames(hamsterMA) %in% goodcells$rn], merge.data = TRUE)
humanMA

## save  human MA cells ----
gx_matrix_humanMA = humanMA@assays$RNA$counts %>% as.data.table(keep.rownames = T)
hh(gx_matrix_humanMA)
fwrite(gx_matrix_humanMA, here("R/results/h227_2_humanMA_rpca_FILTERED_JOINED_gxmatrix.txt.gz"))


cellanno_humanMA = humanMA@meta.data %>% as.data.table(keep.rownames =T)
cellanno_humanMA
cellanno_humanMA[,.N, species]
fwrite(cellanno_humanMA, here("R/results/h227_2_humanMA_rpca_FILTERED_JOINED_cellanno.txt.gz"))

featureanno_humanMA = humanMA@assays$RNA@meta.features %>% as.data.table(keep.rownames =T)
featureanno_humanMA
fwrite(featureanno_humanMA, here("R/results/h227_2_humanMA_rpca_FILTERED_JOINED_featureanno.txt.gz"))

stopifnot(identical(featureanno_humanMA$rn, gx_matrix_humanMA$rn))

saveRDS(humanMA, here("R/results/h227_2_humanMA_rpca_FILTERED_JOINED_seurat.RDS"))


rm(gx_matrix_humanMA, humanMA, hamsterMA)
## merge human PR cells ----

humanPR = merge(human[qlist_humanPR$q1, colnames(human) %in% goodcells$rn],
                hamsterPR[qlist_humanPR$q1, colnames(hamsterPR) %in% goodcells$rn], merge.data = TRUE)
humanPR

## save  human PR cells ----
gx_matrix_humanPR = humanPR@assays$RNA$counts %>% as.data.table(keep.rownames = T)
hh(gx_matrix_humanPR)
fwrite(gx_matrix_humanPR, here("R/results/h227_2_humanPR_rpca_FILTERED_JOINED_gxmatrix.txt.gz"))


cellanno_humanPR = humanPR@meta.data %>% as.data.table(keep.rownames =T)
cellanno_humanPR
cellanno_humanPR[,.N, species]
fwrite(cellanno_humanPR, here("R/results/h227_2_humanPR_rpca_FILTERED_JOINED_cellanno.txt.gz"))

featureanno_humanPR = humanPR@assays$RNA@meta.features %>% as.data.table(keep.rownames =T)
featureanno_humanPR
fwrite(featureanno_humanPR, here("R/results/h227_2_humanPR_rpca_FILTERED_JOINED_featureanno.txt.gz"))

stopifnot(identical(featureanno_humanPR$rn, gx_matrix_humanPR$rn))

saveRDS(humanPR, here("R/results/h227_2_humanPR_rpca_FILTERED_JOINED_seurat.RDS"))

## finalize script ----
finalizeSkript()


## R version 4.2.3 (2023-03-15 ucrt)
## Platform: x86_64-w64-mingw32/x64 (64-bit)
## Running under: Windows 10 x64 (build 19045)
##
## Matrix products: default
##
## locale:
## [1] LC_COLLATE=German_Germany.utf8  LC_CTYPE=German_Germany.utf8
## [3] LC_MONETARY=German_Germany.utf8 LC_NUMERIC=C
## [5] LC_TIME=German_Germany.utf8
##
## attached base packages:
## [1] stats     graphics  grDevices utils     datasets  methods   base
##
## other attached packages:
##  [1] Seurat_4.9.9.9044       SeuratObject_4.9.9.9084 sp_1.6-0
##  [4] here_1.0.1              toolboxH_0.2.17         eulerr_7.0.0
##  [7] testthat_3.1.7          stringr_1.5.0           scales_1.2.1
## [10] readxl_1.4.2            RColorBrewer_1.1-3      png_0.1-8
## [13] fdrtool_1.2.17          R.utils_2.12.2          R.oo_1.25.0
## [16] R.methodsS3_1.8.2       data.table_1.14.8
##
## loaded via a namespace (and not attached):
##   [1] Rtsne_0.16             colorspace_2.1-0       deldir_1.0-6
##   [4] ellipsis_0.3.2         ggridges_0.5.4         rprojroot_2.0.3
##   [7] RcppHNSW_0.4.1         spatstat.data_3.0-1    rstudioapi_0.14
##  [10] leiden_0.4.3           listenv_0.9.0          ggrepel_0.9.3
##  [13] RSpectra_0.16-1        fansi_1.0.4            codetools_0.2-19
##  [16] splines_4.2.3          cachem_1.0.7           knitr_1.42
##  [19] polyclip_1.10-4        spam_2.9-1             jsonlite_1.8.4
##  [22] ica_1.0-3              cluster_2.1.4          uwot_0.1.14
##  [25] spatstat.sparse_3.0-1  sctransform_0.3.5      shiny_1.7.4
##  [28] compiler_4.2.3         httr_1.4.5             Matrix_1.5-3
##  [31] fastmap_1.1.1          lazyeval_0.2.2         cli_3.6.1
##  [34] later_1.3.0            htmltools_0.5.5        tools_4.2.3
##  [37] igraph_1.4.1           dotCall64_1.0-2        gtable_0.3.3
##  [40] glue_1.6.2             reshape2_1.4.4         RANN_2.6.1
##  [43] dplyr_1.1.1            Rcpp_1.0.10            scattermore_1.0
##  [46] cellranger_1.1.0       jquerylib_0.1.4        vctrs_0.6.1
##  [49] nlme_3.1-161           spatstat.explore_3.1-0 progressr_0.13.0
##  [52] lmtest_0.9-40          spatstat.random_3.1-4  xfun_0.38
##  [55] globals_0.16.2         brio_1.1.3             mime_0.12
##  [58] miniUI_0.1.1.1         lifecycle_1.0.3        irlba_2.3.5.1
##  [61] goftest_1.2-3          future_1.32.0          MASS_7.3-58.2
##  [64] zoo_1.8-11             spatstat.utils_3.0-2   promises_1.2.0.1
##  [67] parallel_4.2.3         yaml_2.3.7             gridExtra_2.3
##  [70] reticulate_1.28        pbapply_1.7-0          ggplot2_3.4.2
##  [73] sass_0.4.5             stringi_1.7.12         fastDummies_1.6.3
##  [76] rlang_1.1.0            pkgconfig_2.0.3        matrixStats_0.63.0
##  [79] evaluate_0.20          lattice_0.20-45        tensor_1.5
##  [82] ROCR_1.0-11            purrr_1.0.1            patchwork_1.1.2
##  [85] htmlwidgets_1.6.2      cowplot_1.1.1          tidyselect_1.2.0
##  [88] parallelly_1.35.0      RcppAnnoy_0.0.20       plyr_1.8.8
##  [91] magrittr_2.0.3         R6_2.5.1               generics_0.1.3
##  [94] pillar_1.9.0           fitdistrplus_1.1-8     abind_1.4-5
##  [97] survival_3.5-0         tibble_3.2.1           future.apply_1.10.0
## [100] KernSmooth_2.23-20     utf8_1.2.3             spatstat.geom_3.1-0
## [103] plotly_4.10.1          rmarkdown_2.21         grid_4.2.3
## [106] digest_0.6.31          xtable_1.8-4           tidyr_1.3.0
## [109] httpuv_1.6.9           munsell_0.5.0          viridisLite_0.4.2
## [112] bslib_0.4.2
