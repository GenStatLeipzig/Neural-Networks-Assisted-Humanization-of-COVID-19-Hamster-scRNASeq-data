
.libPaths("~/rpackages/angmar/")
.libPaths()
require(toolboxH)

require(here)
require(Seurat)
require(ggplot2)
require(ggthemes)
require(scales)
require(dplyr)
require(tidyr)
require(patchwork)

require(harmony)

require(future) #https://satijalab.org/seurat/archive/v3.0/future_vignette.html
options(future.globals.maxSize= 30*1024^3) # to prevent the error
future.seed=TRUE # to care for proper random seeds
# The total size of the 15 globals exported for future expression ('FUN()') is 12.32 GiB.. This exceeds the maximum allowed size of 500.00 MiB (option 'future.globals.maxSize'). The three largest globals are 'object' (12.26 GiB of class 'numeric'), 'split.cells' (54.18 MiB of class 'list') and 'rowVars' (226.99 KiB of class 'function')
plan("multisession", workers = 1)
# plan("sequential")
plan()



set.seed(2312)

# load mesaur
hamsterMA = readRDS(here("R/data/seu_blood_new_combined_integrated_annotated.rds"))
hamsterMA = DietSeurat(hamsterMA)
hamsterMA$dataset = "hamsterMA"
hamsterMA

## max expressed pr
mesaur_maxexpressed = AverageExpression(hamsterMA, assays = "SCT")$SCT %>% apply(.,1, max) %>% as.data.table(keep.rownames = T) %>% setnames(., c("seurat_name", "max_per_cluster"))
mesaur_maxexpressed

anno_mesaur = data.table(seurat_name = rownames(hamsterMA),
                         gene_name = rownames(hamsterMA) # modify gene_name if not standard names were used in seuratobject
                         )
anno_mesaur[grep("ENSMAUG", seurat_name)] # category  orthologue_by := '01:gene ID'] not relevant here

# Load gene annotation accross species ----
orthologuesm_pre = fread(here("R/results/h216_2_nXm_orthologues_mus_rat_mesaur.txt.gz"), na.strings = c("", "NA"))


orthologuesm = orthologuesm_pre#[species =="mesaur"]
orthologuesm[, species := factor(species, levels = c("mesaur",
                                                     "rat",
                                                     "mus"))]


orthologuesm[ species=="mesaur" & is.na(human_name)==F, mergeid2 := orthologue_genename]
anno_mesaur$mergeid2 = anno_mesaur$seurat_name

orthologuesm2 = merge(orthologuesm %>% unique(), (anno_mesaur[is.na(mergeid2)==F,.(mergeid2,seurat_name2= seurat_name)] %>% unique()), by ="mergeid2", all.x = T,allow.cartesian=TRUE) %>% unique()

orthologuesm2[is.na(seurat_name2)==F, orthologue_by := '02:gene name of mesaur']

orthologuesm2[, seurat_name := seurat_name2]

orthologuesm2[, uniqueN(seurat_name, na.rm = T), orthologue_by]
orthologuesm2[, uniqueN(human_name, na.rm = T), orthologue_by]


qlist68445b=venn2(anno_mesaur$seurat_name, orthologuesm2$seurat_name)

## merge via human WITH support from rat or mus----
orthologuesm2[is.na(human_name)==F &
                is.na(orthologue_genename)==F &
                species %in% c("mus", "rat")
                , mergeid3 := toupper(paste(orthologue_genename, human_name))]

anno_mesaur[is.na(gene_name)==F, mergeid3 := toupper(paste(gene_name, gene_name))]

orthologuesm3 = merge(orthologuesm2 %>% unique(), (anno_mesaur[is.na(mergeid3)==F,.(mergeid3, seurat_name3 = seurat_name)] %>% unique()), by = "mergeid3", all.x = T,allow.cartesian=TRUE) %>% unique()


orthologuesm3[is.na(seurat_name3)==F & is.na(seurat_name)==T, orthologue_by := '03:same name uppercase in human and [rat or mus]']

orthologuesm3[, seurat_name := ifelse(is.na(seurat_name), seurat_name3, seurat_name)]

orthologuesm3[, uniqueN(seurat_name, na.rm = T), orthologue_by]
orthologuesm3[, uniqueN(human_name, na.rm = T), orthologue_by]


qlist68445c=venn2(anno_mesaur$seurat_name, orthologuesm3$seurat_name)


## merge via human without support from rat or mus----
orthologuesm3[is.na(human_name)==F
              , mergeid4 := toupper(paste(human_name))]

anno_mesaur[is.na(gene_name)==F, mergeid4 := toupper(paste(gene_name))]

orthologuesm4 = merge(orthologuesm3 %>% unique(), (anno_mesaur[is.na(mergeid4)==F,.(mergeid4, seurat_name4 = seurat_name)] %>% unique()), by = "mergeid4", all.x = T,allow.cartesian=TRUE) %>% unique()


orthologuesm4[is.na(seurat_name4)==F & is.na(seurat_name)==T, orthologue_by := '04:same name uppercase in human']

orthologuesm4[, seurat_name := ifelse(is.na(seurat_name), seurat_name4, seurat_name)]

orthologuesm4[, uniqueN(seurat_name, na.rm = T), orthologue_by]
orthologuesm4[, uniqueN(human_name, na.rm = T), orthologue_by]


qlist68445d=venn2(anno_mesaur$seurat_name, orthologuesm4$seurat_name)



# # make unique assignment preferring cluster-specific higher expressed genes ----


orthologuesm4[, mesaur_max_per_cluster := mesaur_maxexpressed[match_hk(orthologuesm4$seurat_name, mesaur_maxexpressed$seurat_name),max_per_cluster]]


setorder(orthologuesm4,
         orthologue_by,


         -orthologue_confidence,

         -mesaur_max_per_cluster,
         human_name,
         na.last = T
)

orthologuesm4[is.na(human_name)==F &is.na(seurat_name)==F][allDuplicatedEntries(seurat_name)]  # multiple human genes for same hamster gene

orthologues_mesaur_pre = orthologuesm4[duplicated(seurat_name)==F ]
orthologues_mesaur_pre[, uniqueN(seurat_name), orthologue_by]
orthologues_mesaur_pre[is.na(human_name)==F][allDuplicatedEntries(human_name)]

qlist10 = venn2(orthologues_mesaur_pre$seurat_name, rownames(hamsterMA))

setorder(orthologues_mesaur_pre,
         orthologue_by,

         -orthologue_confidence,
         species,
         -human_max_per_cluster,
         seurat_name,
         na.last = T
)

orthologues_mesaur_pre[is.na(human_name)==F & is.na(seurat_name )==F][allDuplicatedEntries(human_name)] # multiple human genes for same pr gene

orthologues_mesaur = orthologues_mesaur_pre[duplicated(human_name)==F ]
orthologues_mesaur[, uniqueN(seurat_name), orthologue_by]
orthologues_mesaur[, uniqueN(human_name), orthologue_by]


orthologues_mesaur[is.na(human_name)==F][allDuplicatedEntries(human_name)]

qlist10b = venn2(orthologues_mesaur$seurat_name, rownames(hamsterMA))

anno_mesaur[, human_orthologue_name := orthologues_mesaur[match_hk(anno_mesaur$seurat_name, orthologues_mesaur$seurat_name),human_name]]

anno_mesaur[, human_orthologue_description := orthologues_mesaur[match_hk(anno_mesaur$seurat_name, orthologues_mesaur$seurat_name),description]]

anno_mesaur[, human_orthologue_id := orthologues_mesaur[match_hk(anno_mesaur$seurat_name, orthologues_mesaur$seurat_name),human_id]]

anno_mesaur[, human_orthologue_confidence := orthologues_mesaur[match_hk(anno_mesaur$seurat_name, orthologues_mesaur$seurat_name),orthologue_confidence]]

anno_mesaur[, human_orthologue_by := orthologues_mesaur[match_hk(anno_mesaur$seurat_name, orthologues_mesaur$seurat_name),orthologue_by ]]

anno_mesaur[, .N,human_orthologue_by][order(human_orthologue_by)]

# # table for manual check anno -----

duplicatefailure = orthologues_mesaur_pre[seurat_name %nin% orthologues_mesaur$seurat_name]
duplicatefailure[, seurat_name_used := anno_mesaur[match_hk(duplicatefailure$human_name, anno_mesaur$human_orthologue_name,makeunique = T, importcol = anno_mesaur$seurat_name), seurat_name]]
duplicatefailure[, orthologue_by_used := anno_mesaur[match_hk(duplicatefailure$human_name, anno_mesaur$human_orthologue_name,makeunique = T, importcol = anno_mesaur$seurat_name), human_orthologue_by]]
duplicatefailure2 = duplicatefailure[, .(human_id, human_name, seurat_name, orthologue_by, seurat_name_used, orthologue_by_used)]

no_orthologue = anno_mesaur[is.na(human_orthologue_name) & seurat_name %nin% duplicatefailure$seurat_name & is.na(seurat_name)==F]

no_orthologue2 = no_orthologue[, .(seurat_name, human_orthologue_name)][order(seurat_name)]
no_orthologue2
duplicatefailure2


# # save ----
WriteXLS_hk(c("no_orthologue2", "duplicatefailure2"), here("R/results/h217_2_failed_orthologues_mesaur.xlsx"), AdjWidth = T)


# # save ----
fwrite(anno_mesaur, here("R/results/h217_2_orthologues_measaur_used.txt.gz"), sep = "\t")

fwrite(orthologues_mesaur, here("R/results/h217_2_helperfile_orthologues_mesaur.txt.gz"), sep = "\t")

fwrite(orthologuesm4, here("R/results/h217_2_nXm_orthologues_mesaur.txt.gz"), sep = "\t")





finalizeSkript()

## R version 4.2.2 (2022-10-31)
## Platform: x86_64-suse-linux-gnu (64-bit)
## Running under: openSUSE Leap 15.3
##
## Matrix products: default
## BLAS:   /usr/lib64/R/lib/libRblas.so
## LAPACK: /usr/lib64/R/lib/libRlapack.so
##
## locale:
##  [1] LC_CTYPE=de_DE.UTF-8       LC_NUMERIC=C
##  [3] LC_TIME=de_DE.UTF-8        LC_COLLATE=de_DE.UTF-8
##  [5] LC_MONETARY=de_DE.UTF-8    LC_MESSAGES=de_DE.UTF-8
##  [7] LC_PAPER=de_DE.UTF-8       LC_NAME=C
##  [9] LC_ADDRESS=C               LC_TELEPHONE=C
## [11] LC_MEASUREMENT=de_DE.UTF-8 LC_IDENTIFICATION=C
##
## attached base packages:
## [1] stats     graphics  grDevices utils     datasets  methods   base
##
## other attached packages:
##  [1] future_1.32.0           harmony_0.1.0           Rcpp_1.0.10
##  [4] patchwork_1.1.2         tidyr_1.3.0             dplyr_1.1.1
##  [7] ggthemes_4.2.4          ggplot2_3.4.2           Seurat_4.9.9.9042
## [10] SeuratObject_4.9.9.9084 sp_1.6-0                here_1.0.1
## [13] toolboxH_0.2.17         eulerr_6.1.1            testthat_3.1.4
## [16] stringr_1.5.0           scales_1.2.1            readxl_1.4.1
## [19] RColorBrewer_1.1-3      png_0.1-8               fdrtool_1.2.17
## [22] R.utils_2.12.0          R.oo_1.25.0             R.methodsS3_1.8.2
## [25] data.table_1.14.8
##
## loaded via a namespace (and not attached):
##   [1] Rtsne_0.16             colorspace_2.1-0       deldir_1.0-6
##   [4] ellipsis_0.3.2         ggridges_0.5.4         rprojroot_2.0.3
##   [7] RcppHNSW_0.4.1         spatstat.data_3.0-1    rstudioapi_0.14
##  [10] leiden_0.4.3           listenv_0.9.0          ggrepel_0.9.3
##  [13] RSpectra_0.16-1        fansi_1.0.4            codetools_0.2-18
##  [16] splines_4.2.2          cachem_1.0.7           knitr_1.42
##  [19] polyclip_1.10-4        spam_2.9-1             jsonlite_1.8.4
##  [22] ica_1.0-3              cluster_2.1.4          uwot_0.1.14
##  [25] spatstat.sparse_3.0-1  sctransform_0.3.5      shiny_1.7.4
##  [28] compiler_4.2.2         httr_1.4.5             Matrix_1.5-3
##  [31] fastmap_1.1.1          lazyeval_0.2.2         cli_3.6.1
##  [34] later_1.3.0            htmltools_0.5.5        tools_4.2.2
##  [37] igraph_1.4.1           dotCall64_1.0-2        gtable_0.3.3
##  [40] glue_1.6.2             reshape2_1.4.4         RANN_2.6.1
##  [43] scattermore_1.0        cellranger_1.1.0       jquerylib_0.1.4
##  [46] vctrs_0.6.1            nlme_3.1-159           spatstat.explore_3.1-0
##  [49] progressr_0.13.0       lmtest_0.9-40          spatstat.random_3.1-4
##  [52] xfun_0.38              globals_0.16.2         brio_1.1.3
##  [55] mime_0.12              miniUI_0.1.1.1         lifecycle_1.0.3
##  [58] irlba_2.3.5.1          WriteXLS_6.4.0         goftest_1.2-3
##  [61] MASS_7.3-58.1          zoo_1.8-11             spatstat.utils_3.0-2
##  [64] promises_1.2.0.1       parallel_4.2.2         yaml_2.3.7
##  [67] gridExtra_2.3          reticulate_1.28        pbapply_1.7-0
##  [70] sass_0.4.5             stringi_1.7.12         fastDummies_1.6.3
##  [73] rlang_1.1.0            pkgconfig_2.0.3        matrixStats_0.63.0
##  [76] evaluate_0.20          lattice_0.20-45        tensor_1.5
##  [79] ROCR_1.0-11            purrr_1.0.1            htmlwidgets_1.6.2
##  [82] cowplot_1.1.1          tidyselect_1.2.0       parallelly_1.35.0
##  [85] RcppAnnoy_0.0.20       plyr_1.8.8             magrittr_2.0.3
##  [88] R6_2.5.1               generics_0.1.3         withr_2.5.0
##  [91] pillar_1.9.0           fitdistrplus_1.1-8     abind_1.4-5
##  [94] survival_3.4-0         tibble_3.2.1           future.apply_1.10.0
##  [97] KernSmooth_2.23-20     utf8_1.2.3             spatstat.geom_3.1-0
## [100] plotly_4.10.1          rmarkdown_2.21         grid_4.2.2
## [103] digest_0.6.31          xtable_1.8-4           httpuv_1.6.9
## [106] munsell_0.5.0          viridisLite_0.4.1      bslib_0.4.2
