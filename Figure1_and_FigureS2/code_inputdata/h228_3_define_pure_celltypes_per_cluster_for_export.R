require(toolboxH)
require(here)
require(Seurat)

jpeg2 = function(fn, ...) jpeg(fn, quality = 100, unit ="in", res = 150, ...)

goodcells = fread(here("R/results/h225_5_qc_integrated_rpca_cellanno_FILTERED.txt.gz"))
goodcells[1]
require(ggplot2)
require(patchwork)

## update cell names
goodcells[, uniform_name_overview2 := ifelse(celltype =="Immature neutrophil" | cluster_labels_res.0.8 %in% c("Immature Neutrophils_1", "Immature Neutrophils_2"), "Immature neutrophil", uniform_name_overview)]

plotSankey(goodcells[, .(uniform_name_overview, uniform_name_overview2)], "uniform_name_overview2")


goodcells[,.N, .(species, uniform_name_overview, rpca_clusters.v2)][order(species,uniform_name_overview, rpca_clusters.v2)]

p1a = ggplot(goodcells[grepl("^un", uniform_name_overview)==F], aes(rpca_clusters.v2, fill = factor(uniform_name_overview))) + geom_bar(position = "fill") + facet_grid(species~.)
p2a = ggplot(goodcells[grepl("^un", uniform_name_overview)==F], aes(rpca_clusters.v2, fill = factor(uniform_name_overview))) + geom_bar() + facet_grid(species~.)
p1a
p2a


goodcells[, .N, .(timepoint,species, who_per_sample)]
goodcells[, severity := ifelse(species =="human", who_per_sample, timepoint)]
goodcells[, .N, .(species, severity)]

n_goodcells_pur = goodcells[,.N, .(uniform_name_overview, rpca_clusters.v2, species,severity)]
n_goodcells_pur[, sum_across_species_perCluster := sum(N), .(uniform_name_overview,rpca_clusters.v2)]
n_goodcells_pur[, anteil_celltyp := sum_across_species_perCluster/sum(N), rpca_clusters.v2][order(anteil_celltyp)]
hist(n_goodcells_pur$anteil_celltyp, breaks = 50)

n_goodcells_pur[, uniform_name_overview_keep := anteil_celltyp>0.5]



stats_keep = n_goodcells_pur[, sum(N),uniform_name_overview_keep]
stats_keep


stats_keep_celltype = n_goodcells_pur[, sum(N),.(uniform_name_overview_keep, uniform_name_overview, species, severity)]
stats_keep_celltype
stats_keep_celltype2 = dcast.data.table(stats_keep_celltype,species + severity + uniform_name_overview    ~  uniform_name_overview_keep       , value.var = "V1")
stats_keep_celltype2[, uniform_name_overview_keep_percent := `TRUE`/(`TRUE`+`FALSE`)]

plotly::ggplotly(p1a)
plotly::ggplotly(p2a)


DT::datatable(stats_keep_celltype2)



goodcells[, pasteid := paste(uniform_name_overview , rpca_clusters.v2)]
n_goodcells_pur[, pasteid := paste(uniform_name_overview , rpca_clusters.v2)]

goodcells[, uniform_name_overview_keep := n_goodcells_pur[match_hk(goodcells$pasteid, n_goodcells_pur$pasteid, makeunique = T, importcol = n_goodcells_pur$uniform_name_overview_keep), uniform_name_overview_keep]]

goodcells[,.N, uniform_name_overview_keep]



## update seurat object
if(exists('seurat_filtered')==F) seurat_filtered = readRDS(here("R/results/h225_5_seurat_integrated_rpca_FILTERED.RDS"))
seurat_filtered$severity = ifelse(seurat_filtered$species =="human", seurat_filtered$who_per_sample, seurat_filtered$timepoint)
seurat_filtered$severity %>% table()

seurat_filtered$uniform_name_overview2 = goodcells[match_hk(colnames(seurat_filtered), goodcells$rn),uniform_name_overview2]

seurat_filtered$uniform_name_overview_keep = goodcells[match_hk(colnames(seurat_filtered), goodcells$rn),uniform_name_overview_keep]
seurat_filtered$uniform_name_overview_keep %>% mytable

p3a = DimPlot(seurat_filtered, label = T, repel = T, group.by = "rpca_clusters.v2", raster = F, shuffle = T, reduction = "umap.rpca.v2", split.by = "species")


p4a = DimPlot(seurat_filtered, label = T, repel = T, group.by = "uniform_name_overview2", raster = F, shuffle = T, reduction = "umap.rpca.v2", split.by = "species")

p3a+p4a + plot_annotation(title = "138 698   cells after QC")



goodcells[,.N, .(uniform_name_overview_keep, uniform_name_overview2)][order(uniform_name_overview2,uniform_name_overview_keep)]

seurat_filtered_purecells = seurat_filtered[, seurat_filtered$uniform_name_overview_keep ==T]
seurat_filtered_purecells %>% dim


p3b = DimPlot(seurat_filtered_purecells, label = T, repel = T, group.by = "rpca_clusters.v2", raster = F, shuffle = T, reduction = "umap.rpca.v2", split.by = "species")

p4b = DimPlot(seurat_filtered_purecells, label = T, repel = T, group.by = "uniform_name_overview2", raster = F, shuffle = T, reduction = "umap.rpca.v2", split.by = "species")

p3b+p4b + plot_annotation(title = "7601 5.5% cells with ambiguos annotation excluded")

# plotly::ggplotly( DimPlot(seurat_filtered_purecells[, seurat_filtered_purecells$rpca_clusters.v2==22], label = T,   group.by = "uniform_name_overview2", raster = F, shuffle = T, reduction = "umap.rpca.v2", split.by = "species"))


p5a = ggplot(goodcells[uniform_name_overview_keep==T&grepl("^un", uniform_name_overview)==F], aes(rpca_clusters.v2, fill = factor(uniform_name_overview))) + geom_bar(position = "fill") + facet_grid(species~.)
p6a = ggplot(goodcells[uniform_name_overview_keep==T&grepl("^un", uniform_name_overview)==F], aes(rpca_clusters.v2, fill = factor(uniform_name_overview))) + geom_bar() + facet_grid(species~.)
p5a
p6a

## are the megakaryocytes platelets



VlnPlot(seurat_filtered_purecells[,], features = c("PF4", "PPBP", "CCL5"), group.by = "uniform_name_overview2", split.by = "species")  # https://panglaodb.se/markers/PanglaoDB_markers_27_Mar_2020.tsv.gz


megakaryogenes = c('Vwf', 'Pf4', 'Gata1', 'Selp', 'Gp6',  'Gp1ba') %>% toupper()
VlnPlot(seurat_filtered_purecells[,], features = megakaryogenes, group.by = "uniform_name_overview2", split.by = "species")  # https://www.ahajournals.org/doi/10.1161/ATVBAHA.119.313280

seurat_filtered_purecells$species_celltype = paste(seurat_filtered_purecells$species, seurat_filtered_purecells$uniform_name_overview2)
DimPlot(seurat_filtered_purecells[,seurat_filtered_purecells$uniform_name_overview2 %in% c("Platelet", "Megakaryocyte", "Neutrophils")], group.by = "species_celltype",reduction = "umap.rpca.v2" , shuffle = T)

# p_celltypes = DimPlot(seurat_filtered_purecells, group.by = "uniform_name_overview2",reduction = "umap.rpca.v2" , shuffle = T,  raster = F, label = T, split.by =  "species")
# p_celltypes


unique(seurat_filtered_purecells$uniform_name_overview2)

neutrosubset = seurat_filtered_purecells[,seurat_filtered_purecells$uniform_name_overview2 %in% c("Neutrophils", "Immature neutrophil" )]
neutrosubset



p_celltypes2 = DimPlot(neutrosubset, group.by = c("uniform_name_overview2",'experiment') ,reduction = "umap.rpca.v2" , shuffle = T,  raster = F, label = T, split.by =  "species")
p_celltypes2

p_celltypes3 = DimPlot(neutrosubset, group.by = c("uniform_name_overview2",'rpca_clusters.v2', "cluster_labels_res.0.8") ,reduction = "umap.rpca.v2" , shuffle = T,  raster = F, label = T, split.by =  "species")
p_celltypes3


FeaturePlot(neutrosubset, features = c('CXCR4','CXCR2'), split.by = "species", reduction = "umap.rpca.v2") # ist auch CXCR2 low in https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7442692/ beschrieben fuer immatures , CXCR4 high


# plotly::ggplotly(p_celltypes)

## wie sieht das in den gegebenen Daten aus
human_ori = readRDS(here("R/data/seurat_COVID19_freshWB-PBMC_cohort2_rhapsody_jonas_FG_2020-08-18.rds"))


human_ori_anno = human_ori@meta.data %>% as.data.table(keep.rownames = T)
human_ori_anno[,.N,cluster_labels_res.0.8]

neutrosubset_human = neutrosubset[,neutrosubset$species =="human"]
neutrosubset_human
qlist1 = venn2(human_ori_anno$rn, colnames(neutrosubset_human))

neutrosubset_human$cluster_labels_res.0.8 %>% mytable()

neutro_unintegrated = DimPlot(neutrosubset_human, group.by = c("cluster_labels_res.0.8", "uniform_name_overview2"))


neutro_integrated = DimPlot(neutrosubset_human, group.by = c("cluster_labels_res.0.8", "uniform_name_overview2"),reduction = "umap.rpca.v2"  )

neutro_unintegrated / neutro_integrated

cellanno_neutrosubset_human  =neutrosubset_human@meta.data %>% as.data.table(keep.rownames = T)

cellanno_neutrosubset_human[uniform_name_overview2 %in% "Immature neutrophil",.N,.(severity,cluster_labels_res.0.8 )][order(severity, cluster_labels_res.0.8)]


cellanno_neutrosubset_human[uniform_name_overview2 %in% "Immature neutrophil",.N, .(rpca_clusters.v2, cluster_labels_res.0.8)][order(cluster_labels_res.0.8,N)]

DimPlot(neutrosubset_human, group.by = c("cluster_labels_res.0.8", "rpca_clusters.v2"), reduction = "umap.rpca.v2", label = T)

## ok wir koennen es mit den einsern machen
seurat_filtered_purecells$uniform_name_overview3 = ifelse(seurat_filtered_purecells$rpca_clusters.v2 ==32 &seurat_filtered_purecells$uniform_name_overview2 == "Immature neutrophil", "Immature Neutrophils 2",
                                                          ifelse(seurat_filtered_purecells$uniform_name_overview2== "Immature neutrophil",  "Immature Neutrophils 1",seurat_filtered_purecells$uniform_name_overview2))


p_celltypes0 = DimPlot(seurat_filtered_purecells, group.by = c("uniform_name_overview3"),label = T, repel = T, raster = F, shuffle = T, reduction = "umap.rpca.v2")
p_celltypes0


p_celltypes0

cellanno_filtered_purecells  = seurat_filtered_purecells@meta.data %>% as.data.table(keep.rownames = T)
cellanno_filtered_purecells[,.N, .(uniform_name_overview2,uniform_name_overview3)]


p_species_before = DimPlot(seurat_filtered_purecells, group.by = c("species"), repel = T, raster = F, shuffle = T)

p_species_after = DimPlot(seurat_filtered_purecells, group.by = c("species"), repel = T, raster = F, shuffle = T, reduction = "umap.rpca.v2")


p_species_before + p_species_after + plot_layout(guides = "collect")


## wh plot neutro subset ----
neutrosubset2 = seurat_filtered_purecells[,seurat_filtered_purecells$uniform_name_overview2 %in% c( "Immature neutrophil" )]
neutrosubset2



p_celltypes4 = DimPlot(neutrosubset2, group.by = c("rpca_clusters.v2","uniform_name_overview3") ,reduction = "umap.rpca.v2" , shuffle = T,  raster = F, label = T, split.by =  "species")
p_celltypes4


FeaturePlot(neutrosubset, features = c('ELANE', 'FUT4', "SPN" ,'CEACAM8', "PADI4"), split.by = "species",reduction = "umap.rpca.v2" ) # ELANE FUT4 spn Pro-Neu, PADI4 Pre-Neu in schulte schrepping f cells with a proliferative signature  # Neutrophils from COVID-19, particularly from patients with severe disease, primarily occupied immature pre- and pro-neutrophil-like clusters.  The TF network underlying the transcriptional difference in pro-neutrophils is mainly driven by E2F family members and pre-neutrophils mainly depend on ETS TFs (Figure S6H).
# Pseudotime analysis strongly supported the differentiation trajectory from pro-neutrophils (cluster 8) via pre-neutrophils (cluster 6) to mature neutrophils in cluster 2 and 1 (Figures S6I and S6J). Pseudotime analysis strongly supported the differentiation trajectory from pro-neutrophils (cluster 8) via pre-neutrophils (cluster 6) to mature neutrophils in cluster 2 and 1 (Figures S6I and S6J).


## wh plot CD4  subset ----
cd4subset2 = seurat_filtered_purecells[,seurat_filtered_purecells$uniform_name_overview2 %in% c( "CD4+_T_Cells" )]
cd4subset2

cd4subset2old = seurat_filtered[,seurat_filtered$uniform_name_overview2 %in% c( "CD4+_T_Cells" )]
cd4subset2old



p_celltypes5 = DimPlot(cd4subset2, group.by = c("rpca_clusters.v2","uniform_name_overview3", "severity") ,reduction = "umap.rpca.v2" , shuffle = T,  raster = F, label = T, split.by =  "species")
p_celltypes5

p_celltypes5b = DimPlot(cd4subset2, group.by = c("rpca_clusters.v2","uniform_name_overview3") ,reduction = "umap.rpca.v2" , shuffle = T,  raster = F, label = T, split.by =  "species")
p_celltypes5b


# bei old ist das cluster 15 nicht mehr dabei
p_celltypes5old = DimPlot(cd4subset2old, group.by = c("rpca_clusters.v2","uniform_name_overview2", "severity") ,reduction = "umap.rpca.v2" , shuffle = T,  raster = F, label = T, split.by =  "species")
p_celltypes5old

p_celltypes5bold = DimPlot(cd4subset2old, group.by = c("rpca_clusters.v2","uniform_name_overview2") ,reduction = "umap.rpca.v2" , shuffle = T,  raster = F, label = T, split.by =  "species")
p_celltypes5bold

# das cluster 15 ist ca 10x seltener in hamsters
goodcells[rpca_clusters.v2==15, .N,.(species,uniform_name_overview2)][order(species,-N)]
goodcells[,sum(rpca_clusters.v2==15)/.N, species]

VlnPlot(cd4subset2old[,cd4subset2old$rpca_clusters.v2==15], features = c('RORC', 'NCR2', "IL7R", "KLRF1", "TBX21", "GNLY"), group.by = "species", split.by = "uniform_name_overview2")
# sind das die https://twitter.com/rmassonix/status/1677009892647870464


# speichern ###
jpeg2(here("R/results/h228_3_umap_before_after_integration_species.jpeg"), 11,5)
p_species_before + p_species_after + plot_layout(guides = "collect") + plot_annotation(caption = "h228_3_umap_before_after_integration_species.jpeg")
dev.off()

jpeg2(here("R/results/h228_3_umapcelltypes_allspecies.jpeg"), 9,6)
p_celltypes0 + plot_annotation(caption = "h228_3_umapcelltypes_allspecies.jpeg")
dev.off()

fwrite(cellanno_filtered_purecells, here("R/results/h228_3_qc_integrated_rpca_cellanno_FILTEREDqc_FILTEREDpureCells.txt.gz"), sep = "\t")



saveRDS(seurat_filtered_purecells, here("R/results/h228_3_qc_integrated_rpca_cellanno_FILTEREDqc_FILTEREDpureCells.rds"))
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
##  [1] patchwork_1.1.2         ggplot2_3.4.2           Seurat_4.9.9.9044
##  [4] SeuratObject_4.9.9.9084 sp_1.6-0                here_1.0.1
##  [7] toolboxH_0.2.17         eulerr_7.0.0            testthat_3.1.7
## [10] stringr_1.5.0           scales_1.2.1            readxl_1.4.2
## [13] RColorBrewer_1.1-3      png_0.1-8               fdrtool_1.2.17
## [16] R.utils_2.12.2          R.oo_1.25.0             R.methodsS3_1.8.2
## [19] data.table_1.14.8
##
## loaded via a namespace (and not attached):
##   [1] spam_2.9-1             plyr_1.8.8             igraph_1.4.1
##   [4] lazyeval_0.2.2         splines_4.2.3          RcppHNSW_0.4.1
##   [7] crosstalk_1.2.0        listenv_0.9.0          scattermore_1.0
##  [10] digest_0.6.31          htmltools_0.5.5        fansi_1.0.4
##  [13] magrittr_2.0.3         tensor_1.5             cluster_2.1.4
##  [16] ROCR_1.0-11            globals_0.16.2         matrixStats_0.63.0
##  [19] spatstat.sparse_3.0-1  colorspace_2.1-0       ggrepel_0.9.3
##  [22] xfun_0.38              dplyr_1.1.1            jsonlite_1.8.4
##  [25] progressr_0.13.0       spatstat.data_3.0-1    survival_3.5-0
##  [28] zoo_1.8-11             glue_1.6.2             polyclip_1.10-4
##  [31] gtable_0.3.3           leiden_0.4.3           future.apply_1.10.0
##  [34] abind_1.4-5            spatstat.random_3.1-4  miniUI_0.1.1.1
##  [37] Rcpp_1.0.10            viridisLite_0.4.2      xtable_1.8-4
##  [40] reticulate_1.28        dotCall64_1.0-2        DT_0.27
##  [43] htmlwidgets_1.6.2      httr_1.4.5             ellipsis_0.3.2
##  [46] ica_1.0-3              pkgconfig_2.0.3        farver_2.1.1
##  [49] sass_0.4.5             uwot_0.1.14            deldir_1.0-6
##  [52] utf8_1.2.3             tidyselect_1.2.0       labeling_0.4.2
##  [55] rlang_1.1.0            reshape2_1.4.4         later_1.3.0
##  [58] munsell_0.5.0          cellranger_1.1.0       tools_4.2.3
##  [61] cachem_1.0.7           cli_3.6.1              generics_0.1.3
##  [64] ggridges_0.5.4         evaluate_0.21          fastmap_1.1.1
##  [67] yaml_2.3.7             goftest_1.2-3          knitr_1.42
##  [70] fitdistrplus_1.1-8     purrr_1.0.1            RANN_2.6.1
##  [73] pbapply_1.7-0          future_1.32.0          nlme_3.1-161
##  [76] mime_0.12              ggrastr_1.0.1          brio_1.1.3
##  [79] compiler_4.2.3         rstudioapi_0.14        beeswarm_0.4.0
##  [82] plotly_4.10.1          spatstat.utils_3.0-2   tibble_3.2.1
##  [85] bslib_0.4.2            stringi_1.7.12         RSpectra_0.16-1
##  [88] lattice_0.20-45        Matrix_1.5-3           vctrs_0.6.1
##  [91] pillar_1.9.0           lifecycle_1.0.3        spatstat.geom_3.1-0
##  [94] lmtest_0.9-40          jquerylib_0.1.4        RcppAnnoy_0.0.20
##  [97] cowplot_1.1.1          irlba_2.3.5.1          httpuv_1.6.9
## [100] R6_2.5.1               promises_1.2.0.1       KernSmooth_2.23-20
## [103] gridExtra_2.3          vipor_0.4.5            parallelly_1.35.0
## [106] codetools_0.2-19       fastDummies_1.6.3      MASS_7.3-58.2
## [109] rprojroot_2.0.3        withr_2.5.0            sctransform_0.3.5
## [112] parallel_4.2.3         grid_4.2.3             tidyr_1.3.0
## [115] rmarkdown_2.21         Cairo_1.6-0            Rtsne_0.16
## [118] spatstat.explore_3.1-0 shiny_1.7.4            ggbeeswarm_0.7.1
