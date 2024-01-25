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
