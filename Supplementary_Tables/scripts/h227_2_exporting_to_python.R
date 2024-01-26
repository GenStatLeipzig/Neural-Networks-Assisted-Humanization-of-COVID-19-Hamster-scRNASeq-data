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
