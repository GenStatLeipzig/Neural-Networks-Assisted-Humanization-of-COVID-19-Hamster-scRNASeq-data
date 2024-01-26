require(toolboxH)
require(here)
require(Seurat)



goodcells <- fread(here("R/results/h225_4_qc_integrated_rpca_cellanno_FILTERED.txt.gz"))
require(ggplot2)

## update cell names
goodcells[, uniform_name_overview2 := ifelse(celltype == "Immature neutrophil" | cluster_labels_res.0.8 %in% c("Immature Neutrophils_1", "Immature Neutrophils_2"), "Immature neutrophil", uniform_name_overview)]


p1a <- ggplot(goodcells[grepl("^un", uniform_name_overview) == F], aes(rpca_clusters.v2, fill = factor(uniform_name_overview))) +
  geom_bar(position = "fill") +
  facet_grid(species ~ .)
p2a <- ggplot(goodcells[grepl("^un", uniform_name_overview) == F], aes(rpca_clusters.v2, fill = factor(uniform_name_overview))) +
  geom_bar() +
  facet_grid(species ~ .)


goodcells[, severity := ifelse(species == "human", who_per_sample, timepoint)]

n_goodcells_pur <- goodcells[, .N, .(uniform_name_overview, rpca_clusters.v2, species, severity)]
n_goodcells_pur[, sum_across_species_perCluster := sum(N), .(uniform_name_overview, rpca_clusters.v2)]
n_goodcells_pur[, anteil_celltyp := sum_across_species_perCluster / sum(N), rpca_clusters.v2][order(anteil_celltyp)]

n_goodcells_pur[, uniform_name_overview_keep := anteil_celltyp > 0.5]


goodcells[, pasteid := paste(uniform_name_overview, rpca_clusters.v2)]
n_goodcells_pur[, pasteid := paste(uniform_name_overview, rpca_clusters.v2)]

goodcells[, uniform_name_overview_keep := n_goodcells_pur[match_hk(goodcells$pasteid, n_goodcells_pur$pasteid, makeunique = T, importcol = n_goodcells_pur$uniform_name_overview_keep), uniform_name_overview_keep]]


fwrite(goodcells, here("R/results/h228_1_qc_integrated_rpca_cellanno_FILTEREDqc_FILTEREDpureCells.txt.gz"), sep = "\t")



## finalize script ----
finalizeSkript()
