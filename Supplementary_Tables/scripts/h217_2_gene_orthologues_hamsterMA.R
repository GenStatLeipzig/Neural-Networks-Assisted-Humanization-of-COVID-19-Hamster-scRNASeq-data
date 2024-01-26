
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
