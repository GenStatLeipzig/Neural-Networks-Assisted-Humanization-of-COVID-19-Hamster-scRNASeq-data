rm(list = ls())


.libPaths("/net/ifs1/san_projekte/projekte/genstat/07_programme/rpackages/angmar/")
.libPaths()
packageVersion("Seurat")
require(Seurat)
packageVersion("SeuratObject")
require(SeuratObject)
require(toolboxH)
require(here)
require(stringr)

# # load anno and create flat file ----
anno = fread(here("R/data/phorob_curated_rev3_S2.gtf"))
anno[,.N, .(V2,V3)]
anno[,.N, .(V1)][order(N)][plot(N)]
anno[,.N, .(V6, V7, V8)]
anno[grep("rtrack", V2),unique(V1)] # this is annotation for the Sars-CoV2 virus

anno2 = anno[(grepl("rtrack", V2) & V3=="exon" )==F] # keeping CDS, not exon
anno2[,.N, .(V2, V3)] # now filter on transcripts, only

anno3 = anno2[grepl("transcr", V3)|grepl("rtrack", V2)]

anno3[, transcript_id := str_split(V9, ";") %>% sapply(., function(x) grep("transcript_id", x, value = T)) %>% str_trim()]
anno3[,.N, str_split(transcript_id, " ") %>% sapply("[", 1)]
anno3[, transcript_id:= str_replace_all(transcript_id, 'transcript_id |\\"', "")]

anno3[, gene_id := str_split(V9, ";") %>% sapply(., function(x) grep("gene_id", x, value = T)) %>% str_trim()]
anno3[,.N, str_split(gene_id, " ") %>% sapply("[", 1)]
anno3[, gene_id:= str_replace_all(gene_id, 'gene_id |\\"', "")]


anno3[, refgene := str_split(V9, ";") %>% sapply(., function(x) grep("ref-gene", x, value = T)) %>% str_trim()]
anno3[,.N, str_split(refgene, " ") %>% sapply("[", 1)]
anno3[refgene == "character(0)", refgene := NA]
anno3[, refgene := str_replace_all(refgene, 'ref-gene |\\"', "")]

anno3[, score.2 := str_split(V9, ";") %>% sapply(., function(x) grep("score", x, value = T)) %>% str_trim()]
anno3[,.N, str_split(score.2, " ") %>% sapply("[", 1)]
anno3[score.2 == "character(0)", score.2 := NA]
anno3[,.N, score.2]
anno3[, score.2 := str_replace_all(score.2, 'score |\\"', "") %>% as.numeric()]


anno3[, sumWeight := str_split(V9, ";") %>% sapply(., function(x) grep("sumWeight", x, value = T)) %>% str_trim()]
anno3[,.N, str_split(sumWeight, " ") %>% sapply("[", 1)]
anno3[sumWeight == "character(0)", sumWeight := NA]
anno3[, .N,sumWeight]
anno3[, sumWeight := str_replace_all(sumWeight, 'sumWeight |\\"', "") %>% as.numeric()]

anno3[, evidence := str_split(V9, ";") %>% sapply(., function(x) grep("evidence", x, value = T)) %>% str_trim()]
anno3[,.N, str_split(evidence, " ") %>% sapply("[", 1)]
anno3[evidence == "character(0)", evidence := NA]
anno3[, evidence := str_replace_all(evidence, 'evidence |\\"', "") %>% as.numeric()]
anno3[,.N, .(sumWeight,evidence)] # ok, identisch

anno3[, alternative := str_split(V9, ";") %>% sapply(., function(x) grep("alternative", x, value = T)) %>% str_trim()]
anno3[,.N, str_split(alternative, " ") %>% sapply("[", 1)]
anno3[alternative == "character(0)", alternative := NA]
anno3[, alternative := str_replace_all(alternative, 'alternative |\\"', "") ]



anno3[, transcripts := str_split(V9, ";") %>% sapply(., function(x) grep("transcripts", x, value = T)) %>% str_trim()]
anno3[,.N, str_split(transcripts, " ") %>% sapply("[", 1)]
anno3[,.N, transcripts]
anno3[transcripts == "character(0)", transcripts := NA]
anno3[, transcripts := str_replace_all(transcripts, 'transcripts |\\"', "") %>% as.numeric()]
anno3[transcripts==3][gene_id==gene_id[1]] # todo anguggen, was ein gene mehrere transcripts bedeutet

anno3[, complete := str_split(V9, ";") %>% sapply(., function(x) grep("complete", x, value = T)) %>% str_trim()]
anno3[,.N, str_split(complete, " ") %>% sapply("[", 1)]
anno3[complete == "character(0)", complete := NA]
anno3[, complete := str_replace_all(complete, 'complete |\\"', "") %>% as.numeric()]
anno3[, .N,complete]

anno3[, maxEvidence := str_split(V9, ";") %>% sapply(., function(x) grep("maxEvidence", x, value = T)) %>% str_trim()]
anno3[,.N, str_split(maxEvidence, " ") %>% sapply("[", 1)]
anno3[maxEvidence == "character(0)", maxEvidence := NA]
anno3[, maxEvidence := str_replace_all(maxEvidence, 'maxEvidence |\\"', "") %>% as.numeric()]
anno3[, .N,.(evidence, maxEvidence)]

anno3[, combinedEvidence := str_split(V9, ";") %>% sapply(., function(x) grep("combinedEvidence", x, value = T)) %>% str_trim()]
anno3[,.N, str_split(combinedEvidence, " ") %>% sapply("[", 1)]
anno3[combinedEvidence == "character(0)", combinedEvidence := NA]
anno3[, combinedEvidence := str_replace_all(combinedEvidence, 'combinedEvidence |\\"', "") %>% as.numeric()]
anno3[, .N,.(evidence,maxEvidence,combinedEvidence)]

anno3[, gene_name := str_split(V9, ";") %>% sapply(., function(x) grep("gene_name", x, value = T)) %>% str_trim()]
anno3[,.N, str_split(gene_name, " ") %>% sapply("[", 1)]
anno3[, gene_name := str_replace_all(gene_name, 'gene_name |\\"', "") ]

check2 = unique(anno3[,.(gene_id, gene_name)])
check2[allDuplicatedEntries(gene_id)]
check2[allDuplicatedEntries(gene_name)]#sometimes same gene_name for different gene_id, but rare

setnames(anno3, paste0("V", 1:9), c("seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"))
anno3# from

# https://en.wikipedia.org/wiki/General_feature_format

# 1	seqid	The name of the sequence where the feature is located.
# 2	source	Keyword identifying the source of the feature, like a program (e.g. Augustus or RepeatMasker) or an organization (like TAIR).
# 3	type	The feature type name, like "gene" or "exon". In a well structured GFF file, all the children features always follow their parents in a single block (so all exons of a transcript are put after their parent "transcript" feature line and before any other parent transcript line). In GFF3, all features and their relationships should be compatible with the standards released by the Sequence Ontology Project.
# 4	start	Genomic start of the feature, with a 1-base offset. This is in contrast with other 0-offset half-open sequence formats, like BED.
# 5	end	Genomic end of the feature, with a 1-base offset. This is the same end coordinate as it is in 0-offset half-open sequence formats, like BED.[citation needed]
# 6	score	Numeric value that generally indicates the confidence of the source in the annotated feature. A value of "." (a dot) is used to define a null value.
# 7	strand	Single character that indicates the strand of the feature; it can assume the values of "+" (positive, or 5'->3'), "-", (negative, or 3'->5'), "." (undetermined).
# 8	phase	phase of CDS features; it can be either one of 0, 1, 2 (for CDS features) or "." (for everything else). See the section below for a detailed explanation.
# 9	attributes	All the other information pertaining to this feature. The format, structure and content of this field is the one which varies the most between the three competing file formats.


# Match gene annotatio with gene names from seurat object ----


robo = readRDS(here("R/data/pr_blood_new_combined_integrated_annotated.rds"))
robo


qlist1 = venn2(rownames(robo), anno3$gene_name)
str(qlist1)
# Some genes require renaming to fit provided annotation

robo_renaming = data.table(seurat_name = rownames(robo),
                           anno_name = rownames(robo))
robo_renaming

robo_renaming[,  anno_name := str_replace(anno_name, "unknown-gene-", "unknown_gene_")]  # EMANUELFRAGE
robo_renaming[,  anno_name := str_replace(anno_name, "NEWGENE-", "NEWGENE_")]  # EMANUELFRAGE
robo_renaming[,  anno_name := str_replace(anno_name, "SCoV2-orf1ab-", "SCoV2_orf1ab_")]  # EMANUELFRAGE
robo_renaming[,  anno_name := str_replace(anno_name, "SCoV2-", "SCoV2_")]  # EMANUELFRAGE

qlist2 = venn2(robo_renaming$anno_name, anno3$gene_name)
qlist2$q2 # these seem to be duplicated

# duplicate entries in annotation for duplicated gene names
toduplicate = data.table(dupliname = c("Atxn7l1.1", "Utrn.1", "Plekhn1.1", "Asph.1"),
                         oriname = c("Atxn7l1", "Utrn", "Plekhn1", "Asph"))

robo_renamingb = robo_renaming[seurat_name %in% toduplicate$oriname]
robo_renamingb[,seurat_name := paste0(seurat_name, ".1")]

robo_renamingb


# add seurat name to annotation
anno3[,seurat_name := robo_renaming[match_hk(anno3$gene_name, robo_renaming$anno_name), seurat_name]]

anno3adder = anno3[gene_name %in% robo_renamingb$anno_name]
anno3adder[, seurat_name:= robo_renamingb[match_hk(anno3adder$gene_name, robo_renamingb$anno_name),seurat_name]]
robo_renaming2 = rbind(robo_renaming,
                       robo_renamingb) %>% unique()


anno4 = rbind(anno3, anno3adder) # EMANUELFRAGE

qlist4 = venn2(rownames(robo), anno4$seurat_name)
stopifnot("not all genes from seurat object annotated" = sum(length(c(qlist4$q2, qlist4$q3)))==0)

# split more annotation in subcategories
anno4[, refgene_species := str_split(refgene, "\\_") %>% sapply("[", 1)]
anno4[, .N,refgene_species]

anno4[is.na(refgene_species), .(gene_name, seurat_name, refgene)] %>% data.frame() # ok, all are SARS-CoV-2 genes
anno4[refgene_species=="manual"]  # EMANUELFRAGE Is ths simply https://www.genecards.org/cgi-bin/carddisp.pl?gene=NR3C1#orthologs


anno4[, refgene_geneid := str_split(refgene, "\\_") %>% sapply("[", 2) %>% str_replace("gene\\:", "")%>% str_replace("\\.gene$", "")]
anno4[, .N,str_sub(refgene_geneid, 1,7)]
anno4[str_sub(refgene_geneid, 1,7)=="transcr"]  # EMANUELFRAGE about transcript and why mesaur_ENSMAUT00000000055.gene not ENSMAUG*****

# # LOAD orthogolgues v2 ----
# Ensembl release 109 - Feb 2023 © EMBL-EBI
# http://Feb2023.archive.ensembl.org/biomart/martview/f7671ccf9ae60f59cd8c7520c36c370d
hamster = fread(here("R/data/martquery_0522114054_418_human_MesAurOrthoV2.txt.gz"))
rat = fread(here("R/data/martquery_RAT_V2.txt.gz"))
mouse = fread(here("R/data/martquery_MOUSE_V2.txt.gz"))

qlsit5334 = venn3(hamster$`Gene name`,
                  rat$`Gene name`,
                  mouse$`Gene name`)
orthologues_pre = merge(hamster, rat[,c( 'Gene stable ID',
                                      "Rat gene stable ID",
                                      "Rat gene name",
                                      "Rat orthology confidence [0 low, 1 high]" ), with = F]
                                 , by = 'Gene stable ID', all = T)

orthologues = merge(orthologues_pre, mouse[,c( 'Gene stable ID',
                                     "Mouse gene stable ID",
                                     "Mouse gene name",
                                     "Mouse orthology confidence [0 low, 1 high]" ), with = F]
                    , by = 'Gene stable ID', all = T, allow.cartesian=TRUE) %>% unique()

orthologuesm = rbind(hamster[,.(human_id=`Gene stable ID`,
                                human_idversion = `Gene stable ID version`,
                                human_name = `Gene name`,
                                description =`Gene description`,
                                orthologue_geneid = `Golden Hamster gene stable ID`,
                                orthologue_genename= `Golden Hamster gene name`,
                                orthologue_confidence = `Golden Hamster orthology confidence [0 low, 1 high]`)][,species := "mesaur"],

                     rat[,.(human_id=`Gene stable ID`,
                            human_idversion = `Gene stable ID version`,
                            human_name = `Gene name`,
                            description =`Gene description`,
                            orthologue_geneid = `Rat gene stable ID`,
                                orthologue_genename= `Rat gene name`,
                                orthologue_confidence = `Rat orthology confidence [0 low, 1 high]`)][,species := "rat"],

                     mouse[,.(human_id=`Gene stable ID`,
                              human_idversion = `Gene stable ID version`,
                              human_name = `Gene name`,
                              description =`Gene description`,
                              orthologue_geneid = `Mouse gene stable ID`,
                            orthologue_genename= `Mouse gene name`,
                            orthologue_confidence = `Mouse orthology confidence [0 low, 1 high]`)][,species := "mus"]
                     )


qlist1231 = venn2(orthologuesm$orthologue_genename %>% toupper(), anno4$gene_name %>% toupper())

qlist1231b = venn3(orthologuesm$orthologue_genename %>% toupper(), anno4$gene_name %>% toupper(), orthologuesm$`Gene name` %>% toupper())

str(qlist1231)

### considering source of gene
qlist1231c = venn2(orthologuesm[, paste(orthologue_genename, species)]  %>% toupper(), anno4[, paste(gene_name, refgene_species)] %>% toupper())

qlist1231d = venn2(orthologuesm[, paste(orthologue_geneid)]  %>% toupper(), anno4[, paste(refgene_geneid)] %>% toupper())



#  # replace "" with NA ----

for(i in names(orthologues)) {
  orthologues[get(i)=="", (i) := NA]

}


for(i in names(orthologuesm)) {
  orthologuesm[get(i)=="", (i) := NA]

}

orthologuesm$mergeid1 = orthologuesm$orthologue_geneid
anno4$mergeid1 = anno4$refgene_geneid

orthologuesm2 = merge(orthologuesm %>% unique(), (anno4[is.na(mergeid1)==F,.(mergeid1,seurat_name1= seurat_name)] %>% unique()), by ="mergeid1", all.x = T,allow.cartesian=TRUE) %>% unique()

orthologuesm2[is.na(seurat_name1)==F, orthologue_by := '01:gene ID']

## nun merge nach namen der gene je nach zugeordneten species
orthologuesm2[is.na(orthologue_genename)==F, mergeid2 := paste(orthologue_genename, species)]
anno4[is.na(gene_name)==F, mergeid2 := paste(gene_name, refgene_species)]

orthologuesm3 = merge(orthologuesm2 %>% unique(), (anno4[is.na(mergeid2)==F,.(mergeid2, seurat_name2 = seurat_name)] %>% unique()), by = "mergeid2" , all.x = T,allow.cartesian=TRUE) %>% unique()


orthologuesm3[is.na(seurat_name2)==F & is.na(seurat_name1)==T, orthologue_by := '02:gene name of species']

orthologuesm3[, seurat_name := ifelse(is.na(seurat_name1), seurat_name2, seurat_name1)]

orthologuesm3[, uniqueN(seurat_name, na.rm = T), orthologue_by]



## nun merge nach namen der gene je nach zugeordneten species
orthologuesm3[is.na(human_name)==F, mergeid3 := toupper(human_name)]
anno4[is.na(gene_name)==F, mergeid3 := toupper(gene_name)]

orthologuesm4 = merge(orthologuesm3 %>% unique(), (anno4[is.na(mergeid3)==F,.(mergeid3, seurat_name3 = seurat_name)] %>% unique()), by = "mergeid3", all.x = T,allow.cartesian=TRUE) %>% unique()


orthologuesm4[is.na(seurat_name3)==F & is.na(seurat_name)==T, orthologue_by := '03:same name uppercase in human']

orthologuesm4[, seurat_name := ifelse(is.na(seurat_name), seurat_name3, seurat_name)]

orthologuesm4[, uniqueN(seurat_name, na.rm = T), orthologue_by]


qlist68445=venn2(anno4$seurat_name, orthologuesm4$seurat_name)



# # make unique assignment preferring cluster-specific higher expressed genes ----

 ## max expressed pr
 pr_maxexpressed = AverageExpression(robo, assays = "SCT")$SCT %>% apply(.,1, max) %>% as.data.table(keep.rownames = T) %>% setnames(., c("seurat_name", "max_per_cluster"))
 pr_maxexpressed

 orthologuesm4[, rb_max_per_cluster := pr_maxexpressed[match_hk(orthologuesm4$seurat_name, pr_maxexpressed$seurat_name),max_per_cluster]]


 ## max expressed human
 humanSchulteRhap_pre = readRDS(here("R/data/seurat_COVID19_freshWB-PBMC_cohort2_rhapsody_jonas_FG_2020-08-18.rds"))
 humanSchulteRhap_pre

 humanSchulteRhap_pre$cells %>% mytable()
 human = humanSchulteRhap_pre[, humanSchulteRhap_pre$cells == "Whole_blood"]
 human


DimPlot(human, raster = FALSE)

human
 table(human@active.ident)
 avex_human = AverageExpression(human, assays = "RNA", slot = "data")
 human_maxexpressed = avex_human$RNA %>% apply(.,1, max) %>% as.data.table(keep.rownames = T) %>% setnames(., c("seurat_name", "max_per_cluster"))
 human_maxexpressed

 qlist8 = venn2(human_maxexpressed$seurat_name, orthologuesm4$human_name)

  orthologuesm4[, human_max_per_cluster := human_maxexpressed[match_hk(orthologuesm4$human_name, human_maxexpressed$seurat_name),max_per_cluster]]

  orthologuesm[, human_max_per_cluster := human_maxexpressed[match_hk(orthologuesm$human_name, human_maxexpressed$seurat_name),max_per_cluster]]


orthologuesm4[, species := factor(species, levels = c("mesaur",
                                                       "rat",
                                                       "mus"))]
setorder(orthologuesm4,
         orthologue_by,
         species,

         -orthologue_confidence,

         -rb_max_per_cluster,
         human_name,
         na.last = T
         )

orthologuesm4[is.na(human_name)==F &is.na(seurat_name)==F][allDuplicatedEntries(seurat_name)]  # multiple hamster genes for same human gene

orthologues_pr_pre = orthologuesm4[duplicated(seurat_name)==F ]
orthologues_pr_pre[, uniqueN(seurat_name), orthologue_by]
orthologues_pr_pre[is.na(human_name)==F][allDuplicatedEntries(human_name)]

qlist10 = venn2(orthologues_pr_pre$seurat_name, rownames(robo))
qlist11 = venn2(orthologues_pr_pre$human_name, rownames(human))

setorder(orthologues_pr_pre,
         orthologue_by,
         species,

         -orthologue_confidence,

         -human_max_per_cluster,
         seurat_name,
         na.last = T
)

orthologues_pr_pre[is.na(human_name)==F & is.na(seurat_name )==F][allDuplicatedEntries(human_name)] # multiple human genes for same pr gene

orthologues_pr = orthologues_pr_pre[duplicated(human_name)==F ]
orthologues_pr[, uniqueN(seurat_name), orthologue_by]
orthologues_pr[is.na(human_name)==F][allDuplicatedEntries(human_name)]

qlist10b = venn2(orthologues_pr$seurat_name, rownames(robo))
qlist11b = venn2(orthologues_pr$human_name, rownames(human))

anno4[, human_orthologue_name := orthologues_pr[match_hk(anno4$seurat_name, orthologues_pr$seurat_name),human_name]]

anno4[, human_orthologue_description := orthologues_pr[match_hk(anno4$seurat_name, orthologues_pr$seurat_name),description]]

anno4[, human_orthologue_id := orthologues_pr[match_hk(anno4$seurat_name, orthologues_pr$seurat_name),human_id]]

anno4[, human_orthologue_confidence := orthologues_pr[match_hk(anno4$seurat_name, orthologues_pr$seurat_name),orthologue_confidence]]

anno4[, human_orthologue_by := orthologues_pr[match_hk(anno4$seurat_name, orthologues_pr$seurat_name),orthologue_by ]]

anno4[, .N,human_orthologue_by][order(human_orthologue_by)]


# # table for manual check anno -----
duplicatefailure = orthologues_pr_pre[seurat_name %nin% orthologues_pr$seurat_name]
duplicatefailure[, seurat_name_used := anno4[match_hk(duplicatefailure$human_name, anno4$human_orthologue_name,makeunique = T, importcol = anno4$seurat_name), seurat_name]]
duplicatefailure[, orthologue_by_used := anno4[match_hk(duplicatefailure$human_name, anno4$human_orthologue_name,makeunique = T, importcol = anno4$seurat_name), human_orthologue_by]]
duplicatefailure2 = duplicatefailure[, .(human_id, human_name, seurat_name, orthologue_by, seurat_name_used, orthologue_by_used)]

no_orthologue = anno4[is.na(human_orthologue_name) & seurat_name %nin% duplicatefailure$seurat_name & seurat_name %in% robo_renaming2$seurat_name]

no_orthologue2 = no_orthologue[, .(seurat_name, refgene_species, refgene_geneid, human_orthologue_name)][order(seurat_name)]
no_orthologue2
duplicatefailure2


# # save ----
WriteXLS_hk(c("no_orthologue2", "duplicatefailure2"), here("R/results/h216_2_failed_orthologues_Phodopus_roborovskii.xlsx"), AdjWidth = T)

fwrite(anno4, here("R/results/h216_2_gene_annotation_Phodopus_roborovskii.txt.gz"), sep = "\t")

fwrite(orthologues_pr, here("R/results/h216_2_helperfile_orthologues_Phodopus_roborovskii.txt.gz"), sep = "\t")

fwrite(orthologuesm4, here("R/results/h216_2_nXm_orthologues_Phodopus_roborovskii.txt.gz"), sep = "\t")

orthologuesm$mergeid1=NULL
fwrite(orthologuesm, here("R/results/h216_2_nXm_orthologues_mus_rat_mesaur.txt.gz"), sep = "\t")
# # finalize ----
finalizeSkript()
