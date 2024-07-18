#Figure1 B-G and Supplemental Figure 2

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
library(ggExtra)

#Load Data
Blood_merged <- readRDS(".../h228_3_qc_integrated_rpca_cellanno_FILTEREDqc_FILTEREDpureCells.RDS")

#1B AND 1C
Reductions(Blood_merged)

unintegrated_plot = DimPlot(Blood_merged, group.by = 'species',  raster=FALSE)
unintegrated_plot
unintegrated = unintegrated_plot$data %>% as.data.table(keep.rownames = T)
unintegrated[, species2 := ifelse(species =="human", "H. Sapiens",
                                  ifelse(species =="hamsterMA", "M. auratus",
                                         ifelse(species =="hamsterPR", "P. roborovskii", species)))]

integrated_plot = DimPlot(Blood_merged, group.by = 'species', shuffle =T,  raster=FALSE, reduction = "umap.rpca.v2" )
integrated_plot

integrated = integrated_plot$data %>% as.data.table(keep.rownames = T)

integrated[, species2 := ifelse(species =="human", "H. Sapiens",
                                ifelse(species =="hamsterMA", "M. auratus",
                                       ifelse(species =="hamsterPR", "P. roborovskii", species)))]

nogcolors = c("#368F8B", "#F97E44","#B7245C")

## plot unintegrated----

unintegrated_ggplot <- ggplot(unintegrated, aes(umapunintegrated_1  , umapunintegrated_2      , colour = species2)) +
  geom_point(alpha = 0.3) + theme_classic(base_size = 16) + theme(legend.position = "bottom") + xlab('UMAP 1') + ylab('UMAP 2') +   scale_color_manual(values = nogcolors) + scale_fill_manual(values = nogcolors) +  guides(colour = guide_legend(override.aes = list(size=10, alpha = 1))) + labs(color ="", fill = "")



unintegrated_ggplot2 =ggMarginal(unintegrated_ggplot,type = 'density', groupColour = TRUE, groupFill = TRUE, size = 8, margins = 'both')
unintegrated_ggplot2

## plot integrated----

integrated_ggplot <- ggplot(integrated, aes(umaprpcav2_1  , umaprpcav2_2         , colour = species2)) +
  geom_point(alpha = 0.3) + theme_classic(base_size = 16) + theme(legend.position = "bottom") + xlab('UMAP 1') + ylab('UMAP 2') +   scale_color_manual(values = nogcolors) + scale_fill_manual(values = nogcolors) +  guides(colour = guide_legend(override.aes = list(size=10, alpha = 1))) + labs(color ="", fill = "")

integrated_ggplot2 =ggMarginal(integrated_ggplot,type = 'density', groupColour = TRUE, groupFill = TRUE, size = 8, margins = 'both')
integrated_ggplot2


#1D
Blood_merged@active.ident = factor(Blood_merged$uniform_name_overview3)
Cellorder = c("Classical_Monocytes", "Non_Classical_Monocytes","DC", "Immature Neutrophils 1", "Immature Neutrophils 2", "Neutrophils", "NK_Cells", "CD4+_T_Cells", "CD8+_T_Cells", "B_Cells", "Plasmablasts", "Platelet", "Megakaryocyte",  "Proliferating_Cells")
Nogpalette <-  c("#A85C85","#542E42", "#4F6D7A", "#00ABE7", "#6BE2FF", "#0081AF", "#246A73", "#5CC1BC", "#368F8B", "#62C370", "#7C6A0A", "#F7C548", "#F97E44", "#FB3640")
my_levels <- c(Cellorder)
Blood_merged@active.ident <- factor(x = Blood_merged@active.ident, levels = my_levels)
DimPlot(Blood_merged, reduction = "umap.rpca.v2", label = T, shuffle =T , cols = c(Nogpalette)) + NoLegend()


#1E
T8_sub <- subset(Blood_merged, idents = "CD8+_T_Cells")
T8_sub$severity2 <- T8_sub$severity
val_repl <- c("d0", "D0", "d2", "d3", "d5", "e14", "hd_D2", "hd_D3","ld_D2", "ld_D3")     
T8_sub$severity2 <- sapply(T8_sub$severity2,  # Replace values in certain columns
                           function(x) replace(x, x %in% val_repl, 'Hamster'))
DimPlot(T8_sub, group.by = c("severity2"),label = F, repel = T, raster = F, shuffle = T, cols = c("#1a71b8", "#078d84", "#f7b817", "#e27d26", "#c63735", "darkgray"), reduction = "umap.rpca.v2")


#1F
FeaturePlot(T8_sub, features = "CX3CR1", reduction = "umap.rpca.v2") + coord_fixed(ratio=1)
FeaturePlot(T8_sub, features = "GZMA", reduction = "umap.rpca.v2") + coord_fixed(ratio=1)
FeaturePlot(T8_sub, features = "NKG7", reduction = "umap.rpca.v2") + coord_fixed(ratio=1)
FeaturePlot(T8_sub, features = "PRF1", reduction = "umap.rpca.v2") + coord_fixed(ratio=1)
FeaturePlot(T8_sub, features = "EOMES", reduction = "umap.rpca.v2") + coord_fixed(ratio=1)
FeaturePlot(T8_sub, features = "ZEB2", reduction = "umap.rpca.v2") + coord_fixed(ratio=1)


#1G
# Dataimport and wrangling
cellanno_pre = fread("Data/h228_3_qc_integrated_rpca_cellanno_FILTEREDqc_FILTEREDpureCells.txt.gz")
cellanno_pre[,.N,uniform_name_overview_keep]
cellanno = cellanno_pre[uniform_name_overview_keep==T]
cellanno[,.N, species]
cellanno[,.N, uniform_name_overview3]

## counts 
celltypes_all = cellanno[,  unique(uniform_name_overview3)]
celltypes_all

countis = cellanno[,.(N=.N), .(celltype = uniform_name_overview3, species, severity)]
countis[, proportion := N/sum(N), .(species, severity)]

setorder(countis, -proportion)
celltypelevels = unique(countis$celltype) %>% as.character()
countis[, celltype  := factor(celltype, levels =celltypelevels)]

countis$severity %>% unique() %>% dput()

countis[, celltype  := factor(celltype, levels = unique(celltype))]
countis[, celltype2  := str_replace_all(celltype, "_" ," ") %>% str_wrap(width = 15)]

celltype2levels = unique(countis$celltype2)%>% as.character()
countis[, celltype2  := factor(celltype2, levels =celltype2levels)]

severitylevels  = c("0", "3", "4", "5", "7", "d0", "D0", "d2", "d3", "d5", "e14", "ld_D2","ld_D3","hd_D2",   "hd_D3")
countis[, severity  := factor(severity, levels =severitylevels)]

## ADD INDIVIDUAL POINT DATAFRAME   
cellanno[, individual := paste(severity, donor, hamster)] 

cellanno[,.N, .(individual , donor, hamster, severity)]
countis_ind = cellanno[,.(N=.N), .(celltype = uniform_name_overview3, species, severity,individual)]
countis_ind[, proportion := N/sum(N), .(species, severity,individual)]

countis_ind[, celltype  := factor(celltype, levels = celltypelevels  )]
countis_ind[, celltype2  := str_replace_all(celltype, "_" ," ") %>% str_wrap(width = 15)]
countis_ind[, celltype2  := factor(celltype2, levels = celltype2levels)]
countis_ind[, severity  := factor(severity, levels =severitylevels)]

countis_ind$severity2 <- NA
countis_ind$severity2[countis_ind$severity =="0"] <- "Hcontrol"
countis_ind$severity2[countis_ind$severity =="3"] <- "WOS 3"
countis_ind$severity2[countis_ind$severity =="4"] <- "WOS 4"
countis_ind$severity2[countis_ind$severity =="5"] <- "WOS 5"
countis_ind$severity2[countis_ind$severity =="7"] <- "WOS 7"

countis_ind$severity2[countis_ind$severity =="d0"] <- "0MAcontrol"
countis_ind$severity2[countis_ind$severity =="d2"] <- "2 dpi"
countis_ind$severity2[countis_ind$severity =="d3"] <- "3 dpi"
countis_ind$severity2[countis_ind$severity =="d5"] <- "5 dpi"
countis_ind$severity2[countis_ind$severity =="e14"] <- "x14 dpi"

countis_ind$severity2[countis_ind$severity =="D0"] <- "0PRcontrol"
countis_ind$severity2[countis_ind$severity =="ld_D2"] <- "a2 dpi ld"
countis_ind$severity2[countis_ind$severity =="ld_D3"] <- "b3 dpi ld"
countis_ind$severity2[countis_ind$severity =="hd_D2"] <- "c2 dpi hd"
countis_ind$severity2[countis_ind$severity =="hd_D3"] <- "d3 dpi hd"


ggplot(countis_ind[celltype =="Neutrophils" ], aes(severity2, proportion, fill= celltype2)) +
  stat_summary(fill= "#0081AF", geom = "bar", fun = "mean", na.rm = TRUE,
               position = position_dodge(width = 0.4)) +
  geom_errorbar(stat = "summary", fun.data = "mean_sdl", 
                fun.args = list(mult = 1),
                alpha = 0.9,
                position =  position_dodge(width = 0.9), color = "black", linewidth = 0.7, width = 0.4) +
  geom_point( alpha = 0.7, col = "black", size = 3)  + facet_wrap(~species, scales = "free_x" , nrow = 1) + scale_y_continuous(labels = label_percent(), expand = c(0, 0), breaks = pretty_breaks(8)) + guides(fill = "none") + ylab("Neutrophils") + scale_x_discrete(guide = guide_axis(angle = -45)) + theme_classic(base_size = 14) + coord_cartesian(ylim = c(0, 0.89))

ggplot(countis_ind[celltype =="Immature Neutrophils 1" ], aes(severity2, proportion, fill= celltype2)) +
  stat_summary(fill= "#00ABE7", geom = "bar", fun = "mean", na.rm = TRUE,
               position = position_dodge(width = 0.4)) +
  geom_errorbar(stat = "summary", fun.data = "mean_sdl", 
                fun.args = list(mult = 1),
                alpha = 0.9,
                position =  position_dodge(width = 0.9), color = "black", linewidth = 0.7, width = 0.4) +
  geom_point( alpha = 0.7, col = "black", size = 3)  + facet_wrap(~species, scales = "free_x" , nrow = 1) + scale_y_continuous(labels = label_percent(), expand = c(0, 0), breaks = pretty_breaks(8)) + guides(fill = "none") + ylab("Immature Neutrophils 1") + scale_x_discrete(guide = guide_axis(angle = -45)) + theme_classic(base_size = 14) + coord_cartesian(ylim = c(0, 0.30))

ggplot(countis_ind[celltype =="CD4+_T_Cells" ], aes(severity2, proportion, fill= celltype2)) +
  stat_summary(fill= "#5CC1BC", geom = "bar", fun = "mean", na.rm = TRUE,
               position = position_dodge(width = 0.4)) +
  geom_errorbar(stat = "summary", fun.data = "mean_sdl", 
                fun.args = list(mult = 1),
                alpha = 0.9,
                position =  position_dodge(width = 0.9), color = "black", linewidth = 0.7, width = 0.4) +
  geom_point( alpha = 0.7, col = "black", size = 3)  + facet_wrap(~species, scales = "free_x" , nrow = 1) + scale_y_continuous(labels = label_percent(), expand = c(0, 0), breaks = pretty_breaks(8)) + guides(fill = "none") + ylab("CD4+_T_Cells") + scale_x_discrete(guide = guide_axis(angle = -45)) + theme_classic(base_size = 14) + coord_cartesian(ylim = c(0, 0.30))


#S2

ggplot(countis_ind[celltype =="Classical_Monocytes" ], aes(severity2, proportion, fill= celltype2)) +
  stat_summary(fill= "#A85C85", geom = "bar", fun = "mean", na.rm = TRUE,
               position = position_dodge(width = 0.4)) +
  geom_errorbar(stat = "summary", fun.data = "mean_sdl", 
                fun.args = list(mult = 1),
                alpha = 0.9,
                position =  position_dodge(width = 0.9), color = "black", linewidth = 0.7, width = 0.4) +
  geom_point( alpha = 0.7, col = "black", size = 3)  + facet_wrap(~species, scales = "free_x" , nrow = 1) + scale_y_continuous(labels = label_percent(), expand = c(0, 0), breaks = pretty_breaks(8)) + guides(fill = "none") + ylab("Classical_Monocytes") + scale_x_discrete(guide = guide_axis(angle = -45)) + theme_classic(base_size = 14) + coord_cartesian(ylim = c(0, 0.28))


ggplot(countis_ind[celltype =="Non_Classical_Monocytes" ], aes(severity2, proportion, fill= celltype2)) +
  stat_summary(fill= "#542E42", geom = "bar", fun = "mean", na.rm = TRUE,
               position = position_dodge(width = 0.4)) +
  geom_errorbar(stat = "summary", fun.data = "mean_sdl", 
                fun.args = list(mult = 1),
                alpha = 0.9,
                position =  position_dodge(width = 0.9), color = "black", linewidth = 0.7, width = 0.4) +
  geom_point( alpha = 0.7, col = "black", size = 3)  + facet_wrap(~species, scales = "free_x" , nrow = 1) + scale_y_continuous(labels = label_percent(), expand = c(0, 0), breaks = pretty_breaks(8)) + guides(fill = "none") + ylab("Non_Classical_Monocytes") + scale_x_discrete(guide = guide_axis(angle = -45)) + theme_classic(base_size = 14) + coord_cartesian(ylim = c(0, 0.07))


ggplot(countis_ind[celltype =="DC" ], aes(severity2, proportion, fill= celltype2)) +
  stat_summary(fill= "#4F6D7A", geom = "bar", fun = "mean", na.rm = TRUE,
               position = position_dodge(width = 0.4)) +
  geom_errorbar(stat = "summary", fun.data = "mean_sdl", 
                fun.args = list(mult = 1),
                alpha = 0.9,
                position =  position_dodge(width = 0.9), color = "black", linewidth = 0.7, width = 0.4) +
  geom_point( alpha = 0.7, col = "black", size = 3)  + facet_wrap(~species, scales = "free_x" , nrow = 1) + scale_y_continuous(labels = label_percent(), expand = c(0, 0), breaks = pretty_breaks(8)) + guides(fill = "none") + ylab("DC") + scale_x_discrete(guide = guide_axis(angle = -45)) + theme_classic(base_size = 14) + coord_cartesian(ylim = c(0, 0.015))


ggplot(countis_ind[celltype =="NK_Cells" ], aes(severity2, proportion, fill= celltype2)) +
  stat_summary(fill= "#246A73", geom = "bar", fun = "mean", na.rm = TRUE,
               position = position_dodge(width = 0.4)) +
  geom_errorbar(stat = "summary", fun.data = "mean_sdl", 
                fun.args = list(mult = 1),
                alpha = 0.9,
                position =  position_dodge(width = 0.9), color = "black", linewidth = 0.7, width = 0.4) +
  geom_point( alpha = 0.7, col = "black", size = 3)  + facet_wrap(~species, scales = "free_x" , nrow = 1) + scale_y_continuous(labels = label_percent(), expand = c(0, 0), breaks = pretty_breaks(8)) + guides(fill = "none") + ylab("NK_Cells") + scale_x_discrete(guide = guide_axis(angle = -45)) + theme_classic(base_size = 14) + coord_cartesian(ylim = c(0, 0.09))


ggplot(countis_ind[celltype =="CD8+_T_Cells" ], aes(severity2, proportion, fill= celltype2)) +
  stat_summary(fill= "#368F8B", geom = "bar", fun = "mean", na.rm = TRUE,
               position = position_dodge(width = 0.4)) +
  geom_errorbar(stat = "summary", fun.data = "mean_sdl", 
                fun.args = list(mult = 1),
                alpha = 0.9,
                position =  position_dodge(width = 0.9), color = "black", linewidth = 0.7, width = 0.4) +
  geom_point( alpha = 0.7, col = "black", size = 3)  + facet_wrap(~species, scales = "free_x" , nrow = 1) + scale_y_continuous(labels = label_percent(), expand = c(0, 0), breaks = pretty_breaks(8)) + guides(fill = "none") + ylab("CD8+_T_Cells") + scale_x_discrete(guide = guide_axis(angle = -45)) + theme_classic(base_size = 14) + coord_cartesian(ylim = c(0, 0.20))


ggplot(countis_ind[celltype =="B_Cells" ], aes(severity2, proportion, fill= celltype2)) +
  stat_summary(fill= "#62C370", geom = "bar", fun = "mean", na.rm = TRUE,
               position = position_dodge(width = 0.4)) +
  geom_errorbar(stat = "summary", fun.data = "mean_sdl", 
                fun.args = list(mult = 1),
                alpha = 0.9,
                position =  position_dodge(width = 0.9), color = "black", linewidth = 0.7, width = 0.4) +
  geom_point( alpha = 0.7, col = "black", size = 3)  + facet_wrap(~species, scales = "free_x" , nrow = 1) + scale_y_continuous(labels = label_percent(), expand = c(0, 0), breaks = pretty_breaks(8)) + guides(fill = "none") + ylab("B_Cells") + scale_x_discrete(guide = guide_axis(angle = -45)) + theme_classic(base_size = 14) + coord_cartesian(ylim = c(0, 0.70))


ggplot(countis_ind[celltype =="Immature Neutrophils 2" ], aes(severity2, proportion, fill= celltype2)) +
  stat_summary(fill= "#6BE2FF", geom = "bar", fun = "mean", na.rm = TRUE,
               position = position_dodge(width = 0.4)) +
  geom_errorbar(stat = "summary", fun.data = "mean_sdl", 
                fun.args = list(mult = 1),
                alpha = 0.9,
                position =  position_dodge(width = 0.9), color = "black", linewidth = 0.7, width = 0.4) +
  geom_point( alpha = 0.7, col = "black", size = 3)  + facet_wrap(~species, scales = "free_x" , nrow = 1) + scale_y_continuous(labels = label_percent(), expand = c(0, 0), breaks = pretty_breaks(8)) + guides(fill = "none") + ylab("Immature Neutrophils 2") + scale_x_discrete(guide = guide_axis(angle = -45)) + theme_classic(base_size = 14) + coord_cartesian(ylim = c(0, 0.07))


ggplot(countis_ind[celltype =="Megakaryocyte" ], aes(severity2, proportion, fill= celltype2)) +
  stat_summary(fill= "#F97E44", geom = "bar", fun = "mean", na.rm = TRUE,
               position = position_dodge(width = 0.4)) +
  geom_errorbar(stat = "summary", fun.data = "mean_sdl", 
                fun.args = list(mult = 1),
                alpha = 0.9,
                position =  position_dodge(width = 0.9), color = "black", linewidth = 0.7, width = 0.4) +
  geom_point( alpha = 0.7, col = "black", size = 3)  + facet_wrap(~species, scales = "free_x" , nrow = 1) + scale_y_continuous(labels = label_percent(), expand = c(0, 0), breaks = pretty_breaks(8)) + guides(fill = "none") + ylab("Megakaryocyte") + scale_x_discrete(guide = guide_axis(angle = -45)) + theme_classic(base_size = 14) + coord_cartesian(ylim = c(0, 0.08))


ggplot(countis_ind[celltype =="Plasmablasts" ], aes(severity2, proportion, fill= celltype2)) +
  stat_summary(fill= "#7C6A0A", geom = "bar", fun = "mean", na.rm = TRUE,
               position = position_dodge(width = 0.4)) +
  geom_errorbar(stat = "summary", fun.data = "mean_sdl", 
                fun.args = list(mult = 1),
                alpha = 0.9,
                position =  position_dodge(width = 0.9), color = "black", linewidth = 0.7, width = 0.4) +
  geom_point( alpha = 0.7, col = "black", size = 3)  + facet_wrap(~species, scales = "free_x" , nrow = 1) + scale_y_continuous(labels = label_percent(), expand = c(0, 0), breaks = pretty_breaks(8)) + guides(fill = "none") + ylab("Plasmablasts") + scale_x_discrete(guide = guide_axis(angle = -45)) + theme_classic(base_size = 14) + coord_cartesian(ylim = c(0, 0.02))


ggplot(countis_ind[celltype =="Proliferating_Cells" ], aes(severity2, proportion, fill= celltype2)) +
  stat_summary(fill= "#FB3640", geom = "bar", fun = "mean", na.rm = TRUE,
               position = position_dodge(width = 0.4)) +
  geom_errorbar(stat = "summary", fun.data = "mean_sdl", 
                fun.args = list(mult = 1),
                alpha = 0.9,
                position =  position_dodge(width = 0.9), color = "black", linewidth = 0.7, width = 0.4) +
  geom_point( alpha = 0.7, col = "black", size = 3)  + facet_wrap(~species, scales = "free_x" , nrow = 1) + scale_y_continuous(labels = label_percent(), expand = c(0, 0), breaks = pretty_breaks(8)) + guides(fill = "none") + ylab("Proliferating_Cells") + scale_x_discrete(guide = guide_axis(angle = -45)) + theme_classic(base_size = 14) + coord_cartesian(ylim = c(0, 0.01))


ggplot(countis_ind[celltype =="Platelet" ], aes(severity2, proportion, fill= celltype2)) +
  stat_summary(fill= "#F7C548", geom = "bar", fun = "mean", na.rm = TRUE,
               position = position_dodge(width = 0.4)) +
  geom_errorbar(stat = "summary", fun.data = "mean_sdl", 
                fun.args = list(mult = 1),
                alpha = 0.9,
                position =  position_dodge(width = 0.9), color = "black", linewidth = 0.7, width = 0.4) +
  geom_point( alpha = 0.7, col = "black", size = 3)  + facet_wrap(~species, scales = "free_x" , nrow = 1) + scale_y_continuous(labels = label_percent(), expand = c(0, 0), breaks = pretty_breaks(8)) + guides(fill = "none") + ylab("Platelet") + scale_x_discrete(guide = guide_axis(angle = -45)) + theme_classic(base_size = 14) + coord_cartesian(ylim = c(0, 0.4))

#NEW Supllemental Figure 2

DotPlot(Blood_merged, features = c("CD14", "CCR2", "CX3CR1", "FLT3","CXCR2", "NKG7", "CD3E", "CD4","CD8A", "CD79B", "PRDM1", "PPBP", "PF4", "TOP2A"))+ theme(axis.text.x = element_text(angle = 90, vjust = 0.4, hjust = 1))
                           
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
#   [1] ggExtra_0.10.1          ggforce_0.4.1           reshape2_1.4.4          ggthemes_4.2.4          ggrepel_0.9.3          
# [6] patchwork_1.1.2         cowplot_1.1.1           Seurat_4.9.9.9058       SeuratObject_4.9.9.9091 sp_2.0-0               
# [11] dplyr_1.1.2             ggplot2_3.4.2           toolboxH_0.2.17         eulerr_7.0.0            testthat_3.1.10        
# [16] stringr_1.5.0           scales_1.2.1            readxl_1.4.3            RColorBrewer_1.1-3      png_0.1-8              
# [21] fdrtool_1.2.17          R.utils_2.12.2          R.oo_1.25.0             R.methodsS3_1.8.2       data.table_1.14.8      
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
