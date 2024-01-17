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