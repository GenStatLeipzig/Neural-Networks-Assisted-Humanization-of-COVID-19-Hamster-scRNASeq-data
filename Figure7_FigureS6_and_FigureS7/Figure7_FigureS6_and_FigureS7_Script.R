#Figure 7, Supplemental Figure 6 and Supplemental Figure 7

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

#Load Data
Blood_merged <- readRDS(".../h228_3_qc_integrated_rpca_cellanno_FILTEREDqc_FILTEREDpureCells.RDS")

#7A
diffexpr = read_excel2("Data/h245_7_de_list.xlsx")

GenesInflammation <- c("ACSL3", "AGER", "ALDH2", "APOBEC1", "APOC1", "APOE", "ARG1", "ATOX1", "B2M", "C1R", "C1S", "CCL12", "CCL17", "CCL2", "CCL3", "CCL4", "CCL5", "CCL7", "CCL8", "CCND1", "CCR5", "CD163", "CD274", "CFI", "CLU", "CMPK2", "CSF1", "CTSD", "CXCL10", "CXCL11", "CXCL12", "CXCL13", "CXCL17", "CYP2A5", "DCXR", "DDX60", "EIF2AK2", "FAM89B", "FASLG", "FCGR4", "FGG", "FTL", "GBP2", "GBP5", "GBP7", "GPR84", "GRAMD1C", "H4C9", "HERC6", "HP", "IDO1", "IFI16", "IFI35", "IFI47", "IFIT2", "IFIT3", "IFNAR1", "IFNG", "IGLL1", "IL10", "IL12A", "IL18", "IL1A", "IL1B", "IL1RN", "IL2RA", "IL6", "IL6R", "IRF7", "IRF9", "ISG15", "ISG20", "LBP", "LGALS3BP", "LRG1", "LTA", "LY6E", "LY6E", "MAFB", "MITD1", "MNDA", "MS4A7", "MX1", "MX2", "NLRC5", "NMI", "OAS2", "OASL2", "OLFM4", "PARP9", "PLAAT3", "PLAC8", "PLBD1", "PPBP", "PSME1", "RGCC", "RGS1", "RPS10", "RPS15", "RPS18", "RPS25", "RSAD2", "RTP4", "S100A4", "SCGB2A2", "SEC14L3", "SERPINE1", "SERPING1", "SIGLEC1", "SLAMF7", "SLAMF9", "SLFN11", "SLFN5", "SOCS1", "SP100", "TGFB3", "TM4SF19", "TNF", "TNFRSF11A", "TNFSF10", "TNFSF13B", "TREM1", "TST", "TUBB4B", "UBA52", "UBE2L6", "UBE2L6", "UPP1", "USP18", "XCL1", "YBX3", "TNFSF14", "S100A8", "S100A9", "IL6", "CXCL5", "CXCL2", "CD274", "CD177", "PDCD1")

Goitable = diffexpr[symbol %in% GenesInflammation]
uniqueN(Goitable$symbol)

Goitable$significant <- "no"
Goitable$significant[Goitable$adj.P.Val <= 0.05] <- "yes"
setDT(Goitable)
setorder(Goitable, -logFC)
Goitable[,symbol:= factor(symbol, levels = unique(symbol))]
dput(Goitable$celltype %>% unique())
Goitable[,celltype:=factor(celltype, levels = c("Classical_Monocytes", "Non_Classical_Monocytes","DC", "Immature Neutrophils 1", "Immature Neutrophils 2", "Neutrophils", "NK_Cells", "CD4+_T_Cells", "CD8+_T_Cells", "B_Cells", "Plasmablasts", "Platelet", "Megakaryocyte",  "Proliferating_Cells"))]

GOI_Figure7 <- dplyr::filter(Goitable, (Goitable$significant=="yes"))
range(GOI_Figure7$logFC)


GOI_Sig_MA_2v0 <- dplyr::filter(Goitable, c((Goitable$significant=="yes") & (Goitable$species=="hamsterMA") & (Goitable$contrast1=="d2")))
ggplot(GOI_Sig_MA_2v0, aes(y=symbol,x=celltype,size=-log10(adj.P.Val),color=logFC)) +
  geom_point() +
  scale_size_continuous(name='adj. p-value',
                        breaks=c(2,4,6),
                        labels=c("1e-2","1e-4", "<1e-6"),
                        limits=c(0,NA)) +
  scale_color_gradient2(low='blue3',mid='gray',high='red3',
                        limits=c(-11,12),oob=scales::squish) +
  theme_minimal(base_size=10) +
  guides(color=guide_colorbar(barheight=4)) +
  scale_x_discrete(guide = guide_axis(angle = -45)) + scale_y_discrete(limits=rev) +
  labs(x='',y='')


GOI_Sig_PR_LD_2v0 <- dplyr::filter(Goitable, c((Goitable$significant=="yes") & (Goitable$species=="hamsterPR") & (Goitable$contrast1=="ld_D2")))
ggplot(GOI_Sig_PR_LD_2v0, aes(y=symbol,x=celltype,size=-log10(adj.P.Val),color=logFC)) +
  geom_point() +
  scale_size_continuous(name='adj. p-value',
                        breaks=c(2,4,6),
                        labels=c("1e-2","1e-4", "<1e-6"),
                        limits=c(0,NA)) +
  scale_color_gradient2(low='blue3',mid='gray',high='red3',
                        limits=c(-11,12),oob=scales::squish) +
  theme_minimal(base_size=10) +
  guides(color=guide_colorbar(barheight=4)) +
  scale_x_discrete(guide = guide_axis(angle = -45)) + scale_y_discrete(limits=rev) +
  labs(x='',y='')


GOI_Sig_PR_HD_2v0 <- dplyr::filter(Goitable, c((Goitable$significant=="yes") & (Goitable$species=="hamsterPR") & (Goitable$contrast1=="hd_D2")))
ggplot(GOI_Sig_PR_HD_2v0, aes(y=symbol,x=celltype,size=-log10(adj.P.Val),color=logFC)) +
  geom_point() +
  scale_size_continuous(name='adj. p-value',
                        breaks=c(2,4,6),
                        labels=c("1e-2","1e-4", "<1e-6"),
                        limits=c(0,NA)) +
  scale_color_gradient2(low='blue3',mid='gray',high='red3',
                        limits=c(-11,12),oob=scales::squish) +
  theme_minimal(base_size=10) +
  guides(color=guide_colorbar(barheight=4)) +
  scale_x_discrete(guide = guide_axis(angle = -45)) + scale_y_discrete(limits=rev) +
  labs(x='',y='')


GOI_Sig_HU_3v0 <- dplyr::filter(Goitable, c((Goitable$significant=="yes") & (Goitable$species=="human") & (Goitable$contrast1=="3")))
ggplot(GOI_Sig_HU_3v0, aes(y=symbol,x=celltype,size=-log10(adj.P.Val),color=logFC)) +
  geom_point() +
  scale_size_continuous(name='adj. p-value',
                        breaks=c(2,4,6),
                        labels=c("1e-2","1e-4", "<1e-6"),
                        limits=c(0,NA)) +
  scale_color_gradient2(low='blue3',mid='gray',high='red3',
                        limits=c(-11,12),oob=scales::squish) +
  theme_minimal(base_size=10) +
  guides(color=guide_colorbar(barheight=4)) +
  scale_x_discrete(guide = guide_axis(angle = -45)) + scale_y_discrete(limits=rev) +
  labs(x='',y='')


GOI_Sig_HU_7v0 <- dplyr::filter(Goitable, c((Goitable$significant=="yes") & (Goitable$species=="human") & (Goitable$contrast1=="7")))
ggplot(GOI_Sig_HU_7v0, aes(y=symbol,x=celltype,size=-log10(adj.P.Val),color=logFC)) +
  geom_point() +
  scale_size_continuous(name='adj. p-value',
                        breaks=c(2,4,6),
                        labels=c("1e-2","1e-4", "<1e-6"),
                        limits=c(0,NA)) +
  scale_color_gradient2(low='blue3',mid='gray',high='red3',
                        limits=c(-11,12),oob=scales::squish) +
  theme_minimal(base_size=10) +
  guides(color=guide_colorbar(barheight=4)) +
  scale_x_discrete(guide = guide_axis(angle = -45)) + scale_y_discrete(limits=rev) +
  labs(x='',y='')






#S6A,S6B
GOI_Sig_MA_3v0 <- dplyr::filter(Goitable, c((Goitable$significant=="yes") & (Goitable$species=="hamsterMA") & (Goitable$contrast1=="d3")))
ggplot(GOI_Sig_MA_3v0, aes(y=symbol,x=celltype,size=-log10(adj.P.Val),color=logFC)) +
  geom_point() +
  scale_size_continuous(name='adj. p-value',
                        breaks=c(2,4,6),
                        labels=c("1e-2","1e-4", "<1e-6"),
                        limits=c(0,NA)) +
  scale_color_gradient2(low='blue3',mid='gray',high='red3',
                        limits=c(-11,12),oob=scales::squish) +
  theme_minimal(base_size=10) +
  guides(color=guide_colorbar(barheight=4)) +
  scale_x_discrete(guide = guide_axis(angle = -45)) + scale_y_discrete(limits=rev) +
  labs(x='',y='')


GOI_Sig_MA_5v0 <- dplyr::filter(Goitable, c((Goitable$significant=="yes") & (Goitable$species=="hamsterMA") & (Goitable$contrast1=="d5")))
ggplot(GOI_Sig_MA_5v0, aes(y=symbol,x=celltype,size=-log10(adj.P.Val),color=logFC)) +
  geom_point() +
  scale_size_continuous(name='adj. p-value',
                        breaks=c(2,4,6),
                        labels=c("1e-2","1e-4", "<1e-6"),
                        limits=c(0,NA)) +
  scale_color_gradient2(low='blue3',mid='gray',high='red3',
                        limits=c(-11,12),oob=scales::squish) +
  theme_minimal(base_size=10) +
  guides(color=guide_colorbar(barheight=4)) +
  scale_x_discrete(guide = guide_axis(angle = -45)) + scale_y_discrete(limits=rev) +
  labs(x='',y='')


GOI_Sig_MA_14v0 <- dplyr::filter(Goitable, c((Goitable$significant=="yes") & (Goitable$species=="hamsterMA") & (Goitable$contrast1=="e14")))
ggplot(GOI_Sig_MA_14v0, aes(y=symbol,x=celltype,size=-log10(adj.P.Val),color=logFC)) +
  geom_point() +
  scale_size_continuous(name='adj. p-value',
                        breaks=c(2,4,6),
                        labels=c("1e-2","1e-4", "<1e-6"),
                        limits=c(0,NA)) +
  scale_color_gradient2(low='blue3',mid='gray',high='red3',
                        limits=c(-11,12),oob=scales::squish) +
  theme_minimal(base_size=10) +
  guides(color=guide_colorbar(barheight=4)) +
  scale_x_discrete(guide = guide_axis(angle = -45)) + scale_y_discrete(limits=rev) +
  labs(x='',y='')


GOI_Sig_PR_LD_3v0 <- dplyr::filter(Goitable, c((Goitable$significant=="yes") & (Goitable$species=="hamsterPR") & (Goitable$contrast1=="ld_D3")))
ggplot(GOI_Sig_PR_LD_3v0, aes(y=symbol,x=celltype,size=-log10(adj.P.Val),color=logFC)) +
  geom_point() +
  scale_size_continuous(name='adj. p-value',
                        breaks=c(2,4,6),
                        labels=c("1e-2","1e-4", "<1e-6"),
                        limits=c(0,NA)) +
  scale_color_gradient2(low='blue3',mid='gray',high='red3',
                        limits=c(-11,12),oob=scales::squish) +
  theme_minimal(base_size=10) +
  guides(color=guide_colorbar(barheight=4)) +
  scale_x_discrete(guide = guide_axis(angle = -45)) + scale_y_discrete(limits=rev) +
  labs(x='',y='')


GOI_Sig_PR_HD_3v0 <- dplyr::filter(Goitable, c((Goitable$significant=="yes") & (Goitable$species=="hamsterPR") & (Goitable$contrast1=="hd_D3")))
ggplot(GOI_Sig_PR_HD_3v0, aes(y=symbol,x=celltype,size=-log10(adj.P.Val),color=logFC)) +
  geom_point() +
  scale_size_continuous(name='adj. p-value',
                        breaks=c(2,4,6),
                        labels=c("1e-2","1e-4", "<1e-6"),
                        limits=c(0,NA)) +
  scale_color_gradient2(low='blue3',mid='gray',high='red3',
                        limits=c(-11,12),oob=scales::squish) +
  theme_minimal(base_size=10) +
  guides(color=guide_colorbar(barheight=4)) +
  scale_x_discrete(guide = guide_axis(angle = -45)) + scale_y_discrete(limits=rev) +
  labs(x='',y='')


GOI_Sig_HU_4u5v0 <- dplyr::filter(Goitable, c((Goitable$significant=="yes") & (Goitable$species=="human") & (Goitable$contrast1=="4u5")))
ggplot(GOI_Sig_HU_4u5v0, aes(y=symbol,x=celltype,size=-log10(adj.P.Val),color=logFC)) +
  geom_point() +
  scale_size_continuous(name='adj. p-value',
                        breaks=c(2,4,6),
                        labels=c("1e-2","1e-4", "<1e-6"),
                        limits=c(0,NA)) +
  scale_color_gradient2(low='blue3',mid='gray',high='red3',
                        limits=c(-11,12),oob=scales::squish) +
  theme_minimal(base_size=10) +
  guides(color=guide_colorbar(barheight=4)) +
  scale_x_discrete(guide = guide_axis(angle = -45)) + scale_y_discrete(limits=rev) +
  labs(x='',y='')



#7B and S7A,B
Blood_merged@active.ident = factor(Blood_merged$uniform_name_overview3)
severe_inflammatory <- c("CD177.1", "ZDHHC19.1", "PLAC8.1", "IFI6.2", "MT2A.2", "SERPING1.2", "RSAD2.2", "GBP1.2", "ISG15.2", "LY6E.2", "S100A12.2", "GBP5.2", "XAF1.2", "IFI44L.2", "IFIT3.2", "IFITM3.2", "SELL.2", "TRIM22.2", "EPSTI1.2", "FCER1G.2", "IFI44.2", "CST7.1", "OAS1.2", "TNFSF13B.2", "SLFN5", "IFIT2.2", "MX1.2", "PLSCR1.2", "MMP9.2", "LAP3", "OAS2.2", "IFIT1.2", "TAP1.2", "SAMD9L.2", "HIST1H2AC.1", "APOL6.2", "OASL.1", "DDX58.1", "WSB1.1", "RPL28.1", "IL1RN.2", "MCEMP1.1", "HERC5.1", "IRF7.2", "CKAP4.1", "UBE2L6.1", "RNF213.2", "GBP2.2", "CLEC4D.1", "LILRA6.1", "SHISA5.2", "PHF11.1", "IFI16.2", "OAS3.1", "PIM1.1", "SERPINB1.1", "GYG1.1", "IFITM1.1", "PSTPIP2", "TMSB10.1", "GBP4", "NFKBIA.2", "EIF2AK2.1", "FFAR2.2", "PARP9.1", "NT5C3A", "TNFSF10.2", "MYL12A.1", "PARP14.1", "DDX60", "DDX60L.2", "STAT1.1", "XRN1", "SAMD9.2", "SNX3", "HIST1H2BD", "PGD.1", "FCGR1A", "JUN", "CD44", "PLP2.1", "CARD16.2", "LIMK2.2", "CYSTM1.1", "HMGB2.1", "IFIH1.2", "C1orf162", "KLF4", "STAT2", "ZBP1", "ALPL.1", "GSTK1", "SP110.2", "ANXA1", "ANXA3.1", "METTL9", "PML", "NUCB1", "CEACAM1", "SECTM1", "PSMB9.1", "SAT1", "RNF10.1", "MAPK14.1", "LILRA5", "C3AR1", "ZCCHC2", "SAMHD1", "LGALS9", "CR1", "C4orf3", "IRF1.1", "CD63", "GRINA", "TXN.1", "MX2.1", "GAPDH.1", "RBCK1", "NBN.2", "NMI.1", "CAST", "HIST2H2BE", "NTNG2", "SP100.1", "CASP1.1", "TNFAIP6", "PLEK.1", "LMNB1", "NFIL3", "APOL2", "BST1", "NUB1", "PIK3AP1", "CCR1.2", "ALOX5AP.1", "EMB", "ISG20", "TMEM123.1", "ADAR", "GIMAP4", "SPTLC2", "SAMSN1.1", "SPI1", "GRN.1", "APOBEC3A", "GLRX", "CD53.1", "TRIM38", "CD82", "SH3GLB1", "VIM.1", "ADM.1", "CASP4.1", "S100A6.1", "RAC2.1", "BAZ1A", "ADD3", "ITGAM", "LY96.2", "MSRB1.1", "DYSF.1", "MOB1A", "TPM3.1", "FGR", "ACSL1", "IL2RG.1", "CDKN2D", "FYB1.1", "CLEC4E.1", "HLA-F", "LILRB3", "RBMS1", "PLBD1", "FLOT1", "GNG5", "PTEN.1", "PROK2.1", "UBE2J1", "CD37.1", "KCNJ15.1", "BRI3.1", "HIF1A", "PRR13.1", "LRG1", "MAX", "RHOG", "NFE2", "STXBP2", "B4GALT5", "CAPZA1", "SERPINA1", "ZYX", "RGS19", "FKBP5", "GCA.1", "FKBP1A", "MTPN", "VNN2", "AC245128.3", "CD55", "CFL1.1")
severe_inflammatory2 <- severe_inflammatory[severe_inflammatory %in% rownames(Blood_merged)]

Blood_merged <- AddModuleScore(object = Blood_merged, features = severe_inflammatory2, name = "severe_inflammatory")

#S7A
Blood_merged_neut <- subset(Blood_merged, uniform_name_overview3 == "Neutrophils")
p1_data <- FeaturePlot(object = Blood_merged_neut, features = "severe_inflammatory1", split.by = "species", cols = c("lightblue", "red"), order = T,  raster = F, reduction = "umap.rpca.v2") & coord_cartesian(ylim = c(-6,5), xlim = c(-4,10)) & theme(legend.position = c(0.1,0.3))
p1_data

#S7B
DimPlot(Blood_merged_neut, group.by = c("rpca_clusters.v2"),label = T, repel = T, shuffle = T, reduction = "umap.rpca.v2") + NoLegend() + coord_cartesian(ylim = c(-6,5), xlim = c(-4,10))

#7B
Blood_merged$uniform_name_overview3Neutsplit <- ifelse(Blood_merged$rpca_clusters.v2 %in% c(2,3),  "NeutrophilsB", ifelse(Blood_merged$uniform_name_overview3 == "Neutrophils", "NeutrophilsA", Blood_merged$uniform_name_overview3))
table(Blood_merged$uniform_name_overview3, Blood_merged$uniform_name_overview3Neutsplit)
Blood_merged@active.ident = factor(Blood_merged$uniform_name_overview3Neutsplit)

dput(Blood_merged$severity %>% unique())
Blood_merged$severity[Blood_merged$severity =="D0"] <- "aControl"
Blood_merged$severity[Blood_merged$severity =="ld_D2"] <- "b2 dpi ld"
Blood_merged$severity[Blood_merged$severity =="ld_D3"] <- "c3 dpi ld"
Blood_merged$severity[Blood_merged$severity =="hd_D2"] <- "d2 dpi hd"
Blood_merged$severity[Blood_merged$severity =="hd_D3"] <- "e3 dpi hd"

Blood_merged@active.ident = factor(Blood_merged$uniform_name_overview3)
InnateCells <- subset(Blood_merged, idents = c("Classical_Monocytes", "Non_Classical_Monocytes","DC", "Immature Neutrophils 1", "Immature Neutrophils 2", "Neutrophils", "NK_Cells"))

CellorderVLN = c("Classical_Monocytes", "Non_Classical_Monocytes","DC", "Immature Neutrophils 1", "Immature Neutrophils 2", "NeutrophilsA", "NeutrophilsB", "NK_Cells", "CD4+_T_Cells", "CD8+_T_Cells", "B_Cells", "Plasmablasts", "Platelet", "Megakaryocyte",  "Proliferating_Cells")
my_levels <- c(CellorderVLN)

InnateCells@active.ident = factor(InnateCells$species)
HumanInnate <- subset(InnateCells, idents = "human")
MAInnate <- subset(InnateCells, idents = "hamsterMA")
PRInnate <- subset(InnateCells, idents = "hamsterPR")

HumanInnate@active.ident = factor(HumanInnate$uniform_name_overview3Neutsplit)
MAInnate@active.ident = factor(MAInnate$uniform_name_overview3Neutsplit)
PRInnate@active.ident = factor(PRInnate$uniform_name_overview3Neutsplit)

HumanInnate@active.ident <- factor(x = HumanInnate@active.ident, levels = my_levels)
MAInnate@active.ident <- factor(x = MAInnate@active.ident, levels = my_levels)
PRInnate@active.ident <- factor(x = PRInnate@active.ident, levels = my_levels)

p1_data <- VlnPlot(object = MAInnate, features = "severe_inflammatory1", split.by = "severity", cols = c("#4fb6ca", "#62205f", "#b9563f", "#a9845b", "#ecb27d"), raster = F) & theme(legend.position = c(0.9,1))
p1_data + scale_x_discrete(guide = guide_axis(angle = -45))


p2_data <- VlnPlot(object = PRInnate, features = "severe_inflammatory1", split.by = "severity", cols = c("#508ca7", "#9f2d55", "#86201e", "#341648", "#591d08"), raster = F, alpha = 0.1, pt.size = 0.1) & theme(legend.position = c(0.9,1))
p2_data + scale_x_discrete(guide = guide_axis(angle = -45))

p3_data <- VlnPlot(object = HumanInnate, features = "severe_inflammatory1", split.by = "severity", cols = c("#1a71b8", "#078d84", "#f7b817", "#e27d26", "#c63735"), raster = F , alpha = 0.1, pt.size = 0.1) & theme(legend.position = c(0.9,1))
p3_data + scale_x_discrete(guide = guide_axis(angle = -45))





#S7C
Blood_merged@active.ident = factor(Blood_merged$species)
Human <- subset(Blood_merged, species == "human")
Human@active.ident = factor(Human$uniform_name_overview3)
HumanNeut <- subset(Human, idents = "Neutrophils")

DimPlot(HumanNeut, group.by = c("rpca_clusters.v2"),label = T, repel = T, raster = F, shuffle = T, reduction = "umap.rpca.v2") + NoLegend()

HumanNeut@active.ident = factor(HumanNeut$severity)
HumanNeutSev <- subset(HumanNeut, idents = c("5","7"))
DimPlot(HumanNeutSev, group.by = c("rpca_clusters.v2"),label = T, repel = T, raster = F, shuffle = T, reduction = "umap.rpca.v2")

HumanNeutSev@active.ident = factor(HumanNeutSev$rpca_clusters.v2)
HumanNeutSev <- RenameIdents(HumanNeutSev, "0" = "Classical_Neutrophils", "1" = "Classical_Neutrophils", "2" = "severe_infl_Neutrophils", "3" = "severe_infl_Neutrophils","4" = "Classical_Neutrophils", "8" = "Classical_Neutrophils", "9" = "Classical_Neutrophils", "11" = "Classical_Neutrophils","12" = "Classical_Neutrophils","18" = "Classical_Neutrophils","20" = "Classical_Neutrophils","22" = "Classical_Neutrophils","23" = "Classical_Neutrophils","27" = "Classical_Neutrophils")
DimPlot(HumanNeutSev,label = T, repel = T, raster = F, shuffle = T, reduction = "umap.rpca.v2", cols = c("lightblue", "red"))

HumanNeutSev$NeutSub <- paste(HumanNeutSev@active.ident)
levels(HumanNeutSev@active.ident)
HumanNeutSev@active.ident = factor(HumanNeutSev$NeutSub)
diffExprHUNeut <- FindMarkers(object = HumanNeutSev, ident.1 = "severe_infl_Neutrophils", ident.2 = "Classical_Neutrophils")

diffExprHUNeut <- tibble::rownames_to_column(diffExprHUNeut, "Gene")
diffExprHUNeut$Bonferoni_adj <- "FALSE"
diffExprHUNeut$Bonferoni_adj[diffExprHUNeut$p_val_adj < 0.05] <- "TRUE"

diffExprHUNeut$lable <- ifelse(abs(diffExprHUNeut$avg_log2FC) >= 1 | -log10(diffExprHUNeut$p_val_adj) >= 150, diffExprHUNeut$Gene, "")

diffExprHUNeut2 <- dplyr::filter(diffExprHUNeut, diffExprHUNeut$Bonferoni_adj == "TRUE")

ggplot(data=diffExprHUNeut2, aes(x=avg_log2FC, y=-log10(p_val_adj), label=lable)) + geom_point() + geom_text_repel(size = 10, max.overlaps = 7) + geom_vline(xintercept=c(-1, 1), col="red") + theme_classic()


#S7D
GOI <- diffExprHUNeut2 %>%
  filter(Gene %in% c("ACSL3", "AGER", "ALDH2", "APOBEC1", "APOC1", "APOE", "ARG1", "ATOX1", "B2M", "C1R", "C1S", "CCL12", "CCL17", "CCL2", "CCL3", "CCL4", "CCL5", "CCL7", "CCL8", "CCND1", "CCR5", "CD163", "CD274", "CFI", "CLU", "CMPK2", "CSF1", "CTSD", "CXCL10", "CXCL11", "CXCL12", "CXCL13", "CXCL17", "CYP2A5", "DCXR", "DDX60", "EIF2AK2", "FAM89B", "FASLG", "FCGR4", "FGG", "FTL", "GBP2", "GBP5", "GBP7", "GPR84", "GRAMD1C", "H4C9", "HERC6", "HP", "IDO1", "IFI16", "IFI35", "IFI47", "IFIT2", "IFIT3", "IFNAR1", "IFNG", "IGLL1", "IL10", "IL12A", "IL18", "IL1A", "IL1B", "IL1RN", "IL2RA", "IL6", "IL6R", "IRF7", "IRF9", "ISG15", "ISG20", "LBP", "LGALS3BP", "LRG1", "LTA", "LY6E", "LY6E", "MAFB", "MITD1", "MNDA", "MS4A7", "MX1", "MX2", "NLRC5", "NMI", "OAS2", "OASL2", "OLFM4", "PARP9", "PLAAT3", "PLAC8", "PLBD1", "PPBP", "PSME1", "RGCC", "RGS1", "RPS10", "RPS15", "RPS18", "RPS25", "RSAD2", "RTP4", "S100A4", "SCGB2A2", "SEC14L3", "SERPINE1", "SERPING1", "SIGLEC1", "SLAMF7", "SLAMF9", "SLFN11", "SLFN5", "SOCS1", "SP100", "TGFB3", "TM4SF19", "TNF", "TNFRSF11A", "TNFSF10", "TNFSF13B", "TREM1", "TST", "TUBB4B", "UBA52", "UBE2L6", "UBE2L6", "UPP1", "USP18", "XCL1", "YBX3", "TNFSF14", "S100A8", "S100A9", "IL6", "CXCL5", "CXCL2", "CD274", "CD177", "PDCD1"))

ggplot(data=diffExprHUNeut2, aes(x=avg_log2FC, y=-log10(p_val_adj), label=Gene)) + geom_point() + geom_point(data=GOI, aes(x=avg_log2FC,y=-log10(p_val_adj)), col="blue") + geom_vline(xintercept=c(-1, 1), col="red") + geom_text_repel(size = 10, max.overlaps = 12, data = diffExprHUNeut2 %>%   filter(Gene %in% c("ACSL3", "AGER", "ALDH2", "APOBEC1", "APOC1", "APOE", "ARG1", "ATOX1", "B2M", "C1R", "C1S", "CCL12", "CCL17", "CCL2", "CCL3", "CCL4", "CCL5", "CCL7", "CCL8", "CCND1", "CCR5", "CD163", "CD274", "CFI", "CLU", "CMPK2", "CSF1", "CTSD", "CXCL10", "CXCL11", "CXCL12", "CXCL13", "CXCL17", "CYP2A5", "DCXR", "DDX60", "EIF2AK2", "FAM89B", "FASLG", "FCGR4", "FGG", "FTL", "GBP2", "GBP5", "GBP7", "GPR84", "GRAMD1C", "H4C9", "HERC6", "HP", "IDO1", "IFI16", "IFI35", "IFI47", "IFIT2", "IFIT3", "IFNAR1", "IFNG", "IGLL1", "IL10", "IL12A", "IL18", "IL1A", "IL1B", "IL1RN", "IL2RA", "IL6", "IL6R", "IRF7", "IRF9", "ISG15", "ISG20", "LBP", "LGALS3BP", "LRG1", "LTA", "LY6E", "LY6E", "MAFB", "MITD1", "MNDA", "MS4A7", "MX1", "MX2", "NLRC5", "NMI", "OAS2", "OASL2", "OLFM4", "PARP9", "PLAAT3", "PLAC8", "PLBD1", "PPBP", "PSME1", "RGCC", "RGS1", "RPS10", "RPS15", "RPS18", "RPS25", "RSAD2", "RTP4", "S100A4", "SCGB2A2", "SEC14L3", "SERPINE1", "SERPING1", "SIGLEC1", "SLAMF7", "SLAMF9", "SLFN11", "SLFN5", "SOCS1", "SP100", "TGFB3", "TM4SF19", "TNF", "TNFRSF11A", "TNFSF10", "TNFSF13B", "TREM1", "TST", "TUBB4B", "UBA52", "UBE2L6", "UBE2L6", "UPP1", "USP18", "XCL1", "YBX3", "TNFSF14", "S100A8", "S100A9", "IL6", "CXCL5", "CXCL2", "CD274", "CD177", "PDCD1")), color="blue") + theme_classic()



#S7E
severe_inflammatory <- c("CD177.1", "ZDHHC19.1", "PLAC8.1", "IFI6.2", "MT2A.2", "SERPING1.2", "RSAD2.2", "GBP1.2", "ISG15.2", "LY6E.2", "S100A12.2", "GBP5.2", "XAF1.2", "IFI44L.2", "IFIT3.2", "IFITM3.2", "SELL.2", "TRIM22.2", "EPSTI1.2", "FCER1G.2", "IFI44.2", "CST7.1", "OAS1.2", "TNFSF13B.2", "SLFN5", "IFIT2.2", "MX1.2", "PLSCR1.2", "MMP9.2", "LAP3", "OAS2.2", "IFIT1.2", "TAP1.2", "SAMD9L.2", "HIST1H2AC.1", "APOL6.2", "OASL.1", "DDX58.1", "WSB1.1", "RPL28.1", "IL1RN.2", "MCEMP1.1", "HERC5.1", "IRF7.2", "CKAP4.1", "UBE2L6.1", "RNF213.2", "GBP2.2", "CLEC4D.1", "LILRA6.1", "SHISA5.2", "PHF11.1", "IFI16.2", "OAS3.1", "PIM1.1", "SERPINB1.1", "GYG1.1", "IFITM1.1", "PSTPIP2", "TMSB10.1", "GBP4", "NFKBIA.2", "EIF2AK2.1", "FFAR2.2", "PARP9.1", "NT5C3A", "TNFSF10.2", "MYL12A.1", "PARP14.1", "DDX60", "DDX60L.2", "STAT1.1", "XRN1", "SAMD9.2", "SNX3", "HIST1H2BD", "PGD.1", "FCGR1A", "JUN", "CD44", "PLP2.1", "CARD16.2", "LIMK2.2", "CYSTM1.1", "HMGB2.1", "IFIH1.2", "C1orf162", "KLF4", "STAT2", "ZBP1", "ALPL.1", "GSTK1", "SP110.2", "ANXA1", "ANXA3.1", "METTL9", "PML", "NUCB1", "CEACAM1", "SECTM1", "PSMB9.1", "SAT1", "RNF10.1", "MAPK14.1", "LILRA5", "C3AR1", "ZCCHC2", "SAMHD1", "LGALS9", "CR1", "C4orf3", "IRF1.1", "CD63", "GRINA", "TXN.1", "MX2.1", "GAPDH.1", "RBCK1", "NBN.2", "NMI.1", "CAST", "HIST2H2BE", "NTNG2", "SP100.1", "CASP1.1", "TNFAIP6", "PLEK.1", "LMNB1", "NFIL3", "APOL2", "BST1", "NUB1", "PIK3AP1", "CCR1.2", "ALOX5AP.1", "EMB", "ISG20", "TMEM123.1", "ADAR", "GIMAP4", "SPTLC2", "SAMSN1.1", "SPI1", "GRN.1", "APOBEC3A", "GLRX", "CD53.1", "TRIM38", "CD82", "SH3GLB1", "VIM.1", "ADM.1", "CASP4.1", "S100A6.1", "RAC2.1", "BAZ1A", "ADD3", "ITGAM", "LY96.2", "MSRB1.1", "DYSF.1", "MOB1A", "TPM3.1", "FGR", "ACSL1", "IL2RG.1", "CDKN2D", "FYB1.1", "CLEC4E.1", "HLA-F", "LILRB3", "RBMS1", "PLBD1", "FLOT1", "GNG5", "PTEN.1", "PROK2.1", "UBE2J1", "CD37.1", "KCNJ15.1", "BRI3.1", "HIF1A", "PRR13.1", "LRG1", "MAX", "RHOG", "NFE2", "STXBP2", "B4GALT5", "CAPZA1", "SERPINA1", "ZYX", "RGS19", "FKBP5", "GCA.1", "FKBP1A", "MTPN", "VNN2", "AC245128.3", "CD55", "CFL1.1")
severe_inflammatory2 <- severe_inflammatory[severe_inflammatory %in% rownames(Blood_merged)]
severe_inflammatory2

GenesInflammation_2 <- GenesInflammation[GenesInflammation %in% rownames(Blood_merged)]
GenesInflammation_2

venn2(GenesInflammation, severe_inflammatory)
venn2(GenesInflammation_2, severe_inflammatory2)

