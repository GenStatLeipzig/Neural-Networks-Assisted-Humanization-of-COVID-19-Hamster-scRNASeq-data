#Figure5 A-C and Supplemental Figure 5A-B


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
library(scales)


#Load Data
Blood_merged <- readRDS(".../h228_3_qc_integrated_rpca_cellanno_FILTEREDqc_FILTEREDpureCells.RDS")



#5A

eta1_tab2 = fread("Data/h245_7_eta1_tab.txt")

eta1_tab2[, contrast1v2 := factor(contrast1,
                                  levels = c("3", "4u5", "7", "d2", "d3", "d5", "e14","ld_D2", "ld_D3", "hd_D2", "hd_D3"),
                                  labels = c("WHO3", "WHO4&5", "WHO7", "2 dpi", "3 dpi", "5 dpi", "14 dpi","2 dpi ld", "3 dpi ld", "2 dpi hd", "3 dpi hd"))]




eta1_tab2[, celltype2 := factor(celltype, levels =  c("Classical_Monocytes", "Non_Classical_Monocytes", "Immature Neutrophils 1",  "Neutrophils", "NK_Cells", "CD4+_T_Cells", "CD8+_T_Cells", "B_Cells","DC", "Platelet",  "Plasmablasts","Megakaryocyte", "Immature Neutrophils 2", "Proliferating_Cells") %>% rev())]


p_eta1v2 = ggplot(eta1_tab2, aes(contrast1v2, celltype2,label = label, fill =eta1)) + geom_tile() +
  theme_minimal(base_size = 14) +
  facet_grid(.~species, scales = 'free', space = "free") + geom_text() + scale_fill_gradient_tableau(palette = "Red-Gold") + ylab("") + xlab("") + guides(fill = "none") + ggtitle("Percent transcriptome changed compared with control")
p_eta1v2




#5B

data_fig5b = fread("Data/h265_8_plotdata_FDR20_dotplot_effectsizes.txt.gz")

Cellorder5b = c("Classical\nMonocytes", "Non Classical\nMonocytes","DC", "Immature\nNeutrophils 1", "Neutrophils", "NK Cells", "CD4+ T Cells", "CD8+ T Cells", "B Cells")


data_fig5b$human_contrast1v2[data_fig5b$human_contrast1 =="3"] <- "WOS 3"
data_fig5b$human_contrast1v2[data_fig5b$human_contrast1 =="7"] <- "WOS 7"
data_fig5b$human_contrast1v2[data_fig5b$human_contrast1 =="4u5"] <- "WOS 4&5"

data_fig5b$hamster_contrast1v3[data_fig5b$hamster_contrast1v2 =="d2"] <- "c2 dpi"
data_fig5b$hamster_contrast1v3[data_fig5b$hamster_contrast1v2 =="d3"] <- "b3 dpi"
data_fig5b$hamster_contrast1v3[data_fig5b$hamster_contrast1v2 =="d5"] <- "a5 dpi"
data_fig5b$hamster_contrast1v3[data_fig5b$hamster_contrast1v2 =="e14"] <- "14 dpi"

data_fig5b$hamster_contrast1v3[data_fig5b$hamster_contrast1v2 =="ld_D2"] <- "d2 dpi ld"
data_fig5b$hamster_contrast1v3[data_fig5b$hamster_contrast1v2 =="ld_D3"] <- "c3 dpi ld"
data_fig5b$hamster_contrast1v3[data_fig5b$hamster_contrast1v2 =="hd_D2"] <- "b2 dpi hd"
data_fig5b$hamster_contrast1v3[data_fig5b$hamster_contrast1v2 =="hd_D3"] <- "a3 dpi hd"

data_fig5b$celltype3 <- factor(data_fig5b$celltype3, levels = c(Cellorder5b))


p2=ggplot(data_fig5b[human_contrast1v2 != "WOS 4&5" & n_intersectgenes>0], aes( human_contrast1v2, hamster_contrast1v3, col = -intersectgenes_sign_samedir_proz , size = n_signif_genes , alpha = p_binom<=0.05)) + geom_point() +
  facet_grid(hamster~celltype3, space = "free", scales = "free") +
  scale_size_continuous(range = c(2,9), breaks = c(5, 50, 100,250, 400, 1000))+
  theme_minimal(base_size = 14) +
  labs(size = "Number overlapping significant\ngenes (FDR20%)",
       color = "Percent human effect direction\n same in hamster",
       alpha = "Significant difference from 50% of\n same effect direction in hamster") +
  theme(legend.position = "top",legend.direction = "vertical") +
  scale_color_gradient2_tableau(guide = "legend", breaks = (c(0:5, 2.5)/-5)%>% sort(), labels = (c(0:5, 2.5)/5)%>% sort() %>% rev() *100,limits = c(-1,0)) +
  guides(col = guide_legend(nrow = 1,override.aes = list(size=5)),
         size = guide_legend(ncol = 2),
         alpha = guide_legend(override.aes = list(size=5))) +
  ylab("") +xlab("")

p2


#S5A
p3=ggplot(data_fig5b[human_contrast1v2 == "WOS 4&5" & n_intersectgenes>0], aes( human_contrast1v2, hamster_contrast1v3, col = -intersectgenes_sign_samedir_proz , size = n_signif_genes , alpha = p_binom<=0.05)) + geom_point() +
  facet_grid(hamster~celltype3, space = "free", scales = "free") +
  scale_size_continuous(range = c(2,9), breaks = c(5, 50, 100,250, 400, 1000))+
  theme_minimal(base_size = 14) +
  labs(size = "Number overlapping significant\ngenes (FDR20%)",
       color = "Percent human effect direction\n same in hamster",
       alpha = "Significant difference from 50% of\n same effect direction in hamster") +
  theme(legend.position = "top",legend.direction = "vertical") +
  scale_color_gradient2_tableau(guide = "legend", breaks = (c(0:5, 2.5)/-5)%>% sort(), labels = (c(0:5, 2.5)/5)%>% sort() %>% rev() *100,limits = c(-1,0)) +
  guides(col = guide_legend(nrow = 1,override.aes = list(size=5)),
         size = guide_legend(ncol = 2),
         alpha = guide_legend(override.aes = list(size=5))) +
  ylab("") +xlab("")

p3


#5C
data_fig5c = fread("Data/h245_7_inputdata_plot_h245_7_overlapping_pathways_no4u5_smaller.txt.gz")

data_fig5c$human2[data_fig5c$human =="3"] <- "WOS 3"
data_fig5c$human2[data_fig5c$human =="7"] <- "WOS 7"
data_fig5c$human2[data_fig5c$human =="4u5"] <- "WOS 4&5"

data_fig5c$hamster2[data_fig5c$hamster =="d2"] <- "2 dpi"
data_fig5c$hamster2[data_fig5c$hamster =="d3"] <- "3 dpi"
data_fig5c$hamster2[data_fig5c$hamster =="d5"] <- "5 dpi"
data_fig5c$hamster2[data_fig5c$hamster =="e14"] <- "14 dpi"

data_fig5c$hamster2[data_fig5c$hamster =="ld_D2"] <- "2 dpi ld"
data_fig5c$hamster2[data_fig5c$hamster =="ld_D3"] <- "3 dpi ld"
data_fig5c$hamster2[data_fig5c$hamster =="hd_D2"] <- "2 dpi hd"
data_fig5c$hamster2[data_fig5c$hamster =="hd_D3"] <- "3 dpi hd"

data_fig5c[,human2:=factor(human2)]
dput(data_fig5c$hamster2 %>% unique())
data_fig5c[,hamster2:=factor(hamster2, levels = c("2 dpi", "3 dpi", "5 dpi", "2 dpi ld", "3 dpi ld", "2 dpi hd", "3 dpi hd"))]


p_duplipathways_no4u5 = ggplot(data_fig5c[human2 !="WOS 4&5"], aes(hamster2, human2,  col = variable, size= value)) +
  geom_point( alpha = 0.1) +
  geom_point( alpha = 1, pch = 1) +
  facet_grid(celltype +term_name~., scales= "free_x",switch  = "x") + theme_minimal(base_size = 14) +
  labs(size = "Pathway enrichment factor",
       color = "Species") +
  theme(strip.text.y.left =  element_text(angle = 0, hjust = 0),
        strip.text.y.right =   element_text(angle = 0, hjust = 0),
        legend.position = "top",
        legend.text.align = 0,
        legend.box = "vertical") +
  guides(color = guide_legend(override.aes = list(size=5))) +
  scale_size_continuous(range =  c(1.5,11)) +
  scale_color_manual(values = c("#FB3640", "#0D3B66"))+
  scale_fill_manual(values = c("#FB3640", "#0D3B66"))

p_duplipathways_no4u5


#S5B
p_duplipathways_only4u5 = ggplot(data_fig5c[human2 =="WOS 4&5"], aes(hamster2, human2,  col = variable, size= value)) +
  geom_point( alpha = 0.1) +
  geom_point( alpha = 1, pch = 1) +
  facet_grid(celltype +term_name~., scales= "free_x",switch  = "x") + theme_minimal(base_size = 14) +
  labs(size = "Pathway enrichment factor",
       color = "Species") +
  theme(strip.text.y.left =  element_text(angle = 0, hjust = 0),
        strip.text.y.right =   element_text(angle = 0, hjust = 0),
        legend.position = "top",
        legend.text.align = 0,
        legend.box = "vertical") +
  guides(color = guide_legend(override.aes = list(size=5))) +
  scale_size_continuous(range =  c(1.5,11)) +
  scale_color_manual(values = c("#FB3640", "#0D3B66"))+
  scale_fill_manual(values = c("#FB3640", "#0D3B66"))

p_duplipathways_only4u5
