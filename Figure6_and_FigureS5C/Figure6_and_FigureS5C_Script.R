#Figure 6 and Supplemental Figure 5C

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

#Figure6
data_fig6 = fread("Data/h265_8_plotdata_top5_genes_same_direction_no4u5.txt.gz")
allcelltypes = unique(data_fig6$celltype)

data_fig6[,hamster_contrast1:=factor(hamster_contrast1, levels = c("d2", "d3", "d5", "e14", "ld_D2", "ld_D3", "hd_D2", "hd_D3"))]
data_fig6[,human_contrast1:=factor(human_contrast1)]
dput(data_fig6$hamster_contrast1 %>% unique())

p_similgeneslist_no4u5 = lapply(allcelltypes, function(mycelltype) {
  # mycelltype = allcelltypes[2]
  myplotdat = data_fig6[celltype ==mycelltype & human_contrast1 != "4u5"  & topgen_samedir_nr_no4u5  <=10]
  
  if(nrow(myplotdat)==0) {message("No overlapping genes found for ",mycelltype );return(NULL)}
  myplot = ggplot(myplotdat, aes(hamster_contrast1 , human_contrast1 ,  col = Species, pch = direction , size= abslogFC)) +
    geom_point( alpha = 0.6) +
    facet_grid(features~., switch  = "x") + theme_minimal(base_size = 18) +
    labs(size = "abs. LogFoldChange",
         color = "Species") +
    theme(strip.text.y.left =  element_text(angle = 0, hjust = 0),
          strip.text.y.right =   element_text(angle = 0, hjust = 0),
          legend.position = "top",
          legend.text.align = 0,
          legend.box = "horizontal") +
    scale_shape_manual(values = c(6,2),drop=FALSE)+
    scale_size_continuous(range =  c(1.5,11), limits = c(0.5,11)) +
    scale_color_manual(values = c("#FB3640", "#0D3B66"))+
    scale_fill_manual(values = c("#FB3640", "#0D3B66")) +
    guides(color = guide_legend(override.aes = list(size=5), ncol = 1),
           pch = guide_legend(override.aes = list(size=5), ncol = 2),
           size = guide_legend(ncol = 1, override.aes = list(pch = 2))
    ) +
    scale_x_discrete(guide = guide_axis(angle = -45), drop=FALSE) +
    ggtitle(mycelltype)
  myplot
}
)
p_similgeneslist_no4u5 = p_similgeneslist_no4u5[sapply(p_similgeneslist_no4u5, is.null)==F]
plot_fig6 = patchwork::wrap_plots(p_similgeneslist_no4u5) + plot_layout(guides = "collect",ncol = 3)
plot_fig6


#Figure S5C
p_similgeneslist_just4u5 = lapply(allcelltypes, function(mycelltype) {
  # mycelltype = allcelltypes[6]
  myplotdat = data_fig6[celltype ==mycelltype & human_contrast1 == "4u5"  &topgen_samedir_nr   <=10]
  levels(myplotdat$human_contrast)
  
  if(nrow(myplotdat)==0) {message("No overlapping genes found for ",mycelltype );return(NULL)}
  myplot = ggplot(myplotdat, aes(hamster_contrast1 , human_contrast1,  col = Species, pch = direction , size= abslogFC)) +
    geom_point( alpha = 0.6) +
    facet_grid(features~., scales= "free_y",switch  = "x") + theme_minimal(base_size = 18) +
    labs(size = "abs. LogFoldChange",
         color = "Species") +
    theme(strip.text.y.left =  element_text(angle = 0, hjust = 0),
          strip.text.y.right =   element_text(angle = 0, hjust = 0),
          legend.position = "top",
          legend.text.align = 0,
          legend.box = "horizontal") +
    scale_shape_manual(values = c(6,2),drop=FALSE)+
    scale_size_continuous(range =  c(1.5,11), limits = c(0.5,11)) +
    scale_color_manual(values = c("#FB3640", "#0D3B66"))+
    scale_fill_manual(values = c("#FB3640", "#0D3B66")) +
    guides(color = guide_legend(override.aes = list(size=5), ncol = 1),
           pch = guide_legend(override.aes = list(size=5), ncol = 2),
           size = guide_legend(ncol = 1, override.aes = list(pch = 2))
    ) +
    scale_x_discrete(guide = guide_axis(angle = -45), drop=FALSE)+
    scale_y_discrete(drop=FALSE) +
    ylab("WOS")+
    xlab("")+
    
    ggtitle(mycelltype)
  myplot
}
)
p_similgeneslist_just4u5 = p_similgeneslist_just4u5[sapply(p_similgeneslist_just4u5, is.null)==F]
plot_supplfig5C = patchwork::wrap_plots(p_similgeneslist_just4u5) + plot_layout(guides = "collect",ncol = 3)
plot_supplfig5C
