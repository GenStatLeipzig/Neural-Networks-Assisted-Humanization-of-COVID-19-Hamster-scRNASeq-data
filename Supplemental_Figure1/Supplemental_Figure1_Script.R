#Supplemental Figure 1

require(toolboxH)
require(here)
require(patchwork)
require(ggplot2)
require(UpSetR)

#S1
plotdata = fread("Data/h224_3_overlap_genes_all_datasets_plotdat.txt.gz")
upset(plotdata,nsets = ncol(plotdata)-1,order.by =c("freq" ,"degree"), decreasing = c(TRUE, TRUE), set_size.scale_max =35000, text.scale = 1.3, set_size.show = TRUE,mainbar.y.label = "Overlap Genes", nintersects = ncol(plotdata))
