
library(plyr)
library(dplyr)
library(tibble)
library(ggplot2)
library(stringr)
library(Seurat)
library(HDF5Array)
library(rhdf5)
library(hdf5r)
library(DelayedArray)
library(DelayedMatrixStats)
require(gridExtra)
library(tidyr)
library("dendextend")
library(DESeq2)
source("~/Documents/Largescale-data/notes and scripts/smooth_DimPlot.R")
library(akima)
library(pheatmap)
#see http://www.cookbook-r.com/Graphs/ Plotting_means_and_error_bars_(ggplot2)
#source("~/Documents/Largescale-data/notes and scripts/summarySE.R")
source("/fast/AG_Landthaler/scripts/summarySE.R")
library(cowplot)
library(SeuratObject)
library(lme4)
library(ComplexHeatmap)
library(patchwork)
library(forcats)
library(DoubletFinder)




#For more detailed thresholding, start with combined, non-integrated object with very low threshold (nGene=250)
#This is done in folder /fast/AG_Landthaler/tmp/ML_EW_pr/pr_new_blood
pr_new <- readRDS("../pr_seu_blood_new_combined_250.rds")

pr_new = NormalizeData(pr_new)
pr_new = FindVariableFeatures(pr_new, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pr_new)
pr_new <- ScaleData(pr_new, features = all.genes)
pr_new <- RunPCA(pr_new, features = VariableFeatures(object = pr_new))
pdf("pr_seu_blood_new_combined_250_elbow.pdf")
ElbowPlot(pr_new)
dev.off()
pr_new <- RunUMAP(pr_new, dims = 1:19)
pdf("pr_seu_blood_new_combined_250_UMAPPlot.pdf")
UMAPPlot(object=pr_new)
dev.off()
pr_new <- FindNeighbors(pr_new, dims = 1:19)
pr_new <- FindClusters(pr_new, resolution = 0.9)
pdf("pr_seu_blood_new_combined_250_clusters.pdf")
DimPlot(pr_new, reduction = "umap", label=TRUE)
dev.off()

pr_new@meta.data = cbind(pr_new@meta.data, pr_new@reductions$umap@cell.embeddings)

pr_new@meta.data$timepoint <- gsub("pr_(.*D[0-9])_Z([0-9])_B","\\1",pr_new@meta.data$orig.ident)
pr_new@meta.data$hamster <- gsub("pr_(.*D[0-9])_Z([0-9])_B","Ha\\2",pr_new@meta.data$orig.ident)


ggplot()+geom_violin(data=pr_new@meta.data, aes(x=seurat_clusters, y=log2(nCount_RNA)))
ggsave("pr_blood_new_violin_nUMI_log2_in_clusters.pdf", width=12, height=7)
ggplot()+geom_violin(data=pr_new@meta.data, aes(x=seurat_clusters, y=log2(nFeature_RNA)))
ggsave("pr_blood_new_violin_nGene_log2_in_clusters.pdf", width=12, height=7)
ggplot()+geom_bar(data=pr_new@meta.data, aes(x=seurat_clusters, y=log2(nFeature_RNA)), stat = "summary", fun = "mean")
ggsave("pr_blood_new_barplot_log2_nGene_mean_clusters.pdf", width=12, height=7)
ggplot()+geom_bar(data=pr_new@meta.data, aes(x=seurat_clusters, y=log2(nFeature_RNA)), stat = "summary", fun = "median")
ggsave("pr_blood_new_barplot_log2_nGene_median_clusters.pdf", width=12, height=7)
ggplot()+geom_bar(data=pr_new@meta.data, aes(x=seurat_clusters, y=log2(nCount_RNA)), stat = "summary", fun = "mean")
ggsave("pr_blood_new_barplot_log2_nUMI_mean_clusters.pdf", width=12, height=7)
ggplot()+geom_bar(data=pr_new@meta.data, aes(x=seurat_clusters, y=log2(nCount_RNA)), stat = "summary", fun = "median")
ggsave("pr_blood_new_barplot_log2_nUMI_median_clusters.pdf", width=12, height=7)
ggplot()+geom_bar(data=pr_new@meta.data, aes(x=seurat_clusters, y=nFeature_RNA), stat = "summary", fun = "mean")
ggsave("pr_blood_new_barplot_nGene_mean_clusters.pdf", width=12, height=7)
ggplot()+geom_bar(data=pr_new@meta.data, aes(x=seurat_clusters, y=nFeature_RNA), stat = "summary", fun = "median")
ggsave("pr_blood_new_barplot_nGene_median_clusters.pdf", width=12, height=7)
ggplot()+geom_bar(data=pr_new@meta.data, aes(x=seurat_clusters, y=nCount_RNA), stat = "summary", fun = "mean")
ggsave("pr_blood_new_barplot_nUMI_mean_clusters.pdf", width=12, height=7)
ggplot()+geom_bar(data=pr_new@meta.data, aes(x=seurat_clusters, y=nCount_RNA), stat = "summary", fun = "median")
ggsave("pr_blood_new_barplot_nUMI_median_clusters.pdf", width=12, height=7)




###############
#cell type annotation


markers <- c("Adgre1", "Ccr2", "Cx3cr1", "Tcf4", "Irf8", "Flt3", "Retn", "Jchain", "Slc3a2", "Slc7a5", "Camp", "S100a8", "Cxcr2", "Gng11", "Ppbp", "Cd79b", "Ms4a1", "Nkg7", "Gzma", "Il7r", "Cd8a", "Cd4", "Cd3e", "Cd27", "Tnfrsf17", "Prdm1", "Xbp1", "Irf4", "Sec11c", "Fkbp11", "Mki67", "Top2a")
for (gene in markers) {
  df <- FetchData(pr_new, gene)
  gene <- str_replace_all(gene, "-", "_")
  colnames(df) <- gene
  df = cbind(df, pr_new@reductions$umap@cell.embeddings)
  ggplot()+geom_point(data=subset(df, eval(parse(text = gene))<=0), aes(x=UMAP_1, y=UMAP_2), colour="grey90", shape=16, size=0.5)+
    geom_point(data=subset(df, eval(parse(text = gene)) > 0), aes(x=UMAP_1, y=UMAP_2, colour=eval(parse(text = gene))), shape=16, size=0.5)+
    scale_colour_gradientn(colours=c("#00004C", "#0C005C", "#1D006B", "#31007A", "#49008A", "#660099", "#8700A8", "#AB00B8", "#C700BA", "#D600AB", "#E60099", "#F50083", "#FF0569", "#FD1754", "#FF2441", "#FF3333"), guide="colourbar")+
    theme_bw()+
    theme(axis.ticks=element_blank(), panel.grid.minor=element_blank())+
    coord_fixed(ratio=1)+
    ggtitle(gene)
  ggsave(paste("pr_seu_blood_new_combined_250", gene, "pdf", sep="."), useDingbats=FALSE)
}

avg <- AverageExpression(pr_new, assays=c("RNA"), features = markers, return.seurat = T, slot="data") 
avg2 <- AverageExpression(pr_new, assays=c("RNA"), features = markers, return.seurat = F, slot="data")
avg2 <- as.matrix(avg2$RNA)
hc <- hclust(dist(avg2))
gene.order <- row.names(avg2)[order.hclust(hc)]
DoHeatmap(avg, size=5, features=gene.order)+scale_fill_gradient2(low="midnightblue", mid= "white", high="red", na.value = "white")+NoLegend()
ggsave("pr_seu_blood_new_combined_250_cluster_heatmap.pdf", width=12, height=7)

cluster.markers <- FindAllMarkers(pr_new, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
cluster.markers %>% group_by(cluster) %>% slice_max(n = 2, order_by = avg_log2FC)
write.table(cluster.markers, "./pr_seu_blood_new_combined_250_cluster_markers.txt", sep="\t", quote=FALSE)
top10 <- cluster.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
pdf("pr_seu_blood_new_combined_250_cluster_markers_heatmap.pdf", width=32, height=24)
DoHeatmap(pr_new, features = top10$gene) + NoLegend()
dev.off()

#Include doublet information by reading int the DoubSing file and attach to meta data
pr_new@meta.data$cell_id <- rownames(pr_new@meta.data)
DoubSing <- read.table("../DoubSing_pr_blood.txt", header=TRUE)
pr_new@meta.data <- merge(pr_new@meta.data, DoubSing, by="cell_id")
rownames(pr_new@meta.data) <- pr_new@meta.data$cell_id


Idents(pr_new) <- pr_new@meta.data$seurat_clusters
pr_new <- RenameIdents(pr_new, 
                       '0'='Neutrophil',
                       '1'='B',
                       '2'='Immature neutrophil',
                       '3'='T',
                       '4'='Neutrophil',
                       '5'='B',
                       '6'='B',
                       '7'='B',
                       '8'='B',
                       '9'='mixed1',
                       '10'='mixed2',
                       '11'='mixed3',
                       '12'='mixed4',
                       '13'='Non-classical monocyte',
                       '14'='Classical monocyte',
                       '15'='Classical monocyte',
                       '16'='Neutrophil',
                       '17'='B',
                       '18'='T',
                       '19'='B',
                       '20'='mixed5',
                       '21'='T',
                       '22'='NK',
                       '23'='Neutrophil',
                       '24'='Classical monocyte',
                       '25'='Platelet',
                       '26'='Neutrophil',
                       '27'='Classical monocyte',
                       '28'='mixed6',
                       '29'='Neutrophil',
                       '30'='mDC',
                       '31'='mixed7',
                       '32'='mixed8',
                       '33'='mixed9',
                       '34'='mDC',
                       '35'='mixed10',
                       '36'='B',
                       '37'='mixed11')
pr_new@meta.data$celltype <- Idents(pr_new)



#save at this point
saveRDS(pr_new, "./pr_seu_blood_new_combined_250_ann.rds")




#Order cell types
the_celltypes = c("Classical monocyte",
                  "Non-classical monocyte",
                  "Neutrophil",
                  "Immature neutrophil",
                  "mDC",
                  "pDC",
                  "NK",
                  "T",
                  "Activated T",
                  "B",
                  "Platelet", "mixed1",
                  "mixed2",
                  "mixed3",
                  "mixed4",
                  "mixed5",
                  "mixed6",
                  "mixed7",
                  "mixed8", "mixed9", "mixed10", "mixed11")


celltypecolors<- c("T" = "#368F8B", "Activated T"= "#5CC1BC", "B" = "#62C370", "Classical monocyte" = "#B7245C",
                   "Non-classical monocyte" ="#3E2F5B", "Neutrophil" = "#0081AF",
                   "Immature neutrophil" = "#00ABE7", "NK"= "#246A73", "pDC" = "#7C6A0A",
                   "mDC" = "#4F6D7A", "Platelet" = "#832D52", "unknown" = "#CAD2C5", "mixed1" = "#CAD2C5",
                   "mixed2" = "#CAD2C5",
                   "mixed3" = "#CAD2C5",
                   "mixed4" = "#CAD2C5",
                   "mixed5" = "#CAD2C5",
                   "mixed6" = "#CAD2C5",
                   "mixed7" = "#CAD2C5",
                   "mixed8" = "#CAD2C5",
                   "mixed9" = "#CAD2C5", "mixed10" = "#CAD2C5", "mixed11" = "#CAD2C5")


legendcolors <- c("D0" ="gray80", "hd_D2"="gray65", "hd_D3"="gray50", "ld_D2"="gray35", "ld_D3"="gray20")

expr <- list()
for (ct in names(celltypecolors)) {
  cbright = celltypecolors[[ct]]
  r=(col2rgb(cbright)-40)[[1]]
  r=ifelse(r<0,0,r)
  g=(col2rgb(cbright)-40)[[2]]
  g=ifelse(g<0,0,g)
  b=(col2rgb(cbright)-40)[[3]]
  b=ifelse(b<0,0,b)
  cdark = rgb(r, g, b, maxColorValue = 255)
  the_scale <- scales::seq_gradient_pal(cbright, cdark, "Lab")(seq(0,1,length.out=5))
  expr[[paste(ct, "D0", sep="_")]] <- the_scale[[1]]
  expr[[paste(ct, "hd_D2", sep="_")]] <- the_scale[[2]]
  expr[[paste(ct, "hd_D3", sep="_")]] <- the_scale[[3]]
  expr[[paste(ct, "ld_D2",  sep="_")]] <- the_scale[[4]]
  expr[[paste(ct, "ld_D3",  sep="_")]] <- the_scale[[5]]
}
the_colors = unlist(expr)


#UMAP with cell types
means <- pr_new@meta.data %>%
  group_by(celltype) %>%
  dplyr::summarise(mean_U1 = mean(UMAP_1), mean_U2 = mean(UMAP_2))

ggplot()+
  geom_point(data=pr_new@meta.data, aes(x=UMAP_1, y=UMAP_2, colour=celltype), size=0.5, shape=16)+
  geom_text(data=means, aes(x=mean_U1, y=mean_U2, label=celltype), size=3)+
  coord_fixed(ratio=1)+
  ggtitle("celltype")+
  theme_bw()+
  theme(axis.ticks=element_blank(), panel.grid.minor=element_blank())+
  scale_colour_manual(values = celltypecolors)
ggsave("pr_seu_blood_new_combined_250_celltypes.pdf", useDingbats=FALSE, width=14, height=7)




#doublets by celltype and timepoint
df1 <- cbind.data.frame(pr_new@meta.data$DoubSing,
                        pr_new@meta.data$celltype,
                        pr_new@meta.data$timepoint,
                        pr_new@meta.data$hamster)
colnames(df1) <- c("DoubSing", "celltype", "timepoint", "hamster")

a = df1 %>% group_by(celltype, timepoint, hamster) %>% tally(name="tot") 
b = df1 %>% filter(DoubSing == "Doublet") %>% group_by(celltype, timepoint, hamster) %>% tally(name="pos") 
c =
  left_join(a , b , by = c('celltype', 'timepoint', 'hamster')) %>% 
  replace(., is.na(.), 0) %>%
  mutate(fraction = pos / tot)
#write.table(c, paste("days", gene, "tsv", sep="."), sep="\t", quote=FALSE, row.names = FALSE)
tgc <- summarySE(as.data.frame(c), measurevar="fraction", groupvars=c("celltype", "timepoint"))
tgc = tgc %>% mutate(celltype = forcats::fct_relevel(celltype, the_celltypes))
ggplot(tgc, aes(x=celltype, y=fraction, fill=paste(celltype, timepoint, sep="_")))+
  geom_bar(aes(colour=timepoint), position=position_dodge(.8), stat="identity", size=0.2)+
  geom_point(data=c, aes(x=celltype, y=fraction, fill=paste(celltype, timepoint, sep="_")), position=position_jitterdodge(jitter.width = 0.1, dodge.width = 0.8), size=0.2)+
  geom_errorbar(aes(ymin=fraction-sd, ymax=fraction+sd), width=.2, position=position_dodge(.8), size=0.25)+
  theme_bw()+
  theme(axis.ticks=element_blank(), panel.grid.minor=element_blank())+
  ggtitle(paste("DoubSing per timepoint", sep=" "))+
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=.5))+
  scale_fill_manual(values=the_colors)+
  scale_color_manual(values = legendcolors)+
  guides(fill=FALSE, color = guide_legend(override.aes = list(size=5, fill="white")))
ggsave("pr_seu_blood_new_combined_250_DoubSing.pdf", useDingbats=FALSE, width=14, height=7)



ggplot()+geom_violin(data=pr_new@meta.data %>% mutate(celltype = forcats::fct_relevel(celltype, the_celltypes)), aes(x=celltype, y=log2(nCount_RNA), fill=celltype))+scale_fill_manual(values = celltypecolors)+theme(axis.text.x=element_text(angle=90,hjust=1,vjust=.5, size=10), axis.text.y=element_text(size=10))
ggsave("pr_seu_blood_new_combined_250_violin_nUMI_log2_in_celltype.pdf", width=12, height=7)
ggplot()+geom_violin(data=pr_new@meta.data %>% mutate(celltype = forcats::fct_relevel(celltype, the_celltypes)), aes(x=celltype, y=log2(nFeature_RNA), fill=celltype))+scale_fill_manual(values = celltypecolors)+theme(axis.text.x=element_text(angle=90,hjust=1,vjust=.5, size=10), axis.text.y=element_text(size=10))
ggsave("pr_seu_blood_new_combined_250_violin_nGene_log2_in_celltype.pdf", width=12, height=7)









#############
#final filtering

#Most cell types get normal filtering, mixed without outliers, endothelial and TNK with a subset differentiation
#mixed1, 2, 3, 10 are low quality

expr <- list()
for (ct in c("Neutrophil", 
             "B", "Immature neutrophil", "T",  
             "Non-classical monocyte", "Classical monocyte", 
             "NK", "Platelet", "mDC")) {
  df <- subset(pr_new@meta.data, celltype==ct)
  expr[[ct]] <- as.data.frame(subset(df, nCount_RNA>median(df$nCount_RNA) & nCount_RNA < quantile(df$nCount_RNA)[[4]]+3*IQR(df$nCount_RNA))$cell_id)
  colnames(expr[[ct]]) <- "cell_id"
}  

for (ct in c("mixed4", "mixed5", "mixed6",  "mixed7", "mixed8", "mixed9", "mixed11")) {
  df <- subset(pr_new@meta.data, celltype==ct)
  expr[[ct]] <- as.data.frame(subset(df, nCount_RNA>median(subset(pr_new@meta.data, celltype=="Neutrophils")$nCount_RNA))$cell_id)
  colnames(expr[[ct]]) <- "cell_id"
}  



cells_to_keep <- do.call(rbind,expr)
write.table(cells_to_keep, "../pr_blood_cells_to_keep.txt", sep="\t", quote=FALSE)

write.table(pr_new@meta.data, "./pr_seu_blood_new_combined_250_metadata.txt", sep="\t", quote=FALSE)

#######################
#Run hamster_merging_blood_cellstokeep.R, which does the thresholding using cells_to_keep and then integrates

#read integrated object from cluster
pr.int <- readRDS("/Volumes/fast/AG_Landthaler/tmp/ML_EW_pr/pr_seu_blood_new_combined_integrated.rds")

pr.int@meta.data = cbind(pr.int@meta.data, pr.int@reductions$umap@cell.embeddings)
pr.int@meta.data$timepoint <- gsub("pr_(.*D[0-9])_Z([0-9])_B","\\1",pr.int@meta.data$orig.ident)
pr.int@meta.data$hamster <- gsub("pr_(.*D[0-9])_Z([0-9])_B","Ha\\2",pr.int@meta.data$orig.ident)

DefaultAssay(pr.int) <- 'RNA'
SCoV2_rawcounts <- FetchData(pr.int, grep("SCoV2", pr.int@assays$RNA@counts@Dimnames[[1]], value="TRUE"), slot="counts")
pr.int@meta.data$SCoV2_sum <- rowSums(SCoV2_rawcounts)
SCoV2_rawcounts$SCoV2sum <- rowSums(SCoV2_rawcounts)
pr.int@meta.data$SCoV2_load <- SCoV2_rawcounts$SCoV2sum/pr.int@meta.data$nCount_RNA*100
DefaultAssay(pr.int) <- 'SCT'

DefaultAssay(pr.int) <- 'integrated'
pr.int <- FindClusters(pr.int, resolution = 0.6)
DefaultAssay(pr.int) <- 'SCT'


ggplot()+
  geom_point(data=pr.int@meta.data, aes(x=UMAP_1, y=UMAP_2, colour=timepoint), size=0.5, shape=16, alpha=0.5)+
  ggtitle("timepoint")+
  theme_bw()+
  theme(axis.ticks=element_blank(), panel.grid.minor=element_blank())
ggsave("pr_seu_blood_new_integrated_timepoint.pdf", useDingbats=FALSE)
ggplot()+
  geom_point(data=pr.int@meta.data, aes(x=UMAP_1, y=UMAP_2, colour=hamster), size=0.2, shape=16, alpha=0.5)+
  ggtitle("hamster")+
  theme_bw()+
  theme(axis.ticks=element_blank(), panel.grid.minor=element_blank())
ggsave("pr_seu_blood_new_integrated_hamsters.pdf", useDingbats=FALSE)
UMAPPlot(pr.int, label=TRUE)+
  theme_bw()+
  theme(axis.ticks=element_blank(), panel.grid.minor=element_blank())+
  ggtitle("Clusters")
ggsave("pr_seu_blood_new_integrated_clusters.pdf", useDingbats=FALSE)

ggplot()+geom_point(data=subset(pr.int@meta.data, SCoV2_load<=0), aes(x=UMAP_1, y=UMAP_2), colour="grey90", size=0.5, shape = 16)+
  geom_point(data=subset(pr.int@meta.data, SCoV2_load>0), aes(x=UMAP_1, y=UMAP_2, colour=log10(SCoV2_load)), size=0.5, shape = 16)+
  scale_colour_gradientn(colours=c("#00004C", "#0C005C", "#1D006B", "#31007A", "#49008A", "#660099", "#8700A8", "#AB00B8", "#C700BA", "#D600AB", "#E60099", "#F50083", "#FF0569", "#FD1754", "#FF2441", "#FF3333"), guide="colourbar")+
  theme_bw()+
  theme(axis.ticks=element_blank(), panel.grid.minor=element_blank())+
  ggtitle("log10 SCoV2 load")
ggsave(paste("pr_seu_blood_new_integrated_log10SCoV2_load", "pdf", sep="."), useDingbats=FALSE)

unfilt_metadata <- read.table("/Volumes/fast/AG_Landthaler/tmp/ML_EW_pr/pr_blood_new/pr_seu_blood_new_combined_250_metadata.txt", sep="\t")
pr.int@meta.data$cell_id <- rownames(pr.int@meta.data)
df <- dplyr::left_join(pr.int@meta.data, unfilt_metadata, by="cell_id")
pr.int@meta.data$DoubSing <- df$DoubSing

#cell type annotation
#LOC101833790 = Ccr5

markers <- c("Ccl5", "Adgre1", "Ccr2", "Cx3cr1", "Tcf4", "Irf8", "Flt3", "Retn", "Jchain", "Slc3a2", "Slc7a5", "Camp", "S100a8", "Cxcr2", "Gng11", "Ppbp", "Cd79b", "Ms4a1", "Nkg7", "Gzma", "Il7r", "Cd8a", "Cd4", "Cd3e", "Cd27", "Tnfrsf17", "Prdm1", "Xbp1", "Irf4", "Sec11c", "Fkbp11", "Mki67", "Top2a")
for (gene in c("Ccl5")) {
  df <- FetchData(pr.int, gene)
  gene <- str_replace_all(gene, "-", "_")
  colnames(df) <- gene
  df = cbind(df, pr.int@reductions$umap@cell.embeddings)
  ggplot()+geom_point(data=subset(df, eval(parse(text = gene))<=0), aes(x=UMAP_1, y=UMAP_2), colour="grey90", shape=16, size=0.5)+
    geom_point(data=subset(df, eval(parse(text = gene)) > 0), aes(x=UMAP_1, y=UMAP_2, colour=eval(parse(text = gene))), shape=16, size=0.5)+
    scale_colour_gradientn(colours=c("#00004C", "#0C005C", "#1D006B", "#31007A", "#49008A", "#660099", "#8700A8", "#AB00B8", "#C700BA", "#D600AB", "#E60099", "#F50083", "#FF0569", "#FD1754", "#FF2441", "#FF3333"), guide="colourbar")+
    theme_bw()+
    theme(axis.ticks=element_blank(), panel.grid.minor=element_blank())+
    coord_fixed(ratio=1)+
    ggtitle(gene)
  ggsave(paste("pr_seu_blood_new_integrated", gene, "pdf", sep="."), useDingbats=FALSE)
}

avg <- AverageExpression(pr.int, assays=c("SCT"), features = markers, return.seurat = T, slot="data") 
avg2 <- AverageExpression(pr.int, assays=c("SCT"), features = markers, return.seurat = F, slot="data")
avg2 <- as.matrix(avg2$SCT)
hc <- hclust(dist(avg2))
gene.order <- row.names(avg2)[order.hclust(hc)]
DoHeatmap(avg, size=5, features=gene.order)+scale_fill_gradient2(low="midnightblue", mid= "white", high="red", na.value = "white")+NoLegend()
ggsave("pr_seu_blood_new_integrated_cluster_heatmap.pdf", width=12, height=7)


Idents(pr.int) <- pr.int@meta.data$seurat_clusters

#This is the annotation for the object with more than 250 genes
#pDC are likely mixed into the myDC cluster (24)
pr.int <- RenameIdents(pr.int, 
                       '0'='Immature neutrophil',
                       '1'='Neutrophil',
                       '2'='B',
                       '3'='B',
                       '4'='B',
                       '5'='Classical monocyte',
                       '6'='Neutrophil',
                       '7'='T',
                       '8'='T',
                       '9'='Non-classical monocyte',
                       '10'='B',
                       '11'='T',
                       '12'='NK',
                       '13'='Platelet',
                       '14'='B',
                       '15'='B',
                       '16'='mixed1',
                       '17'='Neutrophil',
                       '18'='mDC',
                       '19'='mixed2')
pr.int@meta.data$celltype <- Idents(pr.int)

saveRDS(pr.int, "/Volumes/fast/AG_Landthaler/tmp/ML_EW_pr/pr_blood_new_combined_integrated_annotated.rds")

#UMAP with cell types
means <- pr.int@meta.data %>%
  group_by(celltype) %>%
  dplyr::summarise(mean_U1 = mean(UMAP_1), mean_U2 = mean(UMAP_2)) 

ggplot()+
  geom_point(data=pr.int@meta.data, aes(x=UMAP_1, y=UMAP_2, colour=celltype), size=0.5, shape=16)+
  geom_text(data=means, aes(x=mean_U1, y=mean_U2, label=celltype), size=3)+
  coord_fixed(ratio=1)+
  ggtitle("celltype")+
  theme_bw()+
  theme(axis.ticks=element_blank(), panel.grid.minor=element_blank())+
  scale_colour_manual(values = celltypecolors)
ggsave("pr_seu_blood_new_integrated_celltypes.pdf", useDingbats=FALSE, width=14, height=7)

#percent cell type per timepoint
df1 <- cbind.data.frame(pr.int@meta.data$SCoV2_load,
                        pr.int@meta.data$celltype,
                        pr.int@meta.data$timepoint,
                        pr.int@meta.data$hamster)
colnames(df1) <- c("SCoV2_load", "celltype", "timepoint", "hamster")

a = df1 %>% group_by(timepoint, hamster) %>% tally(name="tot") 
b = df1 %>% group_by(celltype, timepoint, hamster) %>% tally(name="pos") 
c =
  left_join(a , b , by = c('timepoint', 'hamster')) %>% 
  replace(., is.na(.), 0) %>%
  mutate(fraction = pos / tot)
#write.table(c, "celltypefractions.tsv", sep="\t", quote=FALSE, row.names = FALSE)

tgc <- summarySE(as.data.frame(c), measurevar="fraction", groupvars=c("celltype", "timepoint"))
tgc = tgc %>% mutate(celltype = forcats::fct_relevel(celltype, the_celltypes))
ggplot(tgc, aes(x=celltype, y=fraction, fill=paste(celltype, timepoint, sep="_")))+
  geom_bar(aes(colour=timepoint), position=position_dodge(.9), stat="identity", size=0.2)+
  geom_point(data=c, aes(x=celltype, y=fraction, fill=paste(celltype, timepoint, sep="_")), position=position_jitterdodge(jitter.width = 0.1, dodge.width = 0.9), size=0.2)+
  geom_errorbar(aes(ymin=fraction-sd, ymax=fraction+sd), width=.2, position=position_dodge(.9), size=0.25)+
  theme_bw()+
  theme(axis.ticks=element_blank(), panel.grid.minor=element_blank())+
  ggtitle(paste("celltypes per timepoint", sep=" "))+
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=.5))+
  scale_fill_manual(values=the_colors)+
  scale_color_manual(values = legendcolors)+
  guides(fill=FALSE, color = guide_legend(override.aes = list(size=5, fill="white")))

ggsave("pr_seu_blood_new_integrated_celltypepercentage_pertimepoint.pdf", useDingbats=FALSE)

#percent virus positive by timepoint
a = df1 %>% group_by(celltype, timepoint, hamster) %>% tally(name="tot") 
b = df1 %>% filter(SCoV2_load>0) %>% group_by(celltype, timepoint, hamster) %>% tally(name="pos") 
c =
  left_join(a , b , by = c('celltype', 'timepoint', 'hamster')) %>% 
  replace(., is.na(.), 0) %>%
  mutate(fraction = pos / tot)
#write.table(c, paste("days", gene, "tsv", sep="."), sep="\t", quote=FALSE, row.names = FALSE)
tgc <- summarySE(as.data.frame(c), measurevar="fraction", groupvars=c("celltype", "timepoint"))
tgc = tgc %>% mutate(celltype = forcats::fct_relevel(celltype, the_celltypes))
ggplot(tgc, aes(x=celltype, y=fraction, fill=paste(celltype, timepoint, sep="_")))+
  geom_bar(aes(colour=timepoint), position=position_dodge(.9), stat="identity", size=0.2)+
  geom_point(data=c, aes(x=celltype, y=fraction, fill=paste(celltype, timepoint, sep="_")), position=position_jitterdodge(jitter.width = 0.1, dodge.width = 0.9), size=0.2)+
  geom_errorbar(aes(ymin=fraction-sd, ymax=fraction+sd), width=.2, position=position_dodge(.9), size=0.25)+
  theme_bw()+
  theme(axis.ticks=element_blank(), panel.grid.minor=element_blank())+
  ggtitle(paste("virus positive per timepoint", sep=" "))+
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=.5))+
  scale_fill_manual(values=the_colors)+
  scale_color_manual(values = legendcolors)+
  guides(fill=FALSE, color = guide_legend(override.aes = list(size=5, fill="white")))
ggsave("pr_seu_blood_new_integrated_barplot_viruspositivepercelltype_pertimepoint.pdf", useDingbats=FALSE)
