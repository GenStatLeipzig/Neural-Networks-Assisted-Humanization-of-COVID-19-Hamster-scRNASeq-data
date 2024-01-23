
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
library(dendextend)
library(DESeq2)
source("~/Documents/Largescale-data/notes and scripts/smooth_DimPlot.R")
library(akima)
library(pheatmap)
#see http://www.cookbook-r.com/Graphs/ Plotting_means_and_error_bars_(ggplot2)
source("~/Documents/Largescale-data/notes and scripts/summarySE.R")
source("/fast/AG_Landthaler/scripts/summarySE.R")
library(cowplot)
library(SeuratObject)
library(lme4)
library(ComplexHeatmap)
library(patchwork)
library(forcats)
library(DoubletFinder)
library(dendextend)



#For more detailed thresholding, start with combined, non-integrated object with very low threshold (nGene=250)
#This is done in folder /fast/AG_Landthaler/tmp/ML_EW_096/ma_new_blood
ma_new <- readRDS("../ma_seu_blood_new_combined_250.rds")

ma_new = NormalizeData(ma_new)
ma_new = FindVariableFeatures(ma_new, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(ma_new)
ma_new <- ScaleData(ma_new, features = all.genes)
ma_new <- RunPCA(ma_new, features = VariableFeatures(object = ma_new))
pdf("ma_seu_blood_new_combined_250_elbow.pdf")
ElbowPlot(ma_new)
dev.off()
ma_new <- RunUMAP(ma_new, dims = 1:19)
pdf("ma_seu_blood_new_combined_250_UMAPPlot.pdf")
UMAPPlot(object=ma_new)
dev.off()
ma_new <- FindNeighbors(ma_new, dims = 1:19)
ma_new <- FindClusters(ma_new, resolution = 0.9)
pdf("ma_seu_blood_new_combined_250_clusters.pdf")
DimPlot(ma_new, reduction = "umap", label=TRUE)
dev.off()

ma_new@meta.data = cbind(ma_new@meta.data, ma_new@reductions$umap@cell.embeddings)

ma_new@meta.data$timepoint <- gsub("ma_([de0-9]*)_blood_([0-9])","\\1",ma_new@meta.data$orig.ident)
ma_new@meta.data$hamster <- gsub("ma_([de0-9]*)_blood_([0-9])","Ha\\2",ma_new@meta.data$orig.ident)


ggplot()+geom_violin(data=ma_new@meta.data, aes(x=seurat_clusters, y=log2(nCount_RNA)))
ggsave("ma_blood_new_violin_nUMI_log2_in_clusters.pdf", width=12, height=7)
ggplot()+geom_violin(data=ma_new@meta.data, aes(x=seurat_clusters, y=log2(nFeature_RNA)))
ggsave("ma_blood_new_violin_nGene_log2_in_clusters.pdf", width=12, height=7)
ggplot()+geom_bar(data=ma_new@meta.data, aes(x=seurat_clusters, y=log2(nFeature_RNA)), stat = "summary", fun = "mean")
ggsave("ma_blood_new_barplot_log2_nGene_mean_clusters.pdf", width=12, height=7)
ggplot()+geom_bar(data=ma_new@meta.data, aes(x=seurat_clusters, y=log2(nFeature_RNA)), stat = "summary", fun = "median")
ggsave("ma_blood_new_barplot_log2_nGene_median_clusters.pdf", width=12, height=7)
ggplot()+geom_bar(data=ma_new@meta.data, aes(x=seurat_clusters, y=log2(nCount_RNA)), stat = "summary", fun = "mean")
ggsave("ma_blood_new_barplot_log2_nUMI_mean_clusters.pdf", width=12, height=7)
ggplot()+geom_bar(data=ma_new@meta.data, aes(x=seurat_clusters, y=log2(nCount_RNA)), stat = "summary", fun = "median")
ggsave("ma_blood_new_barplot_log2_nUMI_median_clusters.pdf", width=12, height=7)
ggplot()+geom_bar(data=ma_new@meta.data, aes(x=seurat_clusters, y=nFeature_RNA), stat = "summary", fun = "mean")
ggsave("ma_blood_new_barplot_nGene_mean_clusters.pdf", width=12, height=7)
ggplot()+geom_bar(data=ma_new@meta.data, aes(x=seurat_clusters, y=nFeature_RNA), stat = "summary", fun = "median")
ggsave("ma_blood_new_barplot_nGene_median_clusters.pdf", width=12, height=7)
ggplot()+geom_bar(data=ma_new@meta.data, aes(x=seurat_clusters, y=nCount_RNA), stat = "summary", fun = "mean")
ggsave("ma_blood_new_barplot_nUMI_mean_clusters.pdf", width=12, height=7)
ggplot()+geom_bar(data=ma_new@meta.data, aes(x=seurat_clusters, y=nCount_RNA), stat = "summary", fun = "median")
ggsave("ma_blood_new_barplot_nUMI_median_clusters.pdf", width=12, height=7)




###############
#cell type annotation


markers <- c("Ltf", "Adgre1", "Ccr2", "Cx3cr1", "Tcf4", "Irf8", "Flt3", "Retn", "Jchain", "Sdc1", "Cd38", "Slc3a2", "Slc7a5", "Camp", "S100a8", "Cxcr2", "Gng11", "Ppbp", "Cd79b", "Ms4a1", "Nkg7", "Gzma", "Il7r", "Cd8a", "Cd4", "Cd3e", "Cd27", "Tnfrsf17", "Prdm1", "Xbp1", "Irf4", "Sec11c", "Fkbp11", "Mki67", "Top2a")
for (gene in markers) {
  df <- FetchData(ma_new, gene)
  gene <- str_replace_all(gene, "-", "_")
  colnames(df) <- gene
  df = cbind(df, ma_new@reductions$umap@cell.embeddings)
  ggplot()+geom_point(data=subset(df, eval(parse(text = gene))<=0), aes(x=UMAP_1, y=UMAP_2), colour="grey90", shape=16, size=0.5)+
    geom_point(data=subset(df, eval(parse(text = gene)) > 0), aes(x=UMAP_1, y=UMAP_2, colour=eval(parse(text = gene))), shape=16, size=0.5)+
    scale_colour_gradientn(colours=c("#00004Cx", "#0C005C", "#1D006B", "#31007A", "#49008A", "#660099", "#8700A8", "#AB00B8", "#C700BA", "#D600AB", "#E60099", "#F50083", "#FF0569", "#FD1754", "#FF2441", "#FF3333"), guide="colourbar")+
    theme_bw()+
    theme(axis.ticks=element_blank(), panel.grid.minor=element_blank())+
    coord_fixed(ratio=1)+
    ggtitle(gene)
  ggsave(paste("ma_seu_blood_new_combined_250", gene, "pdf", sep="."), useDingbats=FALSE)
}

avg <- AverageExpression(ma_new, assays=c("RNA"), features = markers, return.seurat = T, slot="data") 
avg2 <- AverageExpression(ma_new, assays=c("RNA"), features = markers, return.seurat = F, slot="data")
avg2 <- as.matrix(avg2$RNA)
hc <- hclust(dist(avg2))
gene.order <- row.names(avg2)[order.hclust(hc)]
DoHeatmap(avg, size=5, features=gene.order)+scale_fill_gradient2(low="midnightblue", mid= "white", high="red", na.value = "white")+NoLegend()
ggsave("ma_seu_blood_new_combined_250_cluster_heatmap.pdf", width=12, height=7)

cluster.markers <- FindAllMarkers(ma_new, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
cluster.markers %>% group_by(cluster) %>% slice_max(n = 2, order_by = avg_log2FC)
top10 <- cluster.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
pdf("ma_seu_blood_new_combined_250_cluster_markers_heatmap.pdf", width=24, height=16)
DoHeatmap(ma_new, features = top10$gene) + NoLegend()
dev.off()

#Include doublet information by reading int the DoubSing file and attach to meta data
ma_new@meta.data$cell_id <- rownames(ma_new@meta.data)
DoubSing <- read.table("../DoubSing_blood.txt", header=TRUE)
ma_new@meta.data <- merge(ma_new@meta.data, DoubSing, by="cell_id")
rownames(ma_new@meta.data) <- ma_new@meta.data$cell_id


Idents(ma_new) <- ma_new@meta.data$seurat_clusters
ma_new <- RenameIdents(ma_new, 
                       '0'='Immature neutrophil',
                       '1'='Neutrophil',
                       '2'='B',
                       '3'='T',
                       '4'='Classical monocyte',
                       '5'='Neutrophil',
                       '6'='T',
                       '7'='Platelet',
                       '8'='Platelet',
                       '9'='Neutrophil',
                       '10'='T',
                       '11'='mixed1',
                       '12'='Neutrophil',
                       '13'='Classical monocyte',
                       '14'='NK',
                       '15'='T',
                       '16'='B',
                       '17'='mixed2',
                       '18'='Non-classical monocyte',
                       '19'='mixed3',
                       '20'='Classical monocyte',
                       '21'='B',
                       '22'='Platelet',
                       '23'='B',
                       '24'='B',
                       '25'='mixed4',
                       '26'='mDC',
                       '27'='mixed5',
                       '28'='B',
                       '29'='Classical monocyte',
                       '30'='mixed6',
                       '31'='mixed7',
                       '32'='mixed8')
ma_new@meta.data$celltype <- Idents(ma_new)



#save at this point
saveRDS(ma_new, "./ma_seu_blood_new_combined_250_ann.rds")

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
                  "Platelet")


celltypecolors<- c("T" = "#368F8B", "Activated T"= "#5CC1BC", "B" = "#62C370", "Classical monocyte" = "#B7245C",
                   "Non-classical monocyte" ="#3E2F5B", "Neutrophil" = "#0081AF",
                   "Immature neutrophil" = "#00ABE7", "NK"= "#246A73", "pDC" = "#7C6A0A",
                   "mDC" = "#4F6D7A", "Platelet" = "#832D52", "unknown" = "#CAD2C5", "mixed1" = "#CAD2C5",
                   "mixed1" = "#CAD2C5",
                   "mixed2" = "#CAD2C5",
                   "mixed3" = "#CAD2C5",
                   "mixed4" = "#CAD2C5",
                   "mixed5" = "#CAD2C5",
                   "mixed6" = "#CAD2C5",
                   "mixed7" = "#CAD2C5",
                   "mixed8" = "#CAD2C5")


legendcolors <- c("d0" ="gray80", "d2"="gray65", "d3"="gray50", "d5"="gray35", "e14"="gray20")

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
  expr[[paste(ct, "d0", sep="_")]] <- the_scale[[1]]
  expr[[paste(ct, "d2", sep="_")]] <- the_scale[[2]]
  expr[[paste(ct, "d3", sep="_")]] <- the_scale[[3]]
  expr[[paste(ct, "d5",  sep="_")]] <- the_scale[[4]]
  expr[[paste(ct, "e14",  sep="_")]] <- the_scale[[5]]
}
the_colors = unlist(expr)


#UMAP with cell types
means <- ma_new@meta.data %>%
  group_by(celltype) %>%
  dplyr::summarise(mean_U1 = mean(UMAP_1), mean_U2 = mean(UMAP_2))

ggplot()+
  geom_point(data=ma_new@meta.data, aes(x=UMAP_1, y=UMAP_2, colour=celltype), size=0.5, shape=16)+
  geom_text(data=means, aes(x=mean_U1, y=mean_U2, label=celltype), size=3)+
  coord_fixed(ratio=1)+
  ggtitle("celltype")+
  theme_bw()+
  theme(axis.ticks=element_blank(), panel.grid.minor=element_blank())+
  scale_colour_manual(values = celltypecolors)
ggsave("ma_seu_blood_new_combined_250_celltypes.pdf", useDingbats=FALSE, width=14, height=7)




#doublets by celltype and treatment
df1 <- cbind.data.frame(ma_new@meta.data$DoubSing,
                        ma_new@meta.data$celltype,
                        ma_new@meta.data$timepoint,
                        ma_new@meta.data$hamster)
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
ggsave("ma_seu_blood_new_combined_250_DoubSing.pdf", useDingbats=FALSE, width=14, height=7)



ggplot()+geom_violin(data=ma_new@meta.data %>% mutate(celltype = forcats::fct_relevel(celltype, the_celltypes)), aes(x=celltype, y=log2(nCount_RNA), fill=celltype))+scale_fill_manual(values = celltypecolors)+theme(axis.text.x=element_text(angle=90,hjust=1,vjust=.5, size=10), axis.text.y=element_text(size=10))
ggsave("ma_seu_blood_new_combined_250_violin_nUMI_log2_in_celltype.pdf", width=12, height=7)
ggplot()+geom_violin(data=ma_new@meta.data %>% mutate(celltype = forcats::fct_relevel(celltype, the_celltypes)), aes(x=celltype, y=log2(nFeature_RNA), fill=celltype))+scale_fill_manual(values = celltypecolors)+theme(axis.text.x=element_text(angle=90,hjust=1,vjust=.5, size=10), axis.text.y=element_text(size=10))
ggsave("ma_seu_blood_new_combined_250_violin_nGene_log2_in_celltype.pdf", width=12, height=7)









#############
#final filtering

#Most cell types get normal filtering, mixed without outliers, endothelial and TNK with a subset differentiation
#omit loq quality mixed1-4

expr <- list()
for (ct in c("Immature neutrophil", "Neutrophil", "B", "T", "Classical monocyte", "Platelet", "NK", "mixed2", "Non-classical monocyte", "mDC")) {
  df <- subset(ma_new@meta.data, celltype==ct)
  expr[[ct]] <- as.data.frame(subset(df, nCount_RNA>median(df$nCount_RNA) & nCount_RNA < quantile(df$nCount_RNA)[[4]]+3*IQR(df$nCount_RNA))$cell_id)
  colnames(expr[[ct]]) <- "cell_id"
}  

for (ct in c("mixed5", "mixed6", "mixed7", "mixed8")) {
  df <- subset(ma_new@meta.data, celltype==ct)
  expr[[ct]] <- as.data.frame(subset(df, nCount_RNA>median(subset(ma_new@meta.data, celltype=="Neutrophils")$nCount_RNA))$cell_id)
  colnames(expr[[ct]]) <- "cell_id"
}  




cells_to_keep <- do.call(rbind,expr)
write.table(cells_to_keep, "../ma_blood_cells_to_keep.txt", sep="\t", quote=FALSE)
write.table(ma_new@meta.data, "./ma_seu_blood_new_combined_250_metadata.txt", sep="\t", quote=FALSE)

#######################
#Run hamster_merging_blood_cellstokeep.R, which does the thresholding using cells_to_keep and then integrates

#read integrated object from cluster
ma.int <- readRDS("/Volumes/fast/AG_Landthaler/tmp/ML_EW_096/seu_blood_new_combined_integrated.rds")

ma.int@meta.data = cbind(ma.int@meta.data, ma.int@reductions$umap@cell.embeddings)
ma.int@meta.data$timepoint <- gsub("ma_([de0-9]*)_blood_([0-9])","\\1",ma.int@meta.data$orig.ident)
ma.int@meta.data$hamster <- gsub("ma_([de0-9]*)_blood_([0-9])","Ha\\2",ma.int@meta.data$orig.ident)

DefaultAssay(ma.int) <- 'RNA'
SCoV2_rawcounts <- FetchData(ma.int, grep("SCoV2", ma.int@assays$RNA@counts@Dimnames[[1]], value="TRUE"), slot="counts")
ma.int@meta.data$SCoV2_sum <- rowSums(SCoV2_rawcounts)
SCoV2_rawcounts$SCoV2sum <- rowSums(SCoV2_rawcounts)
ma.int@meta.data$SCoV2_load <- SCoV2_rawcounts$SCoV2sum/ma.int@meta.data$nCount_RNA*100
DefaultAssay(ma.int) <- 'SCT'

ggplot()+
  geom_point(data=ma.int@meta.data, aes(x=UMAP_1, y=UMAP_2, colour=timepoint), size=0.5, shape=16, alpha=0.5)+
  ggtitle("treatment")+
  theme_bw()+
  theme(axis.ticks=element_blank(), panel.grid.minor=element_blank())
ggsave("ma_seu_blood_new_integrated_timepoint.pdf", useDingbats=FALSE)
ggplot()+
  geom_point(data=ma.int@meta.data, aes(x=UMAP_1, y=UMAP_2, colour=hamster), size=0.2, shape=16, alpha=0.5)+
  ggtitle("hamster")+
  theme_bw()+
  theme(axis.ticks=element_blank(), panel.grid.minor=element_blank())
ggsave("ma_seu_blood_new_integrated_hamsters.pdf", useDingbats=FALSE)
UMAPPlot(ma.int, label=TRUE)+
  theme_bw()+
  theme(axis.ticks=element_blank(), panel.grid.minor=element_blank())+
  ggtitle("Clusters")
ggsave("ma_seu_blood_new_integrated_clusters.pdf", useDingbats=FALSE)

ggplot()+geom_point(data=subset(ma.int@meta.data, SCoV2_load<=0), aes(x=UMAP_1, y=UMAP_2), colour="grey90", size=0.5, shape = 16)+
  geom_point(data=subset(ma.int@meta.data, SCoV2_load>0), aes(x=UMAP_1, y=UMAP_2, colour=log10(SCoV2_load)), size=0.5, shape = 16)+
  scale_colour_gradientn(colours=c("#00004C", "#0C005C", "#1D006B", "#31007A", "#49008A", "#660099", "#8700A8", "#AB00B8", "#C700BA", "#D600AB", "#E60099", "#F50083", "#FF0569", "#FD1754", "#FF2441", "#FF3333"), guide="colourbar")+
  theme_bw()+
  theme(axis.ticks=element_blank(), panel.grid.minor=element_blank())+
  ggtitle("log10 SCoV2 load")
ggsave(paste("ma_seu_blood_new_integrated_log10SCoV2_load", "pdf", sep="."), useDingbats=FALSE)

unfilt_metadata <- read.table("/Volumes/fast/AG_Landthaler/tmp/ML_EW_096/ma_new_blood/ma_seu_blood_new_combined_250_metadata.txt", sep="\t")
ma.int@meta.data$cell_id <- rownames(ma.int@meta.data)
df <- dplyr::left_join(ma.int@meta.data, unfilt_metadata, by="cell_id")
ma.int@meta.data$DoubSing <- df$DoubSing

#cell type annotation
#LOC101833790 = Ccr5

markers <- c("Ltf", "Adgre1", "Ccr2", "Cx3cr1", "Tcf4", "Irf8", "Flt3", "Retn", "Jchain", "Sdc1", "Cd38", "Slc3a2", "Slc7a5", "Camp", "S100a8", "Cxcr2", "Gng11", "Ppbp", "Cd79b", "Ms4a1", "Nkg7", "Gzma", "Il7r", "Cd8a", "Cd4", "Cd3e", "Cd27", "Tnfrsf17", "Prdm1", "Xbp1", "Irf4", "Sec11c", "Fkbp11", "Mki67", "Top2a")
for (gene in markers) {
  df <- FetchData(ma.int, gene)
  gene <- str_replace_all(gene, "-", "_")
  colnames(df) <- gene
  df = cbind(df, ma.int@reductions$umap@cell.embeddings)
  ggplot()+geom_point(data=subset(df, eval(parse(text = gene))<=0), aes(x=UMAP_1, y=UMAP_2), colour="grey90", shape=16, size=0.5)+
    geom_point(data=subset(df, eval(parse(text = gene)) > 0), aes(x=UMAP_1, y=UMAP_2, colour=eval(parse(text = gene))), shape=16, size=0.5)+
    scale_colour_gradientn(colours=c("#00004C", "#0C005C", "#1D006B", "#31007A", "#49008A", "#660099", "#8700A8", "#AB00B8", "#C700BA", "#D600AB", "#E60099", "#F50083", "#FF0569", "#FD1754", "#FF2441", "#FF3333"), guide="colourbar")+
    theme_bw()+
    theme(axis.ticks=element_blank(), panel.grid.minor=element_blank())+
    coord_fixed(ratio=1)+
    ggtitle(gene)
  ggsave(paste("ma_seu_blood_new_integrated", gene, "pdf", sep="."), useDingbats=FALSE)
}

avg <- AverageExpression(ma.int, assays=c("SCT"), features = markers, return.seurat = T, slot="data") 
avg2 <- AverageExpression(ma.int, assays=c("SCT"), features = markers, return.seurat = F, slot="data")
avg2 <- as.matrix(avg2$SCT)
hc <- hclust(dist(avg2))
gene.order <- row.names(avg2)[order.hclust(hc)]
DoHeatmap(avg, size=5, features=gene.order, raster=FALSE)+scale_fill_gradient2(low="midnightblue", mid= "white", high="red", na.value = "white")+NoLegend()
ggsave("ma_seu_blood_new_integrated_cluster_heatmap.pdf", width=12, height=7)

#Annotate cell types
Idents(ma.int) <- ma.int@meta.data$seurat_clusters
ma.int <- RenameIdents(ma.int, 
                       '0'='Immature neutrophil',
                       '1'='Neutrophil',
                       '2'='B',
                       '3'='T',
                       '4'='Classical monocyte',
                       '5'='Neutrophil',
                       '6'='Platelet',
                       '7'='T',
                       '8'='Platelet',
                       '9'='Neutrophil',
                       '10'='Neutrophil',
                       '11'='Classical monocyte',
                       '12'='B',
                       '13'='T',
                       '14'='NK',
                       '15'='Activated T',
                       '16'='B',
                       '17'='Non-classical monocyte',
                       '18'='Classical monocyte',
                       '19'='B',
                       '20'='Classical monocyte',
                       '21'='B',
                       '22'='mDC',
                       '23'='mixed1',
                       '24'='mixed2')
ma.int@meta.data$celltype <- Idents(ma.int)

saveRDS(ma.int, "/Volumes/fast/AG_Landthaler/tmp/ML_EW_096/seu_blood_new_combined_integrated_annotated.rds")
ma.int <- readRDS("/Volumes/fast/AG_Landthaler/tmp/ML_EW_096/seu_blood_new_combined_integrated_annotated.rds")



#UMAP with cell types
means <- ma.int@meta.data %>%
  group_by(celltype) %>%
  dplyr::summarise(mean_U1 = mean(UMAP_1), mean_U2 = mean(UMAP_2)) 

ggplot()+
  geom_point(data=ma.int@meta.data, aes(x=UMAP_1, y=UMAP_2, colour=celltype), size=0.5, shape=16)+
  geom_text(data=means, aes(x=mean_U1, y=mean_U2, label=celltype), size=3)+
  coord_fixed(ratio=1)+
  ggtitle("celltype")+
  theme_bw()+
  theme(axis.ticks=element_blank(), panel.grid.minor=element_blank())+
  scale_colour_manual(values = celltypecolors)
ggsave("ma_seu_blood_new_integrated_celltypes.pdf", useDingbats=FALSE, width=14, height=7)

#percent cell type per timepoint
df1 <- cbind.data.frame(ma.int@meta.data$SCoV2_load,
                        ma.int@meta.data$celltype,
                        ma.int@meta.data$timepoint,
                        ma.int@meta.data$hamster)
colnames(df1) <- c("SCoV2_load", "celltype", "timepoint", "hamster")

a = df1 %>% group_by(timepoint, hamster) %>% tally(name="tot") 
b = df1 %>% group_by(celltype, timepoint, hamster) %>% tally(name="pos") 
c =
  left_join(a , b , by = c('timepoint', 'hamster')) %>% 
  replace(., is.na(.), 0) %>%
  mutate(fraction = pos / tot)
write.table(c, "ma_seu_blood_new_integrated_celltypefractions.tsv", sep="\t", quote=FALSE, row.names = FALSE)

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

ggsave("ma_seu_blood_new_integrated_celltypepercentage_pertimepoint.pdf", useDingbats=FALSE)

