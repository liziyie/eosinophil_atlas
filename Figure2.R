library(Seurat)
library(ComplexHeatmap)
library(circlize)
library(tidyverse)
library(RColorBrewer)
library(scCustomize)
library(harmony)
c_ndcirco <- c("#fb942a","#d945a0","#415f96","#c8ebe8","#f2a7b8",
               "#516dd3","#f5f179","#fbc160","#5cbac8","#97c6df",
               "#f97d3c","#8457c0","#dbc5d5","#7437c7","#e93e83")

#Figure 2a-----
eos_0.8 <- readRDS("./RDS/Fig1e_sc_eosinophil.rds")
DefaultAssay(eos_0.8) <- "AbSeq"
## Not include in
`%ni%` <- Negate(`%in%`)

eos <- eos_0.8 %>% subset(Tissue %ni% "Lung")
mean_count_Ab <- aggregate(t(eos@assays$AbSeq@counts), list(eos@meta.data$Tissue), mean)
row.names(mean_count_Ab) <- mean_count_Ab$Group.1
mean_count_Ab$Group.1 <- c()
Ab_used <- colnames(mean_count_Ab)[colSums(mean_count_Ab == 0) == 0]
Ab_eos_0.8 <- subset(eos, features = Ab_used)
Ab_eos_0.8 <- NormalizeData(Ab_eos_0.8, normalization.method = "CLR")
Ab_eos_0.8 <- Ab_eos_0.8 %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA()
Ab_eos_0.8 <- RunHarmony(Ab_eos_0.8, group.by.vars = "Batch", max.iter.harmony = 20)
Ab_eos_0.8 <- FindNeighbors(Ab_eos_0.8, reduction = "harmony", dims = 1:15)
Ab_eos_0.8 <- FindClusters(Ab_eos_0.8, resolution = 0.4)
Ab_eos_0.8 <- RunUMAP(Ab_eos_0.8, reduction = "harmony", dims = 1:15)

DimPlot_scCustom(Ab_eos_0.8, group.by = "Tissue") + scale_color_manual(values = c_ndcirco)

#Figure 2b-----
library(cowplot)
plot_list <- list()
for(tissue in unique(Ab_eos_0.8@meta.data$Tissue)){
  highlight_cells <- Ab_eos_0.8@meta.data %>% filter(Tissue == tissue) %>% pull(CellName)
  plot_list[[tissue]] <- DimPlot_scCustom(Ab_eos_0.8, colors_use = "lightgrey", 
                                          cells.highlight = highlight_cells, 
                                          cols.highlight = "darkblue", 
                                          label = F, pt.size = .1, 
                                          sizes.highlight = .01) + 
    ggtitle(tissue) + theme_cowplot(font_size = 7) + theme(legend.position = "none")
}
new_order <- c("BM", "Blood", "eWAT", "BAT", "scWAT", "SI", "Colon")
plot_list <- plot_list[new_order]
p <- plot_grid(plotlist = plot_list)

#Figure 2c-----
Ab <- GetAssayData(eos_0.8, assay = "AbSeq", slot = "data")
Ab_used <- rownames(Ab)
metadata<- as.data.frame(eos_0.8@meta.data) %>% 
  dplyr::select(c("CellName","Tissue"))
metadata<-metadata[rownames(metadata) %in% colnames(Ab), ]

markers_mean_exp <- 
  aggregate(
    as.matrix(t(eos_0.8@assays$AbSeq@data[Ab_used, ])),
    list(Tissue = metadata$Tissue),
    mean)
row.names(markers_mean_exp) <- markers_mean_exp$Tissue
markers_mean_exp$Tissue <- c()

## zscore
zscore <- function(x){
  return( (x - mean(x))/sd(x) )
}

markers_plot_matrix <- apply(markers_mean_exp, 2, zscore) %>% t()
max.avg <- apply(markers_plot_matrix, 1, which.max)
#Ab_order
Ab_order <- c()
for(i in 1:ncol(markers_plot_matrix)){
  if(sum(max.avg == i) > 1){
    temp <- data.frame(
      Ab = names(sort(markers_plot_matrix[names(max.avg)[max.avg == i],i], decreasing = T)), 
      Tissue = colnames(markers_plot_matrix)[i], stringsAsFactors = F)
  }else if(sum(max.avg == i) == 1){
    temp <- data.frame(Ab = names(max.avg)[max.avg == i], 
                       Tissue = colnames(markers_plot_matrix)[i])
  }else{
    temp <- c()
  }
  Ab_order <- rbind(Ab_order, temp)
}
Ab_order <- Ab_order[nrow(Ab_order):1,]
color_used <- circlize::colorRamp2(seq(min(markers_plot_matrix), 
                                       max(markers_plot_matrix), 
                                       length = 8), 
                                   RColorBrewer::brewer.pal(9,'BuPu')[1:8])

plotdata <- as.data.frame(markers_plot_matrix[Ab_order$Ab,])

p <- Heatmap(markers_plot_matrix[Ab_order$Ab,],
             cluster_rows = FALSE, 
             cluster_columns = FALSE,
             show_row_names = T,
             column_names_gp = gpar(fontsize = 10),
             row_names_gp = gpar(fontsize = 6),
             col = color_used,
             name = "Exp")
draw(p)

#Figure 2j-----
cytof <- readRDS("./RDS/Fig2i_cytof_AvgExp.rds")
eos_0.8 <- readRDS("./RDS/Fig1e_sc_eosinophil.rds")

cytof[["RNA"]] <- as(object = cytof[["RNA"]], Class = "Assay")
cytof_exp <- cytof@assays$RNA@data
colnames(cytof_exp) <- c("Liver","Lung","Spleen","Thymus","eWAT","BAT","scWAT","Colon","SI","Blood","BM")

eos_gene <- aggregate(t(eos_0.8@assays$RNA@data), list(Tissue = eos_0.8@meta.data[,"Tissue"]), mean)
row.names(eos_gene) <- eos_gene$Tissue
eos_gene$Tissue <- c()
eos_gene <- t(eos_gene)

eos_Abseq <- aggregate(t(eos_0.8@assays$AbSeq@data), list(Tissue = eos_0.8@meta.data[colnames(eos_0.8@assays$AbSeq@data),"Tissue"]), mean)
row.names(eos_Abseq) <- eos_Abseq$Tissue
eos_Abseq$Tissue <- c()
eos_Abseq <- t(eos_Abseq)

AbSeq_used <- c("CD11c:HL3-Itgax-AMM2008-pAbO","Siglec-F-Siglecf-AMM2013-pAbO",
                "CD274-Cd274-AMM2038-pAbO","CD117:2B8-Kit-AMM2023-pAbO",
                "CD16-2-Fcgr4-AMM2126-pAbO","CD90.2:53-2.1-Thy1-AMM2035-pAbO",
                "CD9-Cd9-AMM2096-pAbO","CD45RB-Ptprc-AMM2102-pAbO",
                "CD43-Spn-AMM2041-pAbO","Ly-6G-Ly-6C-Ly6g-Ly6c-AMM2015-pAbO",
                "I-A-I-E-H2-Ab-Ad-Aq-Ed-Ek-AMM2019-pAbO","CD16-CD32-Fcgr3-Fcgr2-AMM2003-pAbO",
                "FcepsilonRIalpha")
Cytof_used <- c("CD11c","SIGLECF","CD274","CD117","CD16-2","CD90","CD9","CD45RB","CD43","Gr1","MHCII","CD16-32","FceRI")
genes_used <- c("Itgax","Siglecf","Cd274","Kit","Fcgr4","Thy1","Cd9","Ptprc","Spn","Ly6g","H2-Ab1","Fcgr3","Fcer1a")

Abseq_exp <- eos_Abseq[AbSeq_used,]
row.names(Abseq_exp) <- Cytof_used
Abseq_exp_norm <- t(apply(Abseq_exp,1,scale))
row.names(Abseq_exp_norm) <- row.names(Abseq_exp)
colnames(Abseq_exp_norm) <- colnames(Abseq_exp)

Cytof_exp <- cytof_exp[Cytof_used,colnames(Abseq_exp)]
Cytof_exp_norm <- t(apply(Cytof_exp,1,scale))
row.names(Cytof_exp_norm) <- row.names(Cytof_exp)
colnames(Cytof_exp_norm) <- colnames(Cytof_exp)

Gene_exp <- eos_gene[genes_used,colnames(Abseq_exp)]
Gene_exp_norm <- t(apply(Gene_exp,1,scale))
row.names(Gene_exp_norm) <- row.names(Gene_exp)
colnames(Gene_exp_norm) <- colnames(Gene_exp)

colnames(Abseq_exp_norm) <- paste0("Abseq_",colnames(Abseq_exp_norm))
colnames(Cytof_exp_norm) <- paste0("Cytof_",colnames(Cytof_exp_norm))
colnames(Gene_exp_norm) <- paste0("Gene_",colnames(Gene_exp_norm))

exp_norm_merged <- cbind(Abseq_exp_norm, Cytof_exp_norm, Gene_exp_norm)
exp_norm_merged <- exp_norm_merged[, paste0(c("Gene_","Abseq_","Cytof_"),rep(colnames(Abseq_exp), each = 3))]

annotation_col <- data.frame(stringr::str_split_fixed(colnames(exp_norm_merged),"_",2)) %>% setNames(c("Source","Tissue"))
color_used <- circlize::colorRamp2(seq(min(exp_norm_merged), 
                                       max(exp_norm_merged), length = 8),
                                   RColorBrewer::brewer.pal(9,'BuPu')[1:8])
ha <- HeatmapAnnotation(df = annotation_col,
                        col = list(Source = c("Gene" = "#90b2d8","Abseq" = "#98c8ba","Cytof" = "#fbc160"),
                                   Tissue = c("BM" = "#8bbdb4", "Blood" = "#92b1d7", "eWAT" = "#d8a1a3", "BAT" = "#e6bf9f", "scWAT" = "#ebd0d1", "SI" = "#a0a0c7", "Colon" = "#957ba6")))
p <- Heatmap(exp_norm_merged,
             cluster_rows = TRUE, 
             cluster_columns = FALSE,
             show_row_names = T,
             row_title_gp = gpar(fontsize = 3),
             column_title_gp = gpar(fontsize = 3),
             row_names_gp = gpar(fontsize = 7),
             column_names_gp = gpar(fontsize = 7),
             col = color_used,
             name = "Exp",
             bottom_annotation = ha)