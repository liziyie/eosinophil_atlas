library(Seurat)
library(tidyverse)
library(cowplot)
library(scCustomize)
library(harmony)
c_ndcirco <- c("#fb942a","#d945a0","#415f96","#c8ebe8","#f2a7b8",
               "#516dd3","#f5f179","#fbc160","#5cbac8","#97c6df",
               "#f97d3c","#8457c0","#dbc5d5","#7437c7","#e93e83")

## Not include in
`%ni%` <- Negate(`%in%`)

#Figure 4b-----
Lung <- readRDS("./RDS/Fig4b_sc_lung.rds")
DimPlot_scCustom(Lung)+ scale_color_manual(values = c_ndcirco)

#Figure 4c-----
Colon <- readRDS("./RDS/Fig4b_sc_colon.rds")
DimPlot_scCustom(Colon_reclustering)+ scale_color_manual(values = c_ndcirco)

#Figure 4d-----
SI <- readRDS("./RDS/Fig4b_sc_SI.rds")
DimPlot_scCustom(SI_reclustering)+ scale_color_manual(values = c_ndcirco)

#Figure 4e-----
eos_sub <- readRDS("./RDS/Fig4e_sc_eosinophil_sub.rds")
DimPlot_scCustom(eos_sub, group.by = "Tissue", colors_use = tissue_color_panel)

plot_list <- list()
for(tissue in unique(eos_sub@meta.data$Tissue)){
  highlight_cells <- eos_sub@meta.data %>% filter(Tissue == tissue) %>% pull(CellName)
  plot_list[[tissue]] <- DimPlot_scCustom(eos_sub, colors_use = "lightgrey", cells.highlight = highlight_cells, cols.highlight = "darkblue", label = F, pt.size = .1, sizes.highlight = .01) + ggtitle(tissue) + theme_cowplot(font_size = 7) + theme(legend.position = "none")
}
p <- plot_grid(plotlist = plot_list)

#Figure 4f-----
expression_data <- eos_sub@assays$RNA@counts
cytotrace2_result <- cytotrace2(expression_data)
eos_sub@meta.data <- cbind(eos_sub@meta.data, cytotrace2_result)
FeaturePlot_scCustom(eos_sub, features = "CytoTRACE2_Relative") +
  scale_color_gradientn(colors = brewer.pal(9, "BuPu"),
                        breaks = c(0, 1),
                        labels = c("More diff.", "Less diff."))

eos.cds <- as.cell_data_set(eos_sub)
eos.cds <- cluster_cells(eos.cds, k = 200, reduction_method = "UMAP")
plot_cells(eos.cds, color_cells_by = "partition")
eos.cds@clusters[["UMAP"]]$clusters <- eos_0.8@active.ident[names(eos.cds@clusters[["UMAP"]]$clusters)]
eos.cds <- learn_graph(eos.cds, 
                       learn_graph_control = list(minimal_branch_len = 8,
                                                  euclidean_distance_ratio = .625), 
                       use_partition = TRUE, close_loop = T)
eos.cds <- order_cells(eos.cds)

myPalette <- colorRampPalette(brewer.pal(9, "BuPu"))
p <- plot_cells(
  eos.cds, 
  color_cells_by = "pseudotime",
  cell_size = 2,
  show_trajectory_graph = TRUE,
  label_cell_groups = FALSE, 
  label_leaves = FALSE, 
  label_branch_points = FALSE) +
  scale_color_gradientn(colours = myPalette(100)) +
  labs(x = "UMAP1", y = "UMAP2") +
  theme(legend.key.width = unit(0.5, "cm"),
        legend.key.height = unit(1, "cm"))

#Figure S7----
plot_cells(eos.cds, label_principal_points = TRUE)

Y61_2_Y122_plot.data <- monocle3_heatmap_1b(eos.cds, start.nodes = "Y_61", end.nodes = "Y_122", add.genes = c("Fcgr4","Clec12a","Cd22","Itgal","Csf2ra"))
write.csv(Y61_2_Y122_plot.data$pr_test_res, file = "./results/Y61_to_Y122_pr_test_res.csv")
Y61_2_Y125_plot.data <- monocle3_heatmap_1b(eos.cds, start.nodes = "Y_61", end.nodes = "Y_125", add.genes = c("Clec12a","Cd22","Itgal","Csf2ra"))
write.csv(Y61_2_Y125_plot.data$pr_test_res, file = "./results/Y61_to_Y125_pr_test_res.csv")

plot.data <- Y61_2_Y122_plot.data$data.plot
plot.data <- pmin(plot.data, quantile(plot.data, .99))
plot.data <- pmax(plot.data, quantile(plot.data, .01))
anno_df <- data.frame(
  Tissue = colData(eos.cds)[Y61_2_Y122_plot.data$cell_order, "Tissue"],
  Pseudotime = pseudotime(eos.cds)[Y61_2_Y122_plot.data$cell_order]
)
ha.col <- HeatmapAnnotation(
  df = anno_df,
  col = list(
    Tissue = tissue_color_panel,
    Pseudotime = colorRamp2(seq(min(anno_df$Pseudotime), max(anno_df$Pseudotime), length = 50), colorRampPalette(brewer.pal(9, "BuPu")[2:8])(50))
  ),
  annotation_height = unit(c(0.5, 0.5), "cm")
)
p <- ComplexHeatmap::Heatmap(
  plot.data,
  col = colorRamp2(seq(min(plot.data), max(plot.data), length = 50), viridis::viridis(50)),
  show_row_names = T,
  show_column_names = F,
  cluster_rows = F,
  cluster_columns = F,
  top_annotation = ha.col
)
draw(p)

plot.data <- Y61_2_Y125_plot.data$data.plot
plot.data <- pmin(plot.data, quantile(plot.data, .99))
plot.data <- pmax(plot.data, quantile(plot.data, .01))
anno_df <- data.frame(
  Tissue = colData(eos.cds)[Y61_2_Y125_plot.data$cell_order, "Tissue"],
  Pseudotime = pseudotime(eos.cds)[Y61_2_Y125_plot.data$cell_order]
)
ha.col <- HeatmapAnnotation(
  df = anno_df,
  col = list(
    Tissue = tissue_color_panel,
    Pseudotime = colorRamp2(seq(min(anno_df$Pseudotime), max(anno_df$Pseudotime), length = 50), colorRampPalette(brewer.pal(9, "BuPu")[2:8])(50))
  ),
  annotation_height = unit(c(0.5, 0.5), "cm")
)
p <- ComplexHeatmap::Heatmap(
  plot.data,
  col = colorRamp2(seq(min(plot.data), max(plot.data), length = 50), viridis::viridis(50)),
  show_row_names = T,
  show_column_names = F,
  cluster_rows = F,
  cluster_columns = F,
  top_annotation = ha.col
)
draw(p)