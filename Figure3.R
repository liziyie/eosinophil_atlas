library(CytoTRACE2)
library(Seurat)
library(scCustomize)
library(RColorBrewer)
library(circlize)
library(monocle3)
library(harmony)
library(cowplot)
library(SeuratWrappers)

eos_0.8 <- readRDS("./RDS/Fig1e_sc_eosinophil.rds")

#Figure 3a----
eos.cds <- as.cell_data_set(eos_0.8)
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
  show_trajectory_graph = FALSE,
  label_cell_groups = FALSE, 
  label_leaves = FALSE, 
  label_branch_points = FALSE) +
  scale_color_gradientn(colours = myPalette(100)) +
  labs(x = "UMAP1", y = "UMAP2") +
  theme(legend.key.width = unit(0.5, "cm"),
        legend.key.height = unit(1, "cm"))

eos_0.8@meta.data$pseudotime <- pseudotime(eos.cds)[row.names(eos_0.8@meta.data)]
ggplot(eos_0.8@meta.data %>% subset(!is.infinite(pseudotime)), aes(x = reorder(Tissue,pseudotime), y = pseudotime)) +
  geom_boxplot(aes(fill = Tissue), outlier.size = .5) +
  scale_fill_manual(values = tissue_color_panel) +
  labs(x = "") +
  theme_cowplot(font_size = 7)

#Figure S4----
plot_cells(eos.cds, label_principal_points = TRUE)

Y46_2_Y64_plot.data <- monocle3_heatmap_1b(eos.cds, start.nodes = "Y_46", end.nodes = "Y_64", add.genes = c("Fcgr4","Clec12a","Cd22","Itgal","Csf2ra","Bcl2l11","Cdkn1a","Morrbid","Cflar","Birc3","Birc2"))
write.csv(Y46_2_Y64_plot.data$pr_test_res, file = "./results/Y46_to_Y64_pr_test_res.csv")
Y46_2_Y186_plot.data <- monocle3_heatmap_1b(eos.cds, start.nodes = "Y_46", end.nodes = "Y_186", add.genes = c("Clec12a","Cd22","Itgal","Csf2ra","Cdkn1a","Birc3"))
write.csv(Y46_2_Y186_plot.data$pr_test_res, file = "./results/Y46_to_Y186_pr_test_res.csv")
Y46_2_Y112_plot.data <- monocle3_heatmap_1b(eos.cds, start.nodes = "Y_46", end.nodes = "Y_112", add.genes = c("Clec12a","Itgal","Csf2ra","Bcl2l11","Cdkn1a","Morrbid","Birc3"))
write.csv(Y46_2_Y112_plot.data$pr_test_res, file = "./results/Y46_to_Y112_pr_test_res.csv")

anno_df <- data.frame(
  Tissue = colData(eos.cds)[Y46_2_Y64_plot.data$cell_order, "Tissue"],
  Pseudotime = pseudotime(eos.cds)[Y46_2_Y64_plot.data$cell_order]
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
  Y46_2_Y64_plot.data$data.plot,
  col = colorRamp2(seq(min(Y46_2_Y64_plot.data$data.plot), max(Y46_2_Y64_plot.data$data.plot), length = 50), viridis::viridis(50)),
  show_row_names = T,
  show_column_names = F,
  cluster_rows = F,
  cluster_columns = F,
  top_annotation = ha.col
)
draw(p)

anno_df <- data.frame(
  Tissue = colData(eos.cds)[Y46_2_Y186_plot.data$cell_order, "Tissue"],
  Pseudotime = pseudotime(eos.cds)[Y46_2_Y186_plot.data$cell_order]
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
  Y46_2_Y186_plot.data$data.plot,
  col = colorRamp2(seq(min(Y46_2_Y186_plot.data$data.plot), max(Y46_2_Y186_plot.data$data.plot), length = 50), viridis::viridis(50)),
  show_row_names = T,
  show_column_names = F,
  cluster_rows = F,
  cluster_columns = F,
  top_annotation = ha.col
)
draw(p)

anno_df <- data.frame(
  Tissue = colData(eos.cds)[Y46_2_Y112_plot.data$cell_order, "Tissue"],
  Pseudotime = pseudotime(eos.cds)[Y46_2_Y112_plot.data$cell_order]
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
  Y46_2_Y112_plot.data$data.plot,
  col = colorRamp2(seq(min(Y46_2_Y112_plot.data$data.plot), max(Y46_2_Y112_plot.data$data.plot), length = 50), viridis::viridis(50)),
  show_row_names = T,
  show_column_names = F,
  cluster_rows = F,
  cluster_columns = F,
  top_annotation = ha.col
)
draw(p)