library(tidyverse)
library(ggrepel)
library(cowplot)
library(ggpubr)
library(Seurat)
library(scCustomize)

eos_0.8 <- readRDS("./RDS/Fig1e_sc_eosinophil.rds")

#Figure 7a-----
my_colors <- get_palette(palette = "BuPu",k=9)
FeaturePlot_scCustom(eos_0.8, features = "Csf2ra", colors_use = my_colors)

#Figure 7b-----
tissue_color_panel <- c("BM" = "#70beb5", "Blood" = "#85b2da",
                        "eWAT" = "#eaa3a5", "scWAT" = "#f5d1d2",
                        "BAT" = "#f4c09f", "SI" = "#a2a2ca",
                        "Colon" = "#9f7ca9", "Lung" = "#72b231")

eos_all_exp <- as.data.frame(eos_0.8@assays$RNA$data)
all_Csf2ra_exp <- eos_all_exp["Csf2ra", ] %>%
  t() %>% as.data.frame() %>% 
  rownames_to_column(var="CellName")

eos_all_metadata <- eos_0.8@meta.data
Csf2ra_group <- eos_all_metadata[ ,c(7,8)]
rownames(Csf2ra_group) <- NULL
plotdata <- merge(all_Csf2ra_exp,Csf2ra_group,by="CellName")

comparisons <- combn(levels(plotdata$Tissue), 2, simplify = FALSE)

ggplot(data = plotdata, aes(x = Tissue, y = Csf2ra, fill= Tissue)) +
  geom_violin(aes(fill = Tissue),color="white",width = 1,scale="width") +  
  geom_boxplot(width = .05, notch = F, lwd=0.25) +
  scale_fill_manual(values = tissue_color_panel) +
  stat_compare_means(comparisons = comparisons,
                     method = "wilcox.test",
                     label = "p.signif",
                     vjust= 0.6) +
  labs(y = "Csf2ra expression level") + 
  scale_y_continuous(expand = expansion(mult = c(0.01, 0.1))) +
  theme_cowplot(font_size = 7)

