library(Seurat)
library(scCustomize)
library(harmony)
library(ggpubr)
c_ndcirco <- c("#fb942a","#d945a0","#415f96","#c8ebe8","#f2a7b8",
               "#516dd3","#f5f179","#fbc160","#5cbac8","#97c6df",
               "#f97d3c","#8457c0","#dbc5d5","#7437c7","#e93e83")

#Figure 1b-----
eos_all <- readRDS("./RDS/Fig1b_sc_all.rds")
DimPlot_scCustom(eos_all)+ scale_color_manual(values = c_ndcirco)

#Figure 1c----
eos_all_markers <- FindAllMarkers(eos_all)
genes_used <- c(eos_all_markers %>% group_by(cluster) %>% filter(!grepl("^Gm|^ENS|Rik|mt|^AI",gene), pct.1 > .3) %>% slice_max(order_by = avg_log2FC, n = 10) %>% pull(gene),"Ccl22","Cd36","Siglecf","Il5ra","Ccr3","Epx") %>% unique()
markers_mean_exp <- 
  aggregate(
    as.matrix(t(eos_all@assays$RNA@data[genes_used,])),
    list(Cluster = eos_all@meta.data[,"Anno_Raw"]),
    mean)
row.names(markers_mean_exp) <- markers_mean_exp$Cluster
markers_mean_exp$Cluster <- c()
markers_plot_matrix <- apply(markers_mean_exp, 2, zscore) %>% t()
max.avg <- apply(markers_plot_matrix, 1, which.max)
gene_order <- c()
for(i in 1:ncol(markers_plot_matrix)){
  if(sum(max.avg == i) > 1){
    temp <- data.frame(
      gene = names(sort(markers_plot_matrix[names(max.avg)[max.avg == i],i], decreasing = T)), 
      cluster = colnames(markers_plot_matrix)[i], stringsAsFactors = F)
  }else if(sum(max.avg == i) == 1){
    temp <- data.frame(gene = names(max.avg)[max.avg == i], 
                       cluster = colnames(markers_plot_matrix)[i])
  }else{
    temp <- c()
  }
  gene_order <- rbind(gene_order, temp)
}
gene_order <- gene_order[nrow(gene_order):1,]
genes_labeled <- c("Pecam1","Gata2","Cpa3","Ighd","Ms4a1","Cd19","Ccr7","Ccl22","Fscn1","Il2ra","Ccr8","Nrgn","Igkc","Igha","Jchain","Cstb","Retnla","Cd36","C1qc","Folr2","Cd163","Ace","Ms4a4c","Chil3","Cxcr2","S100a8","Il36g","Cd3g","Tcf7","Cd3e","Dpp4","Flt3","Kmo","Sema7a","Ccr3","Prg2","Siglecf","Il5ra","Alox15","Epx")
color_used <- circlize::colorRamp2(seq(min(markers_plot_matrix), max(markers_plot_matrix), length = 8), RColorBrewer::brewer.pal(9,'BuPu')[1:8])
p <- Heatmap(markers_plot_matrix[gene_order$gene,],
             cluster_rows = FALSE, 
             cluster_columns = FALSE,
             show_row_names = T,
             column_names_gp = gpar(fontsize = 10),
             row_names_gp = gpar(fontsize = 6),
             col = color_used,
             name = "Exp") +
  rowAnnotation(link = anno_mark(
    at = which(gene_order$gene %in% genes_labeled),
    labels = gene_order$gene[which(gene_order$gene %in% genes_labeled)],
    labels_gp = gpar(fontsize = 10), padding = unit(1, "mm")))
draw(p)

#Figure 1d------
my_colors <- get_palette(palette = "BuPu",k=9)
FeaturePlot_scCustom(eos_all,features = c("Siglecf","Il5ra",
                                      "Epx","Ccr3",
                                      "Retnlg","Alox15"),
                     colors_use = my_colors)

#Figure 1e------
eos_0.8 <- readRDS("./RDS/Fig1e_sc_eosinophil.rds")
DimPlot_scCustom(eos_0.8)+ scale_color_manual(values = c_ndcirco)

#Figure 1g------
library(ComplexHeatmap)
library(circlize)
library(tidyverse)
library(RColorBrewer)

Exp <- as.data.frame(eos_0.8@assays$RNA$data) %>% as.matrix()
Group <- as.data.frame(eos_0.8@meta.data)
Group <- Group[ ,c(1,6,7,8)]
rownames(Group)<-NULL
colnames(Group)[1]<-"Group"
BM <- Group %>% mutate(Tissue = ifelse(Tissue == "BM", "BM", "others"))
Blood <- Group %>% mutate(Tissue = ifelse(Tissue == "Blood", "Blood", "others"))
eWAT <- Group %>% mutate(Tissue = ifelse(Tissue == "eWAT", "eWAT", "others"))
BAT <- Group %>% mutate(Tissue = ifelse(Tissue == "BAT", "BAT", "others"))
scWAT <- Group %>% mutate(Tissue = ifelse(Tissue == "scWAT", "scWAT", "others"))
SI <- Group %>% mutate(Tissue = ifelse(Tissue == "SI", "SI", "others"))
Colon <- Group %>% mutate(Tissue = ifelse(Tissue == "Colon", "Colon", "others"))
Lung <- Group %>% mutate(Tissue = ifelse(Tissue == "Lung", "Lung", "others"))

## LIMMA
LIMMA <- function(expression_matrix, groupid, doAUC = FALSE) {
  ## Differential expressed genes in two groups by Limma.
  ##
  ## Args:
  #' @expression_matrix: Gene*cell matrix.
  #' @groupid: The groupid of each cell, there should be only two groups.
  ##
  ## Returns:
  ## A dataframe with the output of limma and expression percentage.
  library(limma)
  library(ROCR)
  expression_matrix <- as.matrix(expression_matrix)
  groupid <- as.character(groupid)
  groupid_raw <- unique(groupid)
  groupid[groupid == groupid_raw[1]] <- "GroupA"
  groupid[groupid == groupid_raw[2]] <- "GroupB"
  contrast <<- paste0(levels(factor(groupid)), collapse = "-")
  design <- model.matrix( ~ 0 + factor(groupid))
  colnames(design) <- levels(factor(groupid))
  rownames(design) <-
    colnames(expression_matrix)  # design data used in limma
  contrast.matrix <- makeContrasts(contrast, levels = design)
  fit <- lmFit(expression_matrix, design)
  fit2 <- contrasts.fit(fit, contrast.matrix)
  fit2 <- eBayes(fit2)
  tempOutput <-  topTable(fit2, coef = 1, n = Inf)
  nrDEG <-  na.omit(tempOutput)
  nrDEG$Symbol <- row.names(nrDEG)
  positive_group <-
    row.names(fit2$contrasts)[fit2$contrasts == 1]  # high expression when logFC > 0
  negative_group <-
    row.names(fit2$contrasts)[fit2$contrasts == -1]  # low expression when logFC < 0
  nrDEG$Grp <-
    c(negative_group, positive_group)[as.numeric(nrDEG$logFC > 0) + 1]
  nrDEG$Grp[nrDEG$Grp == "GroupA"] <- groupid_raw[1]
  nrDEG$Grp[nrDEG$Grp == "GroupB"] <- groupid_raw[2]
  cell.Grp1 <- which(groupid == levels(as.factor(groupid))[1])
  cell.Grp2 <- which(groupid == levels(as.factor(groupid))[2])
  Exp.Mean.Grp1 <-
    rowMeans(expression_matrix[nrDEG$Symbol, cell.Grp1])  # mean expression in the group
  Exp.Mean.Grp2 <-
    rowMeans(expression_matrix[nrDEG$Symbol, cell.Grp2])  # mean expression out the group
  Exp.Per.Grp1 <-
    apply(expression_matrix[nrDEG$Symbol, cell.Grp1], 1, Expression.Per)  # expression percentage in the group
  Exp.Per.Grp2 <-
    apply(expression_matrix[nrDEG$Symbol, cell.Grp2], 1, Expression.Per)  # expression percentage out the group
  de.genes.all <- cbind(nrDEG, Exp.Mean.Grp1, Exp.Per.Grp1, Exp.Mean.Grp2, Exp.Per.Grp2)
  colnames(de.genes.all)[9:12] <-
    c(paste0(c("Exp.Mean.", "Exp.Per."), groupid_raw[1]),
      paste0(c("Exp.Mean.", "Exp.Per."), groupid_raw[2]))
  if (!is.na(de.genes.all[1, 1])) {
    if(doAUC){
      for (k in 1:nrow(de.genes.all)) {
        category <- as.numeric(groupid_new == de.genes.all$Grp[k])
        pred <-
          prediction(expression_matrix[de.genes.all$Symbol[k], ], category)
        pauc <- performance(pred, measure = "auc")
        de.genes.all[k, "AUC"] <- pauc@y.values[[1]]
      }  # calculate the AUC for each gene
    }
    de.genes.all <- de.genes.all %>% arrange(desc(logFC))
    de.genes.all
    return(de.genes.all)
  } else{
    print("No significant genes!")
  }
  
}
Expression.Per <- function(x, cutoff = 0.1) {
  # percent of gene-expressed cell, cutoff = 0.1
  return(sum(x > cutoff) / length(x))
}

groupid <- BM$Tissue[match(colnames(Exp), BM$CellName)]
DEGs_BM <- LIMMA(Exp,groupid)

groupid <- Blood$Tissue[match(colnames(Exp), Blood$CellName)]
DEGs_Blood<-LIMMA(Exp,groupid)

groupid <- eWAT$Tissue[match(colnames(Exp), eWAT$CellName)]
DEGs_eWAT<-LIMMA(Exp,groupid)

groupid <- BAT$Tissue[match(colnames(Exp), BAT$CellName)]
DEGs_BAT<-LIMMA(Exp,groupid)

groupid <- scWAT$Tissue[match(colnames(Exp), scWAT$CellName)]
DEGs_scWAT<-LIMMA(Exp,groupid)

groupid <- SI$Tissue[match(colnames(Exp), SI$CellName)]
DEGs_SI<-LIMMA(Exp,groupid)

groupid <- Colon$Tissue[match(colnames(Exp), Colon$CellName)]
DEGs_Colon<-LIMMA(Exp,groupid)

groupid <- Lung$Tissue[match(colnames(Exp), Lung$CellName)]
DEGs_Lung<-LIMMA(Exp,groupid)

filter_BAT <- DEGs_BAT[!grepl("^mt-|^Rps|^Rpl|^ENS", rownames(DEGs_BAT)), ] %>% 
  filter(adj.P.Val < 0.05) %>% 
  slice_min(order_by = logFC, n = 10)

filter_Blood <- DEGs_Blood[!grepl("^mt-|^Rps|^Rpl|^ENS", rownames(DEGs_Blood)), ] %>% 
  filter(adj.P.Val < 0.05) %>% 
  slice_min(order_by = logFC, n = 10)

filter_BM <- DEGs_BM[!grepl("^mt-|^Rps|^Rpl|^ENS", rownames(DEGs_BM)), ] %>% 
  filter(adj.P.Val < 0.05) %>% 
  slice_min(order_by = logFC, n = 10)

filter_Colon <- DEGs_Colon[!grepl("^mt-|^Rps|^Rpl|^ENS", rownames(DEGs_Colon)), ] %>% 
  filter(adj.P.Val < 0.05) %>% 
  slice_min(order_by = logFC, n = 10)

filter_eWAT <- DEGs_eWAT[!grepl("^mt-|^Rps|^Rpl|^ENS", rownames(DEGs_eWAT)), ] %>%
  filter(adj.P.Val < 0.05) %>% 
  slice_max(order_by = logFC, n = 10)

filter_Lung <- DEGs_Lung[!grepl("^mt-|^Rps|^Rpl|^ENS", rownames(DEGs_Lung)), ] %>% 
  filter(adj.P.Val < 0.05) %>% 
  slice_min(order_by = logFC, n = 10)

filter_scWAT <- DEGs_scWAT[!grepl("^mt-|^Rps|^Rpl|^ENS", rownames(DEGs_scWAT)), ] %>% 
  filter(adj.P.Val < 0.05) %>% 
  slice_min(order_by = logFC, n = 10)

filter_SI <- DEGs_SI[!grepl("^mt-|^Rps|^Rpl|^ENS", rownames(DEGs_SI)), ] %>% 
  filter(adj.P.Val < 0.05) %>% 
  slice_min(order_by = logFC, n = 10)

top10<-unlist(lapply(list(filter_BAT, filter_Blood, filter_BM,
                          filter_Colon, filter_eWAT, filter_Lung, 
                          filter_scWAT, filter_SI), function(x) x$Symbol))
top10<-unique(top10)

####heatmap
add_gene <- read.csv("./results/tissue based add DEGs.csv")
gene_used <- union(top10,add_gene$gene)

markers_mean_exp <- 
  aggregate(
    as.matrix(t(eos_0.8@assays$RNA@data[gene_used, ])),
    list(Tissue = eos_0.8@meta.data[ ,"Tissue"]),
    mean)
row.names(markers_mean_exp) <- markers_mean_exp$Tissue
markers_mean_exp$Tissue <- c()

## zscore
zscore <- function(x){
  return( (x - mean(x))/sd(x) )
}

markers_plot_matrix <- apply(markers_mean_exp, 2, zscore) %>% t()
max.avg <- apply(markers_plot_matrix, 1, which.max)

gene_order <- c()
for(i in 1:ncol(markers_plot_matrix)){
  if(sum(max.avg == i) > 1){
    temp <- data.frame(
      gene = names(sort(markers_plot_matrix[names(max.avg)[max.avg == i],i], decreasing = T)), 
      Tissue = colnames(markers_plot_matrix)[i], stringsAsFactors = F)
  }else if(sum(max.avg == i) == 1){
    temp <- data.frame(gene = names(max.avg)[max.avg == i], 
                       Tissue = colnames(markers_plot_matrix)[i])
  }else{
    temp <- c()
  }
  gene_order <- rbind(gene_order, temp)
}

color_used <- circlize::colorRamp2(seq(min(markers_plot_matrix), 
                                       max(markers_plot_matrix), 
                                       length = 8), 
                                   RColorBrewer::brewer.pal(9,'BuPu')[1:8])

plotdata <- markers_plot_matrix[gene_order$gene, ]

p <- Heatmap(plotdata, cluster_rows = FALSE, cluster_columns = FALSE,
             show_row_names = T, column_names_gp = gpar(fontsize = 10),
             row_names_gp = gpar(fontsize = 6), col = color_used, name = "Exp")

#Figure 1h------
csv_list <- list(BM = DEGs_BM, Blood = DEGs_Blood, eWAT = DEGs_eWAT,
                 BAT = DEGs_BAT, scWAT = DEGs_scWAT, SI = DEGs_SI,
                 Colon = DEGs_Colon, Lung = DEGs_Lung)

library(clusterProfiler)
library(org.Mm.eg.db)
all_genes <- unique(unlist(lapply(csv_list, rownames)))
universe_entrezid <- bitr(all_genes, fromType = "SYMBOL",
                          toType = c("ENTREZID"), 
                          OrgDb = org.Mm.eg.db)
universe_entrezid <- unique(universe_entrezid$ENTREZID)

gene_list <- list()
for(tissue in names(csv_list)){
  top_gene <-
    csv_list[[tissue]] %>%
    filter(!grepl("^mt-|^Rps|^Rpl|^ENS|Rik", Symbol), 
           adj.P.Val < 0.05, Grp != "others") %>%
    slice_max(order_by = abs(logFC), n = 200) %>%
    pull(Symbol)
  gene.df <- bitr(top_gene,
                  fromType = "SYMBOL",
                  toType = c("ENTREZID"),
                  OrgDb = org.Mm.eg.db)
  gene_list[[tissue]] <- unique(gene.df$ENTREZID)
}

cluster_go <- compareCluster(gene_list, fun = "enrichGO",
                             OrgDb = org.Mm.eg.db, ont = "BP",
                             universe = universe_entrezid,
                             pvalueCutoff = 0.05,
                             qvalueCutoff = 0.05,
                             readable = TRUE, pool = TRUE)

BP_result<-as.data.frame(cluster_go@compareClusterResult)

BP_result<-read.csv("./results/BP_compareCluster_results.csv",row.names = 1)

pathway<-read.csv("./results/chosen BP pathway new.csv")
pathway_used <- pathway %>%
  pull(Description) %>%
  unique()

heatmap_df <- BP_result %>% 
  filter(Description %in% pathway_used) %>%
  dplyr::select(Cluster, Description, p.adjust) %>%
  mutate(logp = -log10(p.adjust)) %>%
  dplyr::select(-p.adjust) %>%
  pivot_wider(names_from = Cluster, values_from = logp, values_fill = 0) %>%
  tibble::column_to_rownames(var = "Description")

Tissue_order <- c("BM","Blood","eWAT","BAT","scWAT","SI","Colon","Lung")
heatmap_df <- heatmap_df[ ,Tissue_order]
markers_plot_matrix <- apply(heatmap_df, 1, zscore) %>% t()
max.avg <- apply(markers_plot_matrix, 1, which.max)

pathway_order <- c()
for(i in 1:ncol(markers_plot_matrix)){
  if(sum(max.avg == i) > 1){
    temp <- data.frame(
      pathway = names(sort(markers_plot_matrix[names(max.avg)[max.avg == i],i], decreasing = T)), 
      tissue = colnames(markers_plot_matrix)[i], stringsAsFactors = F)
  }else if(sum(max.avg == i) == 1){
    temp <- data.frame(pathway = names(max.avg)[max.avg == i], 
                       tissue = colnames(markers_plot_matrix)[i])
  }else{
    temp <- c()
  }
  pathway_order <- rbind(pathway_order, temp)
}
pathway_order <- pathway_order[nrow(pathway_order):1,]
pathways_labeled <- pathway_used

color_used <- circlize::colorRamp2(seq(min(markers_plot_matrix), 
                                       max(markers_plot_matrix), 
                                       length = 8), 
                                   RColorBrewer::brewer.pal(9,'BuPu')[1:8])
plotdata_BP<-as.data.frame(markers_plot_matrix[pathway_order$pathway,])

p <- Heatmap(as.matrix(plotdata_BP),
             cluster_rows = FALSE, 
             cluster_columns = FALSE,
             show_row_names = T,
             column_names_gp = gpar(fontsize = 10),
             row_names_gp = gpar(fontsize = 9),
             row_names_max_width = max_text_width(
               rownames(plotdata_BP)),
             col = color_used,
             show_heatmap_legend = FALSE)

#Figure 1i------
scenic<-read.csv("./results/eos.aucell.csv",row.names = 1,check.names = F)
metadata<-eos_0.8@meta.data
metadata <- metadata[ ,c(7,8)]
rownames(metadata)<- NULL

TF_mean <-aggregate(as.matrix(scenic),
                    list(Tissue = metadata[ ,"Tissue"]), mean)

rownames(TF_mean) <- TF_mean$Tissue
TF_mean$Tissue <- c()
TF_plot_matrix <- apply(TF_mean, 2, zscore) %>% t()
max.avg <- apply(TF_plot_matrix, 1, which.max)

TF_order <- c()
for(i in 1:ncol(TF_plot_matrix)){
  if(sum(max.avg == i) > 1){
    temp <- data.frame(
      TF = names(sort(TF_plot_matrix[names(max.avg)[max.avg == i],i], decreasing = T)), 
      tissue = colnames(TF_plot_matrix)[i], stringsAsFactors = F)
  }else if(sum(max.avg == i) == 1){
    temp <- data.frame(TF = names(max.avg)[max.avg == i], 
                       tissue = colnames(TF_plot_matrix)[i])
  }else{
    temp <- c()
  }
  TF_order <- rbind(TF_order, temp)
}
TF_order <- TF_order[nrow(TF_order):1,]
color_used <- circlize::colorRamp2(seq(min(TF_plot_matrix), 
                                       max(TF_plot_matrix), 
                                       length = 8),
                                   RColorBrewer::brewer.pal(9,'BuPu')[1:8])
p <- Heatmap(TF_plot_matrix[TF_order$TF, ],
             cluster_rows = F, 
             cluster_columns = F,
             show_row_names = T,
             column_names_gp = gpar(fontsize = 10),
             row_names_gp = gpar(fontsize = 6),
             col = color_used,
             name = "Exp")

draw(p)

#Sup Figure 1d-----
DEGs <- FindAllMarkers(eos_0.8, only.pos = TRUE)
add_gene <- read.csv("cluster based add DEGs.csv")
top10 <- DEGs %>% filter(!grepl("^mt-|^Rps|^Rpl|^ENS",gene), pct.1 > .3) %>% 
  filter(p_val_adj < 0.05) %>% group_by(cluster) %>% 
  slice_max(order_by = avg_log2FC, n = 10)

gene_used <- union(unique(top10$gene),add_gene$gene)

markers_mean_exp <- 
  aggregate(
    as.matrix(t(eos_0.8@assays$RNA@data[gene_used, ])),
    list(Cluster = eos_0.8@meta.data[ ,"New_clusters"]),
    mean)
row.names(markers_mean_exp) <- markers_mean_exp$Cluster
markers_mean_exp$Cluster <- c()

## zscore
zscore <- function(x){
  return( (x - mean(x))/sd(x) )
}

markers_plot_matrix <- apply(markers_mean_exp, 2, zscore) %>% t()
max.avg <- apply(markers_plot_matrix, 1, which.max)

gene_order <- c()
for(i in 1:ncol(markers_plot_matrix)){
  if(sum(max.avg == i) > 1){
    temp <- data.frame(
      gene = names(sort(markers_plot_matrix[names(max.avg)[max.avg == i],i], decreasing = T)), 
      cluster = colnames(markers_plot_matrix)[i], stringsAsFactors = F)
  }else if(sum(max.avg == i) == 1){
    temp <- data.frame(gene = names(max.avg)[max.avg == i], 
                       cluster = colnames(markers_plot_matrix)[i])
  }else{
    temp <- c()
  }
  gene_order <- rbind(gene_order, temp)
}
gene_order <- gene_order[nrow(gene_order):1,]
color_used <- circlize::colorRamp2(seq(min(markers_plot_matrix), 
                                       max(markers_plot_matrix), 
                                       length = 8), 
                                   RColorBrewer::brewer.pal(9,'BuPu')[1:8])
plotdata <- markers_plot_matrix[gene_order$gene, ]

Heatmap(plotdata, cluster_rows = FALSE, cluster_columns = FALSE,
        show_row_names = T, column_names_gp = gpar(fontsize = 10),
        row_names_gp = gpar(fontsize = 6), col = color_used, name = "Exp")

#Sup Figure 1e-----
BP_result<-read.csv("BP_compareCluster_results.csv",row.names = 1)
pathway<-read.csv("chosen BP pathway new.csv")
pathway_used <- pathway %>%
  pull(Description) %>%
  unique()
BubblePlot_df <- BP_result %>% 
  filter(Description %in% pathway_used) %>%
  dplyr::select(Cluster, Description, GeneRatio, p.adjust)

BubblePlot_df$GeneRatio <- apply(str_split(BubblePlot_df$GeneRatio, "/", simplify = T), 2,as.numeric)
BubblePlot_df$GeneRatio <- BubblePlot_df$GeneRatio[,1]/BubblePlot_df$GeneRatio[,2]

Tissue_order <- c("BM","Blood","eWAT","BAT","scWAT","SI","Colon","Lung")
BubblePlot_df$Cluster <- factor(BubblePlot_df$Cluster,levels = Tissue_order)

pathway_order<-read.csv("Sup Figure 1e pathway_order.csv",row.names = 1)
pathway_order<-pathway_order$x
BubblePlot_df$Description <- factor(BubblePlot_df$Description,levels = pathway_order)
color_used <- c("#E0ECF4FF","#BFD3E6FF","#9EBCDAFF","#8C96C6FF","#8C6BB1FF","#88419DFF", "#810F7CFF")
ggplot(BubblePlot_df, aes(factor(Cluster),factor(Description))) +
  geom_point(aes(size = GeneRatio, color = -log10(p.adjust)), show.legend = TRUE) +
  scale_color_gradientn(colours = color_used,
                        guide = guide_colorbar(reverse = F,order = 1))+
  scale_size_continuous(range = c(3, 8)) +  
  theme_minimal()+
  labs(x = "", y = "")

#Sup Figure 1f-----
eos_0.8_filtered <- eos_0.8 %>% 
  subset(Tissue %ni% c("BAT","scWAT","eWAT"))
library(GEOquery)
Bulk<-read.csv("GSE211112_RNAseq_SS2.EOS_tissues.filt_tpm.csv",row.names = 1)
gseSet <- getGEO("GSE211112",destdir = '.',getGPL = F)
gseSet <- gseSet[[1]]
bulk_meta <- pData(gseSet)
bulk_meta <- bulk_meta[c(1:24),c(1,8)]
colnames(bulk_meta) <- c("Sample", "Tissue")
bulk_meta <- bulk_meta[!bulk_meta$Tissue %in% c("Spleen","Skin"), ]
bulk_meta <- bulk_meta %>%
  mutate(across(2, ~case_when(
    . == "Bone Marrow" ~ "BM",
    . == "Small Intestine" ~ "SI",
    TRUE ~ .
  )))
rownames(bulk_meta) <- NULL

Bulk<-Bulk[ ,colnames(Bulk) %in% bulk_meta$Sample]
Bulk <- log2(Bulk + 1)
identical(colnames(Bulk),bulk_meta$Sample)

campare_result <- Compare2Clusters(eos_0.8_filtered@assays$RNA@data, 
                                   eos_0.8_filtered@meta.data$Tissue,
                                   Bulk, bulk_meta$Tissue, ngenes = 1000,
                                   name1 = "sc", name2 = "Bulk")

library(ggdendro)
dend_data <- ggdendro::dendro_data(campare_result$hc_fit$hclust, type = "rectangle")
dend_data$labels$dataset <- stringr::str_split_fixed(dend_data$labels$label, "_", 2)[,1]
ggplot() +
  geom_segment(data = dend_data$segments,
               aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_text(data = dend_data$labels,
            aes(x = x, y = y - 0.05 * max(dend_data$segments$y), label = label, color = dataset),
            angle = 90, hjust = 1, size = 3) +
  scale_y_continuous(expand = c(1.2,0)) +
  scale_x_continuous(expand = c(0,1.5)) +
  scale_color_manual(values = c("bulk" = "#33583a", "sc" = "#72477d")) +
  theme_nothing() +
    coord_cartesian(ylim = c(min(dend_data$labels$y) - 0.1 * max(dend_data$segments$y),
                           max(dend_data$segments$y))) +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank())
