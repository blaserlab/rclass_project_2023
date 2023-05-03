# Attach the packages you will need for the analysis.
library(blaseRtools)
library(blaseRdata)
library(monocle3)
library(Seurat)
library(tidyverse)
library(cowplot)
library(ggplot2)
theme_set(theme_cowplot(font_size = 12))



## review laod_process.R

vig_qc_res

# Cell metadata
bb_cellmeta(vignette_cds)

# Gene metadata
bb_rowmeta(vignette_cds)

# Visualizations
bb_var_umap(vignette_cds, var = "sample")

bb_var_umap(vignette_cds, var = "sample",
            alt_dim_x = "prealignment_dim1",
            alt_dim_y = "prealignment_dim2")

# Clustering
##  top markers
vignette_top_markers
bb_var_umap(vignette_cds, var = "partition")
bb_var_umap(vignette_cds, var = "leiden")
bb_var_umap(vignette_cds, var = "louvain")

# gene modules
bb_rowmeta(vignette_cds)

# label transfer using Seurat
bb_var_umap(vignette_cds,
            var = "seurat_celltype_l1",
            alt_dim_x = "seurat_dim1",
            alt_dim_y = "seurat_dim2",
            overwrite_labels = TRUE,
            group_label_size = 4)

bb_var_umap(vignette_cds, var = "seurat_celltype_l1")

# cell assignments
bb_var_umap(vignette_cds, var = "leiden", plot_title = "Leiden Clusters")
leiden_seurat <- bb_cellmeta(vignette_cds) %>%
  group_by(leiden, seurat_celltype_l1) %>%
  summarise(n = n())
leiden_seurat
ggplot(leiden_seurat,
       mapping = aes(x = leiden,
                     y = n,
                     fill = seurat_celltype_l1)) +
  geom_bar(position="fill", stat="identity") +
  scale_fill_brewer(palette = "Set1")

# Recode the leiden clusters with our cell assignments
colData(vignette_cds)$leiden_assignment <- recode(colData(vignette_cds)$leiden,
                                                  "1" = "T/NK",
                                                  "2" = "DC/Mono",
                                                  "3" = "B")
bb_var_umap(vignette_cds, var = "leiden_assignment")


# UMAP plot types

bb_var_umap(vignette_cds,
            var = "leiden_assignment",
            value_to_highlight = "T/NK",
            cell_size = 2,
            foreground_alpha = 0.4)

bb_var_umap(vignette_cds,
            var = "leiden_assignment",
            value_to_highlight = "T/NK",
            palette = "green4",
            cell_size = 2,
            foreground_alpha = 0.4)

bb_var_umap(vignette_cds,
            var = "density",
            facet_by = "equipment",
            cell_size = 2,
            plot_title = "Local Cell Density")

bb_cellmeta(vignette_cds) |> glimpse()

# for demonstration only
colData(vignette_cds)$demo_sample <- paste0(colData(vignette_cds)$sample, "_", 1:3)

bb_cluster_representation2(obj = vignette_cds,
                          sample_var = "demo_sample",
                          cluster_var = "leiden_assignment",
                          comparison_var = "equipment",
                          comparison_levels = c("chromium", "X"),
                          return_val = "plot",
                          sig_val = "PValue")


bb_gene_umap(vignette_cds,
             gene_or_genes = "CD3D")
bb_gene_umap(vignette_cds, gene_or_genes = c("CD19", "CD3D", "CD14"))
agg_genes <-
  bb_rowmeta(vignette_cds) %>%
  select(id, module_labeled) %>%
  filter(module_labeled == "Module 1")


bb_gene_umap(vignette_cds,
             gene_or_genes = agg_genes)

# gene bubbles
bb_genebubbles(vignette_cds,
               genes = c("CD3E", "CD14", "CD19"),
               cell_grouping = "leiden_assignment",
               scale_expr = TRUE)
bb_genebubbles(vignette_cds,
               genes = c("CD3E", "CD14", "CD19"),
               cell_grouping = "leiden_assignment",
               scale_expr = FALSE)
# specify the order of the genes
bb_genebubbles(vignette_cds,
               genes = c("CD3E", "CD14", "CD19"),
               cell_grouping = "leiden_assignment",
               gene_ordering = c("as_supplied"),
               group_ordering = c("as_supplied"),
               scale_expr = FALSE)

