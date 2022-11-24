#Load packages
library(Seurat)
library(SeuratData)
library(dplyr)
library(SPOTlight)
library(igraph)
library(RColorBrewer)

#Load files
DBA2J_relabellednew <- readRDS("C:/Users/lgmheezen/OneDrive - LUMC/Ahmed_Laura_Pietro/Deconvolution/Spotlight/211108_DBA2J_relabellednew.rds")
combined <- readRDS("C:/Users/lgmheezen/OneDrive - LUMC/Ahmed_Laura_Pietro/Deconvolution/Spotlight/combined.rds")
cluster_markers_all <- readRDS("C:/Users/lgmheezen/OneDrive - LUMC/Ahmed_Laura_Pietro/Deconvolution/Spotlight/cluster_markers_all.rds")



spotlight_ls_DBA <- spotlight_deconvolution(
  se_sc = combined,
  counts_spatial = DBA2J_relabellednew@assays$Spatial@counts,
  clust_vr = "cluster", # Variable in sc_seu containing the cell-type annotation
  cluster_markers = cluster_markers_all, # Dataframe with the marker genes
  cl_n = 70, # number of cells per cell type to use
  hvg = 3000, # Number of HVG to use
  ntop = NULL, # How many of the marker genes to use (by default all)
  transf = "uv", # Perform unit-variance scaling per cell and spot prior to factorzation and NLS
  method = "nsNMF", # Factorization method
  min_cont = 0 # Remove those cells contributing to a spot below a certain threshold 
)

saveRDS(object = spotlight_ls_DBA, file = "C:/Users/lgmheezen/OneDrive - LUMC/Ahmed_Laura_Pietro/Deconvolution/Spotlight/Laura/spotlight_ls_DBA2J.rds")
#spotlight_ls <- readRDS("Spotlight/spotlight_ls_DBA2J.rds")

nmf_mod <- spotlight_ls_DBA[[1]]
decon_mtrx <- spotlight_ls_DBA[[2]]

h <- NMF::coef(nmf_mod[[1]])
rownames(h) <- paste("Topic", 1:nrow(h), sep = "_")
topic_profile_plts <- SPOTlight::dot_plot_profiles_fun(
  h = h,
  train_cell_clust = nmf_mod[[2]])

topic_profile_plts[[2]] + ggplot2::theme(
  axis.text.x = ggplot2::element_text(angle = 90), 
  axis.text = ggplot2::element_text(size = 12))

topic_profile_plts[[1]] + theme(axis.text.x = element_text(angle = 90), 
                                axis.text = element_text(size = 12))

basis_spotlight <- data.frame(NMF::basis(nmf_mod[[1]]))

colnames(basis_spotlight) <- unique(stringr::str_wrap(nmf_mod[[2]], width = 30))

basis_spotlight %>%
  dplyr::arrange(desc(FAP)) %>%
  round(., 5) %>% 
  DT::datatable(., filter = "top")

decon_mtrx_sub <- decon_mtrx[, colnames(decon_mtrx) != "res_ss"]
decon_mtrx_sub[decon_mtrx_sub < 0.08] <- 0
decon_mtrx <- cbind(decon_mtrx_sub, "res_ss" = decon_mtrx[, "res_ss"])
rownames(decon_mtrx) <- colnames(DBA2J_relabellednew)

decon_df <- decon_mtrx %>%
  data.frame() %>%
  tibble::rownames_to_column("barcodes")

DBA2J_relabellednew@meta.data <- DBA2J_relabellednew@meta.data %>%
  tibble::rownames_to_column("barcodes") %>%
  dplyr::left_join(decon_df, by = "barcodes") %>%
  tibble::column_to_rownames("barcodes")

Seurat::SpatialFeaturePlot(
  object = DBA2J_relabellednew,
  features = c("IIa", "IIx", "MuSC", "IIb", "FAP",
               "EC", "SMC", "Myob", "MTJ", "IIx_b", "NMJ", "TC", "MPH", "RegMyon"),
  alpha = c(0.1, 1))

cell_types_all <- colnames(decon_mtrx)[which(colnames(decon_mtrx) != "res_ss")]

SPOTlight::spatial_scatterpie(se_obj = DBA2J_relabellednew,
                              cell_types_all = cell_types_all,
                              img_path = "C:/Users/lgmheezen/OneDrive - LUMC/Ahmed_Laura_Pietro/Deconvolution/Spotlight/tissue_lowres_image_DBA2J.png",
                              pie_scale = 0.4)

SPOTlight::spatial_scatterpie(se_obj = DBA2J_relabellednew,
                              cell_types_all = cell_types_all,
                              img_path = "C:/Users/lgmheezen/OneDrive - LUMC/Ahmed_Laura_Pietro/Deconvolution/Spotlight/tissue_lowres_image_DBA2J.png",
                              cell_types_interest = "FAP",
                              pie_scale = 0.8)

SPOTlight::spatial_scatterpie(se_obj = DBA2J_relabellednew,
                              cell_types_all = cell_types_all,
                              img_path = "C:/Users/lgmheezen/OneDrive - LUMC/Ahmed_Laura_Pietro/Deconvolution/Spotlight/tissue_lowres_image_DBA2J.png",
                              cell_types_interest = "MPH", 
                              pie_scale = 0.8)

graph_ntw <- SPOTlight::get_spatial_interaction_graph(decon_mtrx = decon_mtrx[, cell_types_all])
deg <- degree(graph_ntw, mode="all")
# Get color palette for difusion
edge_importance <- E(graph_ntw)$importance

# Select a continuous palette
qual_col_pals <- brewer.pal.info[brewer.pal.info$category == 'seq',]

# Create a color palette
getPalette <- colorRampPalette(brewer.pal(9, "YlOrRd"))

# Get how many values we need
grad_edge <- seq(0, max(edge_importance), 0.1)

# Generate extended gradient palette dataframe
graph_col_df <- data.frame(value = as.character(grad_edge),
                           color = getPalette(length(grad_edge)),
                           stringsAsFactors = FALSE)
# Assign color to each edge
color_edge <- data.frame(value = as.character(round(edge_importance, 1)), stringsAsFactors = FALSE) %>%
  dplyr::left_join(graph_col_df, by = "value") %>%
  dplyr::pull(color)

# Open a pdf file
plot(graph_ntw,
     # Size of the edge
     edge.width = edge_importance,
     edge.color = color_edge,
     # Size of the buble
     vertex.size = deg,
     vertex.color = "#cde394",
     vertex.frame.color = "white",
     vertex.label.color = "black",
     vertex.label.family = "Ubuntu", # Font family of the label (e.g."Times", "Helvetica")
     layout = layout.circle)

# Remove cell types not predicted to be on the tissue
decon_mtrx_sub <- decon_mtrx[, cell_types_all]
decon_mtrx_sub <- decon_mtrx_sub[, colSums(decon_mtrx_sub) > 0]

# Compute correlation
decon_cor <- cor(decon_mtrx_sub)

# Compute correlation P-value
p.mat <- corrplot::cor.mtest(mat = decon_mtrx_sub, conf.level = 0.95)

library(ggcorrplot)
# Visualize
ggcorrplot::ggcorrplot(
  corr = decon_cor,
  p.mat = p.mat[[1]],
  hc.order = TRUE,
  type = "full",
  insig = "blank",
  lab = TRUE,
  outline.col = "lightgrey",
  method = "square",
  # colors = c("#4477AA", "white", "#BB4444"))
  colors = c("#6D9EC1", "white", "#E46726"),
  title = "Predicted cell-cell proportion correlation",
  legend.title = "Correlation\n(Pearson)") +
  ggplot2::theme(
    plot.title = ggplot2::element_text(size = 22, hjust = 0.5, face = "bold"),
    legend.text = ggplot2::element_text(size = 12),
    legend.title = ggplot2::element_text(size = 15),
    axis.text.x = ggplot2::element_text(angle = 90),
    axis.text = ggplot2::element_text(size = 18, vjust = 0.5))



###Added by Laura

#Add in cluster label to compare cell type proportions for each cluster 
decon_df$clusterlabel <- DBA2J_relabellednew@meta.data[["relabelled"]]
write.csv(decon_df, "C:/Users/lgmheezen/OneDrive - LUMC/Ahmed_Laura_Pietro/Deconvolution/Spotlight/Laura/ProportionsClusters_DBA2J.csv") #save as csv


