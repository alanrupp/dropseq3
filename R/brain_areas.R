library(tidyverse)
library(Seurat)


# get scores across all clusters and areas
full_scores <- function(object, clusters = NULL, areas = NULL) {
  # grab the clusters
  if (is.null(clusters)) {
    clusters <- sort(unique(object@active.ident))
  } else {
    clusters <- clusters[clusters %in% object@active.ident]
  }
  
  # grab the area names based on the data
  if (is.null(areas)) {
    grab_areas <- function(areas = NULL, contrast = "HY") {
      # bring in all available areas
      all_areas <- list.files(path = paste0("~/Programs/dropseq3/data/Allen_Institute/", 
                                            contrast), 
                              pattern = "*.csv",
                              full.names = TRUE)
      area_names <- str_extract(all_areas, "[A-Za-z]+(?=\\.csv)")
      
      if (is.null(areas)) {
        areas <- area_names
      }
      
      # catch if entered area is not in names of available areas
      map(areas,
          ~ ifelse(.x %in% area_names, "", 
                   print(paste("Warning:", .x, "not available"))))
      areas <- areas[areas %in% area_names]
      return(areas)
    }
    areas <- grab_areas()
  }
  
  # 1) get average scaled difference for all genes within cluster
  enrichment_score <- function(object, cluster, area = NULL) {
    
    # grab enrichment scores
    grab_area_enrichment <- function(area, contrast = "HY") {
      df <- read_csv(paste0("~/Programs/dropseq3/data/Allen_Institute/", 
                      contrast, "/", area, ".csv"),
               progress = FALSE, col_types = cols())
      return(df)
      
    }
    
    # grab area enrichment scores
    enrichment <- grab_area_enrichment(area)
    
    # genes that are enriched
    genes <- unique(enrichment$"gene-symbol")
    genes <- genes[genes %in% rownames(object@assays$RNA@scale.data)]
    
    # grab cells
    cells <- names(object@active.ident)[object@active.ident == cluster]
    
    # calculate mean expression
    return(rowMeans(object@assays$RNA@scale.data[genes, cells]))
  }
  
  # score for multiple areas
  score_multiple <- function(object, cluster, areas) {
    score <- map(areas, ~ enrichment_score(object, cluster, .x))
    score <- map(score, ~ as.data.frame(.x))
    score <- map(score, ~ rename(.x, "score" = ".x"))
    return(map(seq(length(score)), ~mutate(score[[.x]], "area" = areas[.x])) %>%
      bind_rows()
    )
  }
  
  result <- map(clusters, ~ score_multiple(object, .x, areas))
  result <- map(1:length(result), ~ mutate(result[[.x]], "cluster" = clusters[.x]))
  return(result)
}

# summarize the scorse by mean z-score
summarize_scores <- function(scores) {
  scores %>%
    bind_rows() %>%
    group_by(area, cluster) %>%
    summarize(score = mean(score)) %>%
    arrange(desc(score)) %>%
    ungroup()
}

# get list of brain areas above 0 z-score
predict_brain_area <- function(scores_summary) {
  scores_summary %>%
    filter(score > 0) %>%
    group_by(cluster) %>%
    summarize(areas = list(area)) %>%
    mutate(areas = unlist(areas))
}

# plot heatmap
heatmap_plot <- function(object, genes = NULL, cells = NULL, scores,
                         scale_clip_low = -3, scale_clip_high = 3) {
  if (is.null(genes)) {
    genes <- object@var.genes
  }
  if (is.null(cells)) {
    cells <- sample(colnames(object@data), 1000)
  }
  
  # calculate mean z-score by gene and cluster --------------------------------
  cells_use <- function(object, cluster) {
    cell_names <- names(object@ident)[object@ident == cluster]
    cell_names <- cell_names[cell_names %in% cells]
  }
  
  # find mean scaled expression for genes by clusters
  mean_scale_data <-
    map(unique(object@ident),
        ~ rowMeans(object@scale.data[genes, cells_use(object, .x)])
        ) %>%
    as.data.frame() %>%
    set_names(unique(object@ident))
  
  # clip data with big outliers -----------------------------------------------
  mean_scale_data <- apply(mean_scale_data, c(1,2),
                           function(x) ifelse(x > scale_clip_high, 
                                              scale_clip_high, x))
  mean_scale_data <- apply(mean_scale_data, c(1,2),
                           function(x) ifelse(x < scale_clip_low, 
                                              scale_clip_low, x))
  
  # cluster using euclidean distance -----------------------------------------
  cluster_distances <- dist(t(mean_scale_data))
  cluster_clusters <- hclust(cluster_distances)
  cluster_order <- cluster_clusters$labels[cluster_clusters$order]
  
  gene_distances <- dist(mean_scale_data)
  gene_clusters <- hclust(gene_distances)
  gene_order <- gene_clusters$labels[gene_clusters$order]

  # make dendrogram
  dendro <- 
    ggdendro::ggdendrogram(cluster_clusters, rotate = TRUE) +
    theme_void()
  
  # make heatmap
  tiles <- 
    mean_scale_data %>%
    as.data.frame() %>%
    rownames_to_column("gene") %>%
    gather(-gene, key = "cluster", value = "expr") %>%
    mutate(cluster = factor(cluster, levels = cluster_order)) %>%
    mutate(gene = factor(gene, levels = gene_order)) %>%
    ggplot(aes(x = gene, y = cluster, fill = expr)) +
    geom_tile(show.legend = FALSE) +
    scale_fill_gradient2(low = "#A50026", mid = "#FFFFBF", high = "#313695") +
    theme_void() +
    theme(axis.text.y = element_text())
  
  # add score plot
  score_plot <-
    scores %>%
    group_by(cluster) %>%
    mutate("rank" = seq(n())) %>%
    slice(1:3) %>%
    ungroup() %>%
    complete(cluster) %>%
    mutate(cluster = factor(cluster, levels = cluster_order)) %>%
    ggplot(aes(x = rank, y = cluster, size = score, label = area)) +
    geom_text(show.legend = FALSE, hjust = "inward") +
    scale_x_reverse() +
    theme_void() +
    theme(plot.margin = unit(c(0,1,0,1), "lines"))

  cowplot::plot_grid(score_plot, tiles, dendro, 
                     ncol = 3,
                     rel_widths = c(0.3, 0.6, 0.1))
}
