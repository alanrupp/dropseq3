library(Seurat)
library(tidyverse)

# - Grab cell names by cluster membership or dismembership -------------------
get_cells <- function(object, cluster, not = FALSE) {
  if (not) {
    names(object@active.ident)[object@active.ident != cluster]
  } else {
    names(object@active.ident)[object@active.ident == cluster]
  }
}

# - Cluster means -------------------------------------------------------------
cluster_means <- function(object, genes = NULL, assay = "RNA", 
                          data_slot = "scale.data") {
  clusters <- sort(unique(object@active.ident))
  mtx <- slot(slot(object, "assays")[[assay]], data_slot)
  if (is.null(genes)) {
    genes <- rownames(mtx)
  } else {
    genes <- genes[genes %in% rownames(mtx)]
  }
  df <-
    map(clusters, ~ Matrix::rowMeans(mtx[genes, get_cells(object, .x)])) %>% 
    set_names(clusters) %>%
    bind_cols() %>%
    as.data.frame()
  rownames(df) <- genes
  return(df)
}

# - Order clusters by gene list -----------------------------------------------
order_clusters <- function(object, genes = NULL, scale = TRUE) {
  mean_values <- cluster_means(object, genes)
  tree <- hclust(dist(t(mean_values)))
  return(tree$labels[tree$order])
}

# - Choose PCs --------------------------------------------------------------
pca <- function(object) {
  # run PCA (adding 50 PCs each time) to capture all significant PCs
  pcs <- 0
  while (TRUE) {
    pcs <- pcs + 50
    object <- RunPCA(object, verbose = FALSE, npcs = pcs)
    object <- JackStraw(object, dims = pcs)
    object <- ScoreJackStraw(object, dims = 1:pcs)
    insig <- which(
      object@reductions$pca@jackstraw@overall.p.values[,"Score"] >= 0.05
    )
    if (length(insig) != 0) {
      break()
    }
  }
  pcs <- insig[1] - 1
  print(paste("Using", pcs, "PCs"))
  object@reductions$pca@misc$sig_pcs <- pcs
  return(object)
}

# - Choose resolution --------------------------------------------------------
choose_resolution <- function(object, resolutions = NULL, 
                              assay = "integrated", seed = NA) {
  dims <- 1:object@reductions$pca@misc$sig_pcs
  
  # cluster at each resolution, finding mean silhouette width for each
  while (TRUE) {
    if (is.null(resolutions)) {
      return(NULL)
    }
    print(paste("Finding clusters for resolutions", 
            paste(resolutions, collapse = ", ")))
    clusters <- 
      map(resolutions, 
          ~ FindClusters(object, resolution = .x, 
                         verbose = FALSE)@active.ident)
    distances <- dist(object@reductions$pca@cell.embeddings[, dims])
    
    # calc silhouettes for all resolutions
    calc_silhouettes <- function(cluster_results, assay = "integrated") {
      sil <- cluster::silhouette(as.numeric(cluster_results), distances)
      return(mean(sil[, 3]))
    }
    widths <- map_dbl(clusters, calc_silhouettes)
    names(widths) <- resolutions
    
    # choose resolution with highest mean silhouette width
    best_width <- sort(widths, decreasing = TRUE)[1] %>% names()
    if (best_width != max(resolutions)) {
      break()
    } else {
      resolutions <- seq(max(resolutions), max(resolutions) + 1, by = 0.2)
    }
  }
  print(paste("Using resolution", best_width, "for clustering."))
  return(as.numeric(best_width))
}

# - Cluster -------------------------------------------------------------------
cluster <- function(object, assay = "RNA", seed = NA) {
  if (assay == "integrated") {
    DefaultAssay(object) <- "integrated"
  } else {
    DefaultAssay(object) <- "RNA"
  }
  pcs <- object@reductions$pca@misc$sig_pcs
  # Find neighbors and cluster and different resolutions
  print("Finding nearest neighbors")
  object <- FindNeighbors(object, dims = 1:pcs, verbose = FALSE)
  res <- choose_resolution(object, assay = assay,
                           resolutions = seq(0.2, 1, by = 0.2))
  object <- FindClusters(object, resolution = res, random.seed = seed)
  gc(verbose = FALSE)
  return(object)
}

# - Optimize UMAP ------------------------------------------------------------
# Run UMAP for different neighbors and distance parameters
get_umap_coordinates <- function(mtx, neighbors, distance) {
  as.data.frame(uwot::umap(mtx, n_neighbors = neighbors, min_dist = distance))
}

# Get Dunn Index from embeddings and clusters
get_dunn_index <- function(embeddings, clusters) {
  distances <- dist(embeddings) %>% as.matrix()
  result <- map_dbl(
    clusters, 
    ~ min(distances[clusters == .x, clusters != .x]) / 
      max(distances[clusters == .x, clusters == .x])
  )
  return(mean(result))
}

# Get Silhouette width from embeddings and clusters
get_silhouette_width <- function(embeddings, clusters) {
  distances <- dist(embeddings)
  sil <- cluster::silhouette(clusters, distances)
  return(mean(sil[, 3]))
}

# Find optimized UMAP parameters
optimize_umap <- function(object, method = "dunn") {
  # set up parameters
  cat("Setting up parameter space\n")
  pcs <- object@reductions$pca@misc$sig_pcs
  mtx <- object@reductions$pca@cell.embeddings[, 1:pcs]
  neighbors <- seq(5, 50, by = 5)
  dists <- c(0.001, 0.005, 0.01, 0.05, 0.1, 0.5)
  search_grid <- expand.grid(neighbors, dists)
  
  # get UMAP coordinates
  cat("Calculating UMAP dimensions for each set of UMAP parameters\n")
  results <- apply(search_grid, 1, function(x) 
    get_umap_coordinates(mtx, x["Var1"], x["Var2"])
    )
  
  # find correlation between umap result and distance
  cat("Finding UMAP parameters that are highly correlated with PC distance\n")
  pc_distances <- dist(mtx) %>% as.matrix() %>% as.vector()
  correlations <- map_dbl(
    results, ~ cor(pc_distances, dist(.x) %>% as.matrix %>% as.vector)
    )
  
  cat("Choosing result based on highest correlation\n")
  search_grid$score <- correlations
  search_grid <- rename(search_grid, "n_neighbors" = Var1, "min_dist" = Var2)
  search_grid <- arrange(search_grid, desc(score))
  return(search_grid)
}

# - Run UMAP ------------------------------------------------------------------
run_umap <- function(object, method = "dunn") {
  parameters <- optimize_umap(object, method = method)
  object <- RunUMAP(object, 
                    dims = 1:object@reductions$pca@misc$sig_pcs,
                    n.neighbors = parameters[1, ]$n_neighbors,
                    min.dist = parameters[1, ]$min_dist,
                    verbose = FALSE)
  return(object)
}


# find highly expressed genes -----------------------------------------------
find_expressed <- function(object, min_cells = 10, min_counts = 1) {
  calc <- rowSums(object@raw.data >= min_counts) >= min_cells
  calc <- calc[calc == TRUE]
  return(names(calc))
}

# - Find markers ------------------------------------------------------------
find_markers <- function(object, cluster, other = NULL, genes = NULL,
                         remove_insig = TRUE) {
  if (is.null(genes)) {
    genes <- rownames(object@assays$RNA@counts)
  }
  cells_in <- names(object@active.ident)[object@active.ident %in% cluster]
  if (is.null(other)) {
    print(paste("Calculating cluster", paste(cluster, collapse = ",")))
    cells_out <- names(object@active.ident)[!object@active.ident %in% cluster]
  } else {
    print(paste("Calculating cluster", paste(cluster, collapse = ","), 
                "vs.", paste(other, collapse = ",")))
    cells_out <- names(object@active.ident)[object@active.ident %in% other]
  }
  
  # function to calculate gene attributes
  gene_test <- function(gene) {
    c(wilcox.test(
      object@assays$RNA@data[gene, cells_in], 
      object@assays$RNA@data[gene, cells_out], 
      alternative = "greater")$p.value,
      round(log(mean(object@assays$RNA@data[gene, cells_in])) - 
              log(mean(object@assays$RNA@data[gene, cells_out])), 3),
      round(sum(object@assays$RNA@counts[gene, cells_in] > 0) / 
              length(cells_in), 3),
      round(sum(object@assays$RNA@counts[gene, cells_out] > 0) / 
              length(cells_out), 3)
    )
  }
  
  # run on all genes
  mtx <- pbapply::pbsapply(genes, gene_test)
  mtx <- t(mtx)
  
  # adjust P value and export
  mtx <- as.data.frame(mtx) %>% 
    set_names("p_val", "avg_logFC", "pct.1", "pct.2") %>%
    mutate("p_val_adj" = p.adjust(p_val, method = "BH"),
           "cluster" = paste(cluster, collapse = ", "),
           "gene" = genes) %>%
    arrange(p_val_adj)
  if (remove_insig) {
    mtx <- filter(mtx, p_val_adj < 0.05)
  }
  gc(verbose = FALSE)
  return(mtx)
}

# - Find All Markers --------------------------------------------------------
find_all_markers <- function(object, genes = NULL, remove_insig = TRUE) {
  if (is.null(genes)) {
    genes <- rownames(object@assays$RNA@counts)
  }
  clusters <- sort(unique(object@active.ident))
  result <- map(
    clusters, 
    ~ find_markers(object, .x, genes = genes, remove_insig = FALSE)) %>%
    bind_rows() %>%
    mutate(p_val_adj = p.adjust(p_val, method = "BH"))
  if (remove_insig) {
    result <- filter(result, p_val_adj < 0.05)
  }
  gc(verbose = FALSE)
  return(result)
}

# - Cluster identity all cells -----------------------------------------------
get_gsea_scores <- function(markers, only_classes = NULL) {
  classes <- read_csv("~/Programs/dropseq3/data/celltype_markers.csv", 
                      col_types = c("ccc"))
  
  markers <- markers %>%
    group_by(cluster) %>%
    filter(!duplicated(gene))
  
  # get each score
  get_each_gsea_score <- function(marker_genes, class) {
    genes <- filter(classes, cells == class)$gene
    other_genes <- filter(classes, cells != class)$gene
    marker_genes <- marker_genes %>%
      mutate(score = ifelse(gene %in% genes, avg_logFC,
                            ifelse(gene %in% other_genes, -avg_logFC, 0)))
    score <- cumsum(marker_genes$score)
    return(max(score, na.rm = TRUE))
  }
  
  # for each cluster
  unique_classes <- sort(unique(classes$cells))
  if (!is.null(only_classes)) {
    unique_classes <- unique_classes[unique_classes %in% only_classes]
  }
  get_cluster_gsea_scores <- function(clstr) {
    marker_genes <- filter(markers, cluster == clstr) %>%
      filter(p_val_adj < 0.05) %>%
      arrange(desc(avg_logFC))
    scores <- map_dbl(unique_classes, ~ get_each_gsea_score(marker_genes, .x))
    # normalize
    scores <- scores - min(scores)
    scores <- scores / max(scores)
    return(scores)
  }
  
  # for all clusters
  clusters <- sort(unique(markers$cluster))
  result <- map(clusters, get_cluster_gsea_scores)
  result <- bind_cols(result) %>% set_names(clusters) %>% as.data.frame()
  rownames(result) <- unique_classes
  return(result)
}

# - Simplify conserved markers -----------------------------------------------
simplify_conserved_markers <- function(markers) {
  if (!"gene" %in% colnames(markers)) {
    markers <- rownames_to_column(markers, "gene")
  }
  
  # run logitp function from metap package to combine p values
  combine_p <- function(markers) {
    p_vals <- select(markers, ends_with("p_val_adj"))
    p_vals <- as.data.frame(p_vals)
    # convert any 1 to 1-10^-9
    p_vals <- sapply(p_vals, function(x) ifelse(x == 1, 1-10^-6, x))
    p_vals <- apply(p_vals, 1, function(x) metap::logitp(x)$p)
    return(p_vals)
  }
  p_vals <- combine_p(markers)
  
  # calculate pct cells expressing gene and fold change
  pct1 <- rowMeans(select(markers, ends_with("pct.1")))
  pct2 <- rowMeans(select(markers, ends_with("pct.2")))
  fc <- rowMeans(select(markers, ends_with("avg_logFC")))
  
  # select only relevant columns
  if ("cluster" %in% colnames(markers)) {
    markers <- select(markers, cluster, gene)
  } else {
    markers <- select(markers, gene)
  }
  markers <- markers %>% mutate(
    "avg_logFC" = fc,
    "pct.1" = pct1,
    "pct.2" = pct2,
    "p_val_adj" = p_vals)
  return(markers)
}

# - Find all conserved markers -----------------------------------------------
FindAllConservedMarkers <- function(object, ident2 = NULL, 
                                    groupby = NULL, clusters = NULL,
                                    verbose = FALSE) {
  if (is.null(clusters)) {
    clusters <- sort(unique(object@active.ident))
  }
  
  # if there are < 3 cells per cluster per sample, run standard FindMarkers
  too_few <- table(object@active.ident, object@meta.data[,groupby]) %>%
    as.data.frame() %>%
    filter(Freq < 3) %>%
    .$Var1
  too_few <- too_few[too_few %in% clusters]
  if (length(too_few) > 0) {
    few_markers <- map(
      too_few, ~ FindMarkers(object, .x, only.pos = TRUE, verbose = verbose)
      )
    few_markers <- map(few_markers, as.data.frame)
    few_markers <- map(few_markers, ~ rownames_to_column(.x, "gene"))
    few_markers <- map(seq(length(few_markers)),
                       ~ mutate(few_markers[[.x]], "cluster" = too_few[.x]))
    few_markers <- bind_rows(few_markers)
    few_markers <- select(few_markers, -p_val)
    few_markers <- mutate(few_markers, "note" = 'standard')
    clusters <- clusters[!clusters %in% too_few]
  }
  
  # Run FindConservedMarkers on each cluster
  markers <- map(clusters, 
                 ~ FindConservedMarkers(object, .x, ident.2 = ident2,
                                        grouping.var = groupby,
                                        only.pos = TRUE,
                                        verbose = verbose
                                        ))
  markers <- map(markers, as.data.frame)
  markers <- map(markers, ~ rownames_to_column(.x, "gene"))
  markers <- map(clusters,
                 ~ mutate(markers[[which(clusters == .x)]], "cluster" = .x))
  markers <- bind_rows(markers)
  
  # run logitp function from metap package to combine p values
  markers <- simplify_conserved_markers(markers)
  
  # add standard FindMarkers results
  if (length(too_few) > 0) {
    markers <- bind_rows(markers, few_markers)
    markers <- mutate(markers, note = ifelse(is.na(note), 'conserved', note))
  }
  
  # filter, arrange, and return final result
  markers <- filter(markers, p_val_adj < 0.05)
  markers <- markers %>%
    group_by(cluster) %>%
    arrange(desc(pct.1 - pct.2))
  return(markers)
}


# - Summarize markers -------------------------------------------------------
summarize_markers <- function(markers) {
  df <- markers %>%
    group_by(cluster) %>%
    summarize(total = sum(p_val_adj < 0.05)) %>%
    complete(cluster) %>%
    mutate(total = ifelse(is.na(total), 0, total))
  if (length(unique(markers$note)) > 1) {
    df <- left_join(df, select(filter(markers, !duplicated(cluster)), 
                               cluster, note), by = "cluster")
  }
  return(df)
}


# - Merge markerless --------------------------------------------------------
# merge markerless clusters with nearest cluster by correlation
merge_markerless <- function(object, markers) {
  marker_summary <- summarize_markers(markers)
  markerless <- filter(marker_summary, total == 0)$cluster
  
  # print clusters that don't have markers
  if (length(markerless) >= 1) {
    print(paste("Cluster(s)", paste(markerless, collapse = ", "), 
                "have no significantly enriched genes. Merging with neighbors or dropping. ")
    )
  } else {
    return(object)
  }
  
  # merge markerless clusters with nearest neighbor
  clusters <- unique(object@active.ident)
  mean_pcs <- map(clusters,
                  ~ colMeans(object@reductions$pca@cell.embeddings[
                    names(object@active.ident)[object@active.ident == .x],
                    1:object@reductions$pca@misc$sig_pcs
                    ])) %>%
    bind_cols() %>%
    set_names(clusters) %>%
    as.data.frame()
  neighbors <- cor(mean_pcs) %>%
    as.data.frame() %>%
    rownames_to_column("cluster") %>%
    gather(-cluster, key = "neighbor", value = "corr") %>%
    filter(cluster != neighbor) %>%
    group_by(cluster) %>%
    arrange(desc(corr)) %>%
    slice(1) %>%
    ungroup()
  
  neighbors <- neighbors %>% mutate(
    "cluster" = factor(cluster, levels = levels(object@active.ident)),
    "neighbor" = factor(neighbor, levels = levels(object@active.ident))
  )
  markerless_neighbors <- filter(neighbors, cluster %in% markerless)

  # merge or drop cells based on distance to nearest neighbor
  merge_cells <- function(cluster, neighbor) {
    print(paste("Cluster", cluster, "and cluster", neighbor,
              "are neighbors.",
              "Merging cluster", cluster, "into cluster", neighbor, '.'))
    object@active.ident[object@active.ident == cluster] <- neighbor
    return(object)
  }
  drop_cells <- function(cluster) {
    print(paste("Dropping cluster", cluster, "because it has no clear other",
              "cluster to merge into."))
    object <- subset(
      object, cells = names(object@active.ident)[object@active.ident != cluster]
      )
    return(object)
  }
  mean_corr <- mean(neighbors$corr)
  for (i in 1:nrow(markerless_neighbors)) {
    # if both umap and corr are the same, merge clusters
    if (markerless_neighbors[i, ]$corr > mean_corr) {
      object <- merge_cells(markerless_neighbors[i, ]$cluster, 
                            markerless_neighbors[i, ]$neighbor)
    } else {
      object <- drop_cells(markerless_neighbors[i, ]$cluster)
    }
  }
  # reset factor levels
  object@active.ident <- factor(object@active.ident, 
                                levels = sort(unique(object@active.ident)))
  return(object)
}

# - Find unique genes -------------------------------------------------------
find_unique_genes <- function(object, genes = NULL, clusters = NULL,
                              top_n = 1) {
  mtx <- object@assays$RNA@scale.data
  if (is.null(clusters)) {
    clusters <- sort(unique(object@active.ident))
  }
  if (is.null(genes)) {
    genes <- rownames(mtx)
  }
  
  cluster_expression <- function(mtx) {
    df <-
      map(clusters,
          ~ Matrix::rowMeans(mtx[genes, get_cells(object, .x)] > 0)
      ) %>% 
      set_names(clusters)
    return(df)
  }
  df <- cluster_expression(mtx) %>% bind_cols()
  
  # top 1 and 2 clusters for each gene
  first <- apply(df, 1, function(x) names(sort(x, decreasing = TRUE)[1]))
  second <- apply(df, 1, function(x) names(sort(x, decreasing = TRUE)[2]))
  diff <- apply(df, 1, function(x) sort(x, decreasing = TRUE)[1] - 
                  sort(x, decreasing = TRUE)[2])
  
  # result
  result <- data.frame(
    "gene" = genes,
    "first" = first,
    "second" = second,
    "diff" = diff
  )
  
  # set factor levels
  result <- result %>% 
    mutate_at(vars(first, second), ~ factor(.x, levels = levels(object@active.ident)))
  
  # keep top n
  result <- result %>%
    arrange(desc(diff)) %>%
    group_by(first) %>%
    slice(1:top_n) %>%
    ungroup() %>%
    arrange(first)
  
  return(result)
}

# - Eigengene ---------------------------------------------------------------
eigengene <- function(object, genes) {
  genes <- genes[genes %in% rownames(object@assays$RNA@scale.data)]
  pca <- prcomp(t(object@assays$RNA@scale.data[genes, ]))
  return(pca$x[,1])
}

# - deScore -------------------------------------------------------------------
# calculate "deScore" as in Tasic et al. 2018
deScore <- function(object, ident1, ident2) {
  n_samples <- table(object@active.ident, object$mouse)
  n_samples <- as.data.frame(n_samples) %>% filter(Var1 %in% c(ident1, ident2))
  if (sum(n_samples$Freq < 3) > 0) {
    result <- FindMarkers(object, ident.1 = ident1, ident.2 = ident2,
                          verbose = FALSE)
  } else {
    result <- FindConservedMarkers(object, ident.1 = ident1,
                                   ident.2 = ident2, grouping.var = "mouse",
                                   verbose = FALSE)
    result <- simplify_conserved_markers(result)
  }
  result <- mutate(result, p_val_adj = -log10(p_val_adj))
  result <- mutate(result, p_val_adj = ifelse(p_val_adj > 20, 20, p_val_adj))
  result <- sum(result$p_val_adj)
  return(result)
}

# - DE genes ------------------------------------------------------------------
# calculate total DE genes between 2 clusters
de_genes <- function(object, ident1, ident2) {
  n_samples <- table(object@active.ident, object$mouse)
  n_samples <- as.data.frame(n_samples) %>% filter(Var1 %in% c(ident1, ident2))
  if (sum(n_samples$Freq < 3) > 0) {
    result <- FindMarkers(object, ident.1 = ident1, ident.2 = ident2,
                          verbose = FALSE)
  } else {
    result <- FindConservedMarkers(object, ident.1 = ident1,
                                   ident.2 = ident2, grouping.var = "mouse",
                                   verbose = FALSE)
    result <- simplify_conserved_markers(result)
  }
  return(sum(result$p_val_adj < 0.05, na.rm = TRUE))
}

# - Merging clusters --------------------------------------------------------
merge_clusters <- function(object, markers) {
  while (TRUE) {
    # find 2 neighbors for every cluster based UMAP distance
    clusters <- unique(object@active.ident)
    umap <- object@reductions$umap@cell.embeddings
    cluster_means <- umap %>% 
      mutate("cluster" = object@active.ident) %>%
      group_by(cluster) %>%
      summarize("UMAP1" = median(UMAP_1), "UMAP2" = median(UMAP_2)) %>%
      as.data.frame() %>%
      column_to_rownames(cluster)
    distances <- dist(cluster_means) %>% as.matrix() %>%
      as.data.frame() %>%
      rownames_to_column("cluster") %>%
      gather(-cluster, key = "neighbor", value = "dist") %>%
      filter(cluster != neighbor)
    neighbors <- distances %>%
      group_by(cluster) %>%
      arrange(dist) %>%
      slice(1:2) %>%
      ungroup()
    
    # calculate deScore for all neighbors and merge most similar
    neighbors$deScore <- apply(neighbors, 1, 
                               function(x) de_genes(object, x[1], x[2]))
    neighbors <- arrange(neighbors, deScore)
    neighbors <- neighbors %>% mutate(
      cluster = factor(cluster, levels = levels(object@active.ident)),
      neighbor = factor(neighbor, levels = levels(object@active.ident))
    )
    
    if (sum(neighbors$deScore < 10) == 0) {
      break()
    } else {
      # merge 2 most similar clusters
      to_merge <- slice(neighbors, 1) %>% as.data.frame()
      print(paste("Merging clusters", to_merge$cluster, "and", to_merge$neighbor,
                  "| de genes:", as.integer(to_merge$deScore)))
      object@active.ident[object@active.ident == to_merge$neighbor] <-
        to_merge$cluster
    }
  }
  return(object)
}


# - Finding doublets ---------------------------------------------------------
# Finding doublets based on the Tasic et al. Nature 2018 criteria
find_doublets <- function(object, markers) {
  # find eigengene for each cell for each set of cluster markers
  eigengenes <-
    map(unique(object@active.ident),
        ~ eigengene(object, filter(markers, cluster == .x)$gene))
  eigengenes <- map(eigengenes, scale)
  members <- map(eigengenes, ~ .x[.x > 3, ])
  members <- unlist(members)
  doublets <- names(members)[duplicated(members)]
  return(doublets)
}

# - Removing clusters of mostly doublets ------------------------------------
find_doublet_clusters <- function(object, doublets) {
  if (length(doublets) == 0) {
    return(NULL)
  }
  # find expected doublet rate based on cells per sample
  doublet_rates <- read_csv("~/Programs/dropseq3/data/10x_doublet_rate.csv")
  model <- lm(rate ~ cells, data = doublet_rates)
  expected_doublets <- predict(model,
    data.frame("cells" = ncol(object@assays$RNA@data)/
                 length(unique(object$mouse))))
  
  # find actual doublet rates
  cluster_rates <- table(object@active.ident, 
                         names(object@active.ident) %in% doublets
                         )
  above_threshold <- as.data.frame(cluster_rates) %>%
    group_by(Var1) %>%
    summarize("Rate" = Freq[Var2 == TRUE]/sum(Freq)) %>%
    filter(Rate > expected_doublets)
  
  # if clusters have rates above expected, run a statistical test
  if (nrow(above_threshold) == 0) {
    cat("No clusters have more doublets than expected")
    return(NULL)
  } else {
    cluster_rates <- as.matrix(cluster_rates)
    cluster_rates <- cluster_rates[above_threshold$Var1, ]
    result <- map(cluster_rates, 
                  ~ prop.test(.x, p = expected_doublets)$p.value)
    result <- unlist(result)
    names(result) <- rownames(cluster_rates)
    result <- result[result < 0.05]
    if (length(result > 0)) {
      cat(paste(length(result), "clusters have higher frequency of doublets than expected"))
      return(names(result))
    } else {
      cat("No clusters have a significantly elevated frequency of doublets")
      return(NULL)
    }
  }
}


# - edgeR for treatment effects -----------------------------------------------
edgeR_test <- function(object, treatment) {
  # get data
  mtx <- object@assays$RNA@counts
  clusters <- object@active.ident
  tx <- paste(object@meta.data[, treatment], clusters, sep = "_")
  
  # set up edgeR
  library(edgeR)
  y <- DGEList(counts = mtx, group = tx)
  y <- calcNormFactors(y)
  design <- model.matrix(~ 0 + tx)
  colnames(design) <- levels(tx)
  y <- estimateDisp(y, design)
  
  # make contrasts
  contrast_list <- c()
  contrasts <- makeContrasts(
    "0" = CNO_0 - Saline_0,
    "1" = CNO_1 - Saline_1,
    "2" = CNO_2 - Saline_2,
    "3" = CNO_3 - Saline_3,
    "4" = CNO_4 - Saline_4,
    "5" = CNO_5 - Saline_5,
    "6" = CNO_6 - Saline_6,
    "7" = CNO_7 - Saline_7,
    "8" = CNO_8 - Saline_8,
    "9" = CNO_9 - Saline_9,
    "10" = CNO_10 - Saline_10,
    "11" = CNO_11 - Saline_11,
    levels = design
  )
  
  # perform LRT
  fit <- glmFit(y, design)
  de <- function(coef) {
    lrt <- glmLRT(fit, contrast = contrasts[, coef])
    topTags(lrt, n = nrow(mtx))
  }
  result <- map(colnames(contrasts), de)
  names(result) <- colnames(contrasts)
  result <- map(result, as.data.frame)
  result <- map(result, ~ rownames_to_column(.x, "gene"))
  
  map(result, ~ sum(.x$FDR < 0.05, na.rm = TRUE))
  result <- map(seq(length(result)),
                ~ mutate(result[[.x]], "cluster" = names(result)[.x]))
  result <- bind_rows(result)
}

# - Grab clusters by type ------------------------------------------------------
grab_type <- function(object, classes, type) {
  # construct a tree based on marker genes for that cell type
  class_genes <- read_csv("~/Programs/dropseq3/data/celltype_markers.csv")
  means <- cluster_means(object, genes = class_genes$gene)
  tree <- hclust(dist(t(means)))
  
  # cut the tree at different levels and choose cut based on Youden's J
  youden <- function(cut) {
    new_clusters <- cutree(tree, cut)
    confusion <- table(new_clusters, classes$class == type) %>% as.data.frame()
    real_type <- confusion %>%
      filter(Var2 == TRUE) %>%
      filter(Freq == max(Freq)) %>%
      slice(1) %>%
      .$new_clusters
    tp <- confusion %>% filter(new_clusters == real_type & Var2 == TRUE) %>%
      .$Freq
    tn <- confusion %>% filter(new_clusters != real_type & Var2 == TRUE) %>%
      .$Freq %>% sum()
    fn <- confusion %>% filter(new_clusters != real_type & Var2 == FALSE) %>%
      .$Freq %>% sum()
    fp <- confusion %>% filter(new_clusters == real_type & Var2 == FALSE) %>%
      .$Freq
    j <- (tp/(tp+fn)) + (tn/(tn+fp)) - 1
    return(j)
  }
  js <- map_dbl(2:ncol(means), youden)
  names(js) <- 2:ncol(means)
  cut_to_use <- names(js)[js == max(js)][1]
  
  # return cluster IDs associated with that cut
  new_clusters <- cutree(tree, cut_to_use)
  confusion <- table(new_clusters, classes$class == type) %>% as.data.frame()
  real_type <- confusion %>%
    filter(Var2 == TRUE) %>%
    filter(Freq == max(Freq)) %>%
    slice(1) %>%
    .$new_clusters
  return(names(new_clusters)[new_clusters == real_type])
}

# - Co-clustering -----------------------------------------------------------
# recluster a subset of samples
cluster_sample <- function(object, iter, prop = 0.8) {
  set.seed(iter)
  all_cells <- names(object@active.ident)
  cells <- sample(all_cells, length(all_cells) * prop)
  object <- subset(object, cells = cells)
  object <- cluster(object, seed = iter)
  return(object@active.ident)
}

# run reclustering N times
cocluster_frequency <- function(object, iterations = 100) {
  clusters <- map(seq(iterations), ~ cluster_sample(object, .x))
  clusters <- map(clusters, ~ data.frame("cell" = names(.x), "cluster" = .x))
  clusters <- map(1:length(clusters), ~
                    mutate(clusters[[.x]], "iter" = .x)) %>% 
    bind_rows()
  write.csv(clusters, "clusters.csv", row.names = FALSE)
  # analyze results in python for speed
  source_python("/home/alanrupp/Programs/dropseq3/python/cocluster.py")
  d <- read_clusters("clusters.csv")
  scores <- score(d)
  rm("clusters.csv")
  return(scores)
}

# choose clusters based on co-clustering score
choose_coclustering_groups <- function(co, max_clusters = 50) {
  distances <- dist(co)
  tree <- hclust(distances)
  distances <- as.matrix(distances)
  
  # get J scores
  calc_youden <- function(cut) {
    clusters <- cutree(tree, cut)
    tp <- sum(map_dbl(clusters, ~ sum(co[clusters == .x, clusters == .x] > 0)))
    tn <- sum(map_dbl(clusters, ~ sum(co[clusters == .x, clusters != .x] == 0)))
    fp <- sum(map_dbl(clusters, ~ sum(co[clusters == .x, clusters == .x] == 0)))
    fn <- sum(map_dbl(clusters, ~ sum(co[clusters == .x, clusters != .x] > 0)))
    j <- (tp/(tp+fn)) + (tn/(tn+fp)) - 1
    return(j)
  }
  cuts <- seq(2, max_clusters)
  j <- map_dbl(cuts, calc_youden)
  names(j) <- cuts
  
  # choose clustering with highest J
  clusters <- cutree(tree, as.integer(names(j)[j == max(j)]))
  return(clusters)
}

# name clusters
name_clusters <- function(object, markers) {
  if (length(unique(markers$cluster)) != length(unique(object@active.ident))) {
    stop("Mismatched Seurat object and markers data.frame")
  }
  
}

# - Traverse tree ------------------------------------------------------------
traverse_tree <- function(object, genes = NULL) {
  
}

# - Calculate p_val, avg_logFC, pct.1, and pct.2 for a gene -------------------
gene_test <- function(object, gene, cells_in, cells_out) {
  if (length(cells_in) == 0) {
    warning("Not enough cells")
    return(c("p_val" = NA, "avg_logFC" = NA, "pct.1" = NA, "pct.2" = NA))
  }
  if (length(cells_out) == 0) {
    warning("Not enough cells")
    return(c("p_val" = NA, "avg_logFC" = NA, "pct.1" = NA, "pct.2" = NA))
  }
  p_val <- wilcox.test(
    object@assays$RNA@data[gene, cells_in], 
    object@assays$RNA@data[gene, cells_out], 
    alternative = "greater")$p.value
  if (p_val == 1) {
    p_val <- 1-10^-9
  }
  c("p_val" = p_val,
    "avg_logFC" = round(log(mean(object@assays$RNA@data[gene, cells_in])) -
                        log(mean(object@assays$RNA@data[gene, cells_out])), 3),
    "pct.1" = round(sum(object@assays$RNA@counts[gene, cells_in] > 0) / 
                    length(cells_in), 3),
    "pct.2" = round(sum(object@assays$RNA@counts[gene, cells_out] > 0) / 
                    length(cells_out), 3)
  )
}

# - Find conserved markers ---------------------------------------------------
find_conserved_markers <- function(object, cluster, groupby, 
                                   other = NULL, remove_insig = TRUE,
                                   genes = NULL, progress_bar = TRUE) {
  groups <- unique(object@meta.data[, groupby])
  if (is.null(genes)) {
    genes <- rownames(object@assays$RNA@counts)
  }
  print(paste0("Testing cluster ", cluster, ": ", 
               length(genes), " genes from ",
               sum(object@active.ident == cluster), " cells across ",
               length(groups), " groups (", 
               format(Sys.time(), '%H:%M'), ")"))
  
  
  # functions to grab cells by cluster and group membership
  get_cells <- function(group) {
    names(object@active.ident)[object@active.ident == cluster & 
                                 object@meta.data[,groupby] == group]
  }
  get_other <- function(group) {
    if (is.null(other)) {
      names(object@active.ident)[object@active.ident != cluster & 
                                 object@meta.data[, groupby] == group]
    } else {
      names(object@active.ident)[object@active.ident %in% other & 
                                 object@meta.data[, groupby] == group]
    }
  }
  
  # run on all genes
  group_test <- function(gene) {
    result <- sapply(groups, function(x) 
      gene_test(object, gene, get_cells(x), get_other(x))
    )
    if (sum(is.na(result[1,]) == length(result[1,]))) {
      return(c("p_val" = NA, rowMeans(result[2:4,], na.rm = TRUE)))
    } else {
      return(c("p_val" = metap::logitp(result[1, ][!is.na(result[1, ])])$p,
               rowMeans(result[2:4,], na.rm = TRUE)))
    }
  }
  if (progress_bar) {
    mtx <- pbapply::pbsapply(genes, group_test)
  } else {
    mtx <- sapply(genes, group_test)
  }
  mtx <- t(mtx)
  
  # adjust P value and export
  mtx <- as.data.frame(mtx) %>% 
    mutate("p_val_adj" = p.adjust(p_val, method = "BH"),
           "cluster" = paste(cluster, collapse = ", "),
           "gene" = genes) %>%
    arrange(p_val_adj)
  if (remove_insig) {
    mtx <- filter(mtx, p_val_adj < 0.05)
  }
  gc(verbose = FALSE)
  return(mtx)
}

# - Find all conserved markers -----------------------------------------------
find_all_conserved_markers <- function(object, groupby, remove_insig = TRUE,
                                       genes = NULL, progress_bar = FALSE) {
  if (is.null(genes)) {
    genes <- rownames(object@assays$RNA@counts)
  }
  clusters <- sort(unique(object@active.ident))
  result <- map(
    clusters, 
    ~ find_conserved_markers(object, .x, groupby, 
                             genes = genes, remove_insig = FALSE,
                             progress_bar = progress_bar)) %>%
    bind_rows() %>%
    mutate(p_val_adj = p.adjust(p_val, method = "BH"))
  if (remove_insig) {
    result <- filter(result, p_val_adj < 0.05)
  }
  gc(verbose = FALSE)
  return(result)
}
