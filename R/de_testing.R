library(Seurat)
library(tidyverse)
library(reticulate)

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
                          data = "scale.data") {
  clusters <- sort(unique(object@active.ident))
  mtx <- slot(object@assays[[assay]], data)
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
order_clusters <- function(object, genes = NULL, max_expr = 2, scale = TRUE) {
  mean_values <- cluster_means(object, genes)
  mean_values <- clip_matrix(mean_values, max_expr)
  if (scale) {
    tree <- hclust(dist(t(mean_values)))
  } else {
    tree <- hclust(dist(mean_values))
  }
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
      object@reductions$pca@jackstraw@overall.p.values[, "Score"] >= 0.05
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

# - Choose k neighbors -------------------------------------------------------
choose_neighbors <- function(object, klim = c(20, NA), assay = "RNA") {
  if (is.na(klim[2])) { klim[2] <- round(
    sqrt(ncol(slot(object@assays[[assay]], "data"))), 0) }
  if (klim[2] < klim[1]) { klim[1] <- klim[2] }
  k_vals <- seq(klim[1], klim[2], by = 10)
  pcs <- object@reductions$pca@misc$sig_pcs
  pipeline <- function(k) {
    object <- FindNeighbors(object, dims = 1:pcs, verbose = FALSE, k.param = k,
                            assay = assay)
    FindClusters(object, resolution = 1, verbose = FALSE)@active.ident
  }
  clusters <- map(k_vals, pipeline)
  # calc silhouettes for all k values
  distances <- dist(object@reductions$pca@cell.embeddings[, 1:pcs])
  calc_silhouettes <- function(cluster_results) {
    sil <- cluster::silhouette(as.numeric(cluster_results), distances)
    return(mean(sil[, 3]))
  }
  widths <- map_dbl(clusters, calc_silhouettes)
  gc(verbose = FALSE)
  return(k_vals[widths == max(widths)])
}


# - Choose resolution --------------------------------------------------------
choose_resolution <- function(object, resolutions = NULL, 
                              assay = "RNA", seed = NA) {
  dims <- 1:object@reductions$pca@misc$sig_pcs
  print("Calculating distance matrix ...")
  distances <- dist(object@reductions$pca@cell.embeddings[, dims])
  
  # cluster at each resolution, finding mean silhouette width for each
  while (TRUE) {
    if (is.null(resolutions)) { return(NULL) }
    print(paste("Finding clusters for resolutions", 
            paste(resolutions, collapse = ", ")))
    clusters <- map(
      resolutions, 
      ~ FindClusters(object, resolution = .x, verbose = FALSE)@active.ident
      )
    
    # calc silhouettes for all resolutions
    calc_silhouettes <- function(cluster_results) {
      sil <- cluster::silhouette(as.numeric(cluster_results), distances)
      return(mean(sil[, 3]))
    }
    widths <- map_dbl(clusters, calc_silhouettes)
    names(widths) <- resolutions
    
    # choose resolution with highest mean silhouette width
    best_width <- sort(widths, decreasing = TRUE)[1] %>% names()
    # if best resolution is the max resolution tested, increase resolution
    if (best_width != max(resolutions)) {
      break()
    } else {
      resolutions <- seq(max(resolutions), max(resolutions) + 1, by = 0.2)
    }
    gc(verbose = FALSE)
  }
  # return clusters from max resolution 
  print(paste("Using resolution", best_width, "for clustering."))
  return(clusters[[which(resolutions == best_width)]])
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
  print(paste("Choosing optimal k for k-nearest neighbors"))
  k <- choose_neighbors(object)
  print(paste("Finding nearest neighbors using k =", k))
  object <- FindNeighbors(object, dims = 1:pcs, verbose = FALSE, k.param = k)
  clusters <- choose_resolution(object, assay = assay,
                                resolutions = seq(0.2, 1, by = 0.2))
  object@active.ident <- clusters
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
                         remove_insig = TRUE, progress_bar = TRUE,
                         only_pos = TRUE) {
  if (is.null(genes)) { genes <- rownames(object@assays$RNA@counts) }
  cells_in <- names(object@active.ident)[object@active.ident %in% cluster]
  if (is.null(other)) {
    print(paste("Finding markers for cluster", paste(cluster, collapse = ",")))
    cells_out <- names(object@active.ident)[!object@active.ident %in% cluster]
  } else {
    print(paste("Finding markers for cluster", paste(cluster, collapse = ","), 
                "vs.", paste(other, collapse = ",")))
    cells_out <- names(object@active.ident)[object@active.ident %in% other]
  }
  
  # function to calculate gene attributes
  gene_test <- function(gene) {
    result <- c(
      "avg_logFC" = round(log(mean(object@assays$RNA@data[gene, cells_in])) -
                          log(mean(object@assays$RNA@data[gene, cells_out])), 
                          3),
      "pct.1" = round(sum(object@assays$RNA@counts[gene, cells_in] > 0) / 
                        length(cells_in), 3),
      "pct.2" = round(sum(object@assays$RNA@counts[gene, cells_out] > 0) /
                        length(cells_out), 3)
    )
    if (only_pos) {
      if (result["avg_logFC"] > 0 & result["pct.1"] > result["pct.2"]) {
        p_val <- wilcox.test(
          object@assays$RNA@data[gene, cells_in], 
          object@assays$RNA@data[gene, cells_out], 
          alternative = "greater")$p.value
      } else {
        p_val <- 1-10^-9
      }
    } else {
      p_val <- wilcox.test(
        object@assays$RNA@data[gene, cells_in], 
        object@assays$RNA@data[gene, cells_out], 
        alternative = "greater")$p.value
    }
    return(c("p_val" = p_val, result))
  }
  
  # run on all genes
  if (progress_bar) { mtx <- pbapply::pbsapply(genes, gene_test) }
  else { mtx <- sapply(genes, gene_test) }
  
  # adjust P value and export
  mtx <- as.data.frame(t(mtx)) %>% 
    mutate("p_val_adj" = p.adjust(p_val, method = "BH"),
           "cluster" = paste(cluster, collapse = ", "),
           "gene" = genes) %>%
    arrange(p_val_adj)
  if (remove_insig) { mtx <- filter(mtx, p_val_adj < 0.05) }
  gc(verbose = FALSE)
  return(mtx)
}

# - Find All Markers --------------------------------------------------------
find_all_markers <- function(object, genes = NULL, remove_insig = FALSE,
                             progress_bar = TRUE, only_pos = TRUE) {
  if (is.null(genes)) { genes <- rownames(object@assays$RNA@data) }
  clusters <- sort(unique(object@active.ident))
  result <- map(
    clusters, 
    ~ find_markers(object, .x, genes = genes, remove_insig = FALSE,
                   progress_bar = progress_bar, only_pos = only_pos)) %>%
    bind_rows() %>%
    mutate(p_val_adj = p.adjust(p_val, method = "BH"))
  if (remove_insig) { result <- filter(result, p_val_adj < 0.05) }
  gc(verbose = FALSE)
  return(result)
}

# - Cluster identity all cells -----------------------------------------------
get_gsea_scores <- function(markers, only_classes = NULL) {
  classes <- read_csv("~/Programs/dropseq3/data/celltype_markers.csv", 
                      col_types = c("ccc"))
  
  markers <- markers %>%
    group_by(cluster) %>%
    filter(!duplicated(gene)) %>%
    filter(pct.1 > pct.2)
  
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
    arrange(desc(pct.1 - pct.2)) %>%
    ungroup()
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


# - Eigengene ---------------------------------------------------------------
eigengene <- function(object, genes) {
  # keep genes that have been scaled and are not NA in the scale.data slot
  genes <- genes[genes %in% rownames(object@assays$RNA@scale.data)]
  genes <- genes[
    sapply(genes, function(x) sum(is.na(object@assays$RNA@scale.data[x, ])) == 0)
    ]
  pca <- prcomp(t(object@assays$RNA@scale.data[genes, ]))
  return(pca$x[, 1])
}

# - Finding doublets ---------------------------------------------------------
# Finding doublets based on the Tasic et al. Nature 2018 criteria
find_doublets <- function(object, markers, sd_threshold = 1) {
  # find eigengene for each cell for each set of cluster markers
  eigengenes <- map(
    sort(unique(object@active.ident)),
    ~ eigengene(object, filter(markers, cluster == .x)$gene)
    )
  # find cells that are >`threshold` SD above mean eigengene for >=2 clusters
  eigengenes <- map(eigengenes, scale)
  members <- map(eigengenes, ~ .x[.x > sd_threshold, ]) %>% unlist()
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

# - Calculate Youden's J ------------------------------------------------------
calc_youden <- function(contingency_table, return_all = FALSE) {
  TP <- contingency_table[2, 2]
  TN <- contingency_table[1, 1]
  FP <- contingency_table[2, 1]
  FN <- contingency_table[1, 2]
  J <- (TP/(TP+FN)) + (TN/(TN+FP)) - 1
  if (return_all) {
    return(c("TP" = TP, "TN" = TN, "FP" = FP, "FN" = FN, "J" = J))
  } else {
    return(c("J" = J))
  }
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

# - Check that a gene is in a given assay & dataset ---------------------------
check_gene <- function(object, gene, assay = "RNA", data = "data") {
  gene %in% rownames(slot(object@assays[[assay]], data))
}

# - Calculate p_val, avg_logFC, pct.1, and pct.2 for a gene -------------------
gene_test <- function(object, gene, cells_in, cells_out, test = "wilcox") {
  if (length(cells_in) == 0 | length(cells_out) == 0) {
    warning("Not enough cells")
    return(c("p_val" = NA, "avg_logFC" = NA, "pct.1" = NA, "pct.2" = NA))
  }
  if (test == "wilcox") {
    if (check_gene(object, gene)) {
      p_val <- wilcox.test(
        object@assays$RNA@data[gene, cells_in], 
        object@assays$RNA@data[gene, cells_out], 
        alternative = "greater")$p.value
    } else {
      warning(paste(gene, "not in dataset"))
      return(c("p_val" = NA, "avg_logFC" = NA, "pct.1" = NA, "pct.2" = NA))
    }
  }
  if (p_val == 1) { p_val <- 1-10^-9 }
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
  if (is.null(genes)) { genes <- rownames(object@assays$RNA@counts) }
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
    if (sum(is.na(result[1, ]) == length(result[1, ]))) {
      return(c("p_val" = NA, rowMeans(result[2:4,], na.rm = TRUE)))
    } else if (sum(!is.na(result[1,])) == 1) {
      return(c("p_val" = result[1,][!is.na(result[1,])], 
               rowMeans(result[2:4,], na.rm = TRUE)))
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
  if (remove_insig) { mtx <- filter(mtx, p_val_adj < 0.05) }
  gc(verbose = FALSE)
  return(mtx)
}

# - Find all conserved markers -----------------------------------------------
find_all_conserved_markers <- function(object, groupby, remove_insig = TRUE,
                                       genes = NULL, progress_bar = FALSE) {
  if (is.null(genes)) { genes <- rownames(object@assays$RNA@counts) }
  clusters <- sort(unique(object@active.ident))
  result <- map(
    clusters, 
    ~ find_conserved_markers(object, .x, groupby, 
                             genes = genes, remove_insig = FALSE,
                             progress_bar = progress_bar)) %>%
    bind_rows() %>%
    mutate(p_val_adj = p.adjust(p_val, method = "BH"))
  if (remove_insig) { result <- filter(result, p_val_adj < 0.05) }
  gc(verbose = FALSE)
  return(result)
}

# - CELLEX --------------------------------------------------------------------
cellex <- function(counts, clusters) {
  clusters <- data.frame("cluster" = clusters)
  source_python("/home/alanrupp/Programs/dropseq3/python/CELLEX.py")
  esmu <- run_cellex(counts, clusters)
  return(esmu)
}

# - Traverse tree -------------------------------------------------------------
traverse_tree <- function(object, genes = NULL, assay = "RNA",
                          data = "data", groupby = NULL,
                          max_expr = 2) {
  # grab mean data for each cluster
  if (is.null(genes)) { genes <- rownames(slot(object[[assay]], data)) }
  means <- cluster_means(object, genes, assay = assay, data = data)
  means <- clip_matrix(means, max_expr)
  # generate tree
  tree <- hclust(dist(t(means)))
  max_clusters <- length(unique(object@active.ident))
  # find max "DE-ness" for every cut
  # this is the max -log10 p val for every gene
  get_de <- function(cut) {
    print(paste("--- Testing", cut, "clusters of", max_clusters, "---"))
    clusters <- cutree(tree, cut)
    if (is.null(groupby)) {
      markers <- map(
        unique(clusters),
        ~ FindMarkers(object, 
                      ident.1 = names(clusters)[clusters == .x],
                      ident.2 = names(clusters)[clusters != .x],
                      features = genes, only.pos = TRUE, verbose = FALSE) %>%
          as.data.frame() %>%
          rownames_to_column("gene") %>%
          mutate("cluster" = .x)
      )
    } else {
      markers <- map(
        unique(clusters),
        ~ FindConservedMarkers(object, 
                               ident.1 = names(clusters)[clusters == .x],
                               ident.2 = names(clusters)[clusters != .x],
                               grouping.var = groupby,
                               features = genes, only.pos = TRUE,
                               meta.method = metap::logitp,
                               verbose = FALSE)  %>%
          as.data.frame() %>%
          rownames_to_column("gene") %>%
          mutate("cluster" = .x) %>%
          simplify_conserved_markers()
      )
    }
    markers <- bind_rows(markers)
    # keep only the lowest p value for each gene
    result <- markers %>% 
      filter(p_val_adj < 0.05) %>%
      arrange(desc(pct.1 - pct.2)) %>%
      group_by(cluster) %>%
      slice(1:20) %>% ungroup() %>%
      summarize("score" = median(pct.1 - pct.2) * 100) %>%
      .$score
    print(paste("Median pct de:", result))
    return(result)
  }
  # find DE score for all levels of tree
  result <- map_dbl(2:max_clusters, get_de)
  names(result) <- 2:max_clusters
  print(paste("Scores:", result))
  return(cutree(tree, names(result)[result == max(result)]))
}

# - Extend table --------------------------------------------------------------
# extend table that drops row/column of all 0 values
extend_table <- function(ct, n = 2) {
  to_add <- rep(0, n)
  if (ncol(ct) < n) {
    if ("TRUE" %in% colnames(ct)) {
      ct <- cbind(to_add, as.matrix(ct))
    } else {
      ct <- cbind(as.matrix(ct), to_add)
    }
  } else if (nrow(ct) < n) {
    if ("TRUE" %in% rownames(ct)) {
      ct <- rbind(to_add, as.matrix(ct))
    } else {
      ct <- rbind(as.matrix(ct), to_add)
    }
  }
  return(ct)
}

# - Get informative genes -----------------------------------------------------
get_informative_genes <- function(object, markers, n = 1, clusters = NULL,
                                  p_val_max = 0.05, n_genes = 50,
                                  other_pct = 0.2) {
  # filter the markers data.frame by user input to get top `n_genes` genes
  cluster_levels <- levels(markers$cluster)
  if (!is.null(clusters)) { markers <- filter(markers, cluster %in% clusters) }
  markers <- filter(markers, pct.1 > pct.2)
  markers <- filter(markers, pct.2 < other_pct)
  markers <- filter(markers, p_val_adj < p_val_max)
  markers <- markers %>% 
    group_by(cluster) %>%
    arrange(desc(pct.1 - pct.2)) %>% 
    dplyr::slice(1:n_genes)
  if (nrow(markers) == 0) { stop("Not enough markers") }
  # function to return gene TP, TN, FP, FN, and J score
  gene_test <- function(gene, cluster) {
    if (length(gene) == 1) {
      ct <- table(object@assays$RNA@counts[gene, ] > 0, 
                  object@active.ident == cluster)
    } else if (length(gene) > 1) {
      ct <- table(
        Matrix::colSums(object@assays$RNA@counts[gene, ] > 0) == length(gene), 
        object@active.ident == cluster
      )
    }
    # catch if one category has all zeros
    if (nrow(ct) != 2 | ncol(ct) != 2) { ct <- extend_table(ct, n = 2) }
    # calculate true/false positives and negatives
    TP <- ct[2, 2]; TN <- ct[1, 1]; FP <- ct[2, 1]; FN <- ct[1, 2]
    J <- (TP/(TP+FN)) + (TN/(TN+FP)) - 1
    return(data.frame(
      "gene" = paste(gene, collapse = ", "), "cluster" = cluster,
      "TP" = TP, "TN" = TN, "FP" = FP, "FN" = FN, "J" = J
    ))
  }
  # function to generate gene lists by cluster and n
  cluster_test <- function(clstr, n) {
    print(paste("Testing cluster", clstr, "|", n, "gene combinations"))
    genes <- expand.grid(rep(list(dplyr::filter(markers, cluster == clstr)$gene), n))
    if (n > 1) { # remove rows with identical values & duplicate values
      genes <- genes[apply(genes, 1, function(x) length(unique(x)) == length(x)), ]
      genes <- genes[!duplicated(t(apply(genes, 1, sort))), ]
    }
    if (nrow(genes) > 0) {
      result <- apply(genes, 1, function(x) gene_test(x, clstr))
      result <- map(result, ~ mutate(.x, gene = as.character(gene)))
      return(bind_rows(result))
    } else { 
      return(data.frame("gene" = NA, "cluster" = cluster,
                        "TP" = NA, "TN" = NA, "FP" = NA, "FN" = NA, "J" = NA
      )) 
    }
  }
  # run on all clusters for all n values
  result <- expand.grid(unique(markers$cluster), seq(n))
  result <- apply(result, 1, function(x) cluster_test(x["Var1"], x["Var2"]))
  result <- bind_rows(result)
  result <- dplyr::arrange(result, desc(J))
  result <- result %>% complete(cluster, fill = list("gene" = ""))
  return(result)
}

# - Name clusters -------------------------------------------------------------
name_clusters <- function(object, markers, other_pct = 0.2, n = 1,
                          only_annotated = FALSE) {
  # get informative genes based on TP/TN/FP/FN metrics
  cluster_levels <- levels(markers$cluster)
  info_genes <- get_informative_genes(object, markers, n = n,
                                      other_pct = other_pct)
  if (only_annotated) {
    annotation <- read_csv("~/Programs/dropseq3/data/annotation.csv")
    annotation <- dplyr::filter(annotation, rowSums(annotation[,2:ncol(annotation)]) > 0)
    info_genes <- dplyr::filter(info_genes, gene %in% annotation$gene)
  }
  info_genes <- mutate(info_genes, cluster = factor(cluster, levels = cluster_levels))
  info_genes <- info_genes %>% group_by(cluster) %>% dplyr::slice(1) %>% ungroup()
  info_genes <- select(info_genes, cluster, gene)
  info_genes <- mutate(info_genes, "name" = paste(cluster, gene, sep = "."))
  info_genes <- mutate(info_genes, name = factor(name, levels = name))
  if (n > 1) {
    info_genes <- info_genes %>% mutate(
      "name" = ifelse(str_detect(name, "\\,"), str_replace(name, "\\, ", "_"),
                      name))
  }
  return(info_genes)
}

# - Make pseudobulk -----------------------------------------------------------
make_pseudobulk <- function(object, treatment, batch) {
  clusters <- sort(unique(object@active.ident))
  get_cells <- function(tx, btch = NULL) {
    rownames(object@meta.data)[object@meta.data[, treatment] == tx &
                                 object@meta.data[, batch] == btch]
  }
  generate_matrix <- function(cluster) {
    groups <- expand.grid(unique(object@meta.data[, treatment]), 
                          unique(object@meta.data[, batch]))
    mtx <- apply(groups, 1, function(x) 
      Matrix::rowSums(object@assays$RNA@counts[ 
        , intersect(names(object@active.ident)[object@active.ident == cluster],
                    get_cells(x["Var1"], x["Var2"]))
        ])
    )
    colnames(mtx) <- paste(groups[,1], groups[,2], sep = "_")
    return(mtx)
  }
  mtx <- map(clusters, generate_matrix)
  names(mtx) <- clusters
  return(mtx)
}

# - Make pseudobulk metadata --------------------------------------------------
make_pseudobulk_metadata <- function(pseudobulk) {
  # make sure all colnames in pseudobulk matrix are the same
  if (!all(apply(sapply(pseudobulk, colnames), 1, 
                 function(x) length(unique(x)) == 1) == TRUE)) {
    stop("Cluster matrices do not have consistent column names")
  } else {
    df <- data.frame("Group" = colnames(pseudobulk[[1]])) %>%
      separate(Group, into = c("Treatment", "Batch"), sep = "_")
  }
}