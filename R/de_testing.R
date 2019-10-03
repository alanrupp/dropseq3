library(Seurat)
library(tidyverse)

# - Grab cell names by cluster membership or dismembership -------------------
cluster_cells <- function(object, cluster) {
  names(object@active.ident)[object@active.ident == cluster]
}
other_cells <- function(object, cluster) {
  names(object@active.ident)[object@active.ident != cluster]
}

# - Cluster means -------------------------------------------------------------
cluster_means <- function(object, genes = NULL, assay = "RNA", 
                          data_slot = "scale.data") {
  clusters <- sort(unique(object@active.ident))
  mtx <- slot(slot(object, "assays")[[assay]], data_slot)
  if (is.null(genes)) {
    genes <- rownames(mtx)
  } else {
    genes <- genes[genes %in% rownames(object@assays$RNA@data)]
  }
  df <-
    map(clusters,
        ~ Matrix::rowMeans(mtx[genes, cluster_cells(object, .x)])
    ) %>% 
    set_names(clusters) %>%
    bind_cols() %>%
    as.data.frame()
  rownames(df) <- genes
  return(df)
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
                              assay = "integrated") {
  dims <- 1:object@reductions$pca@misc$sig_pcs
  
  while (TRUE) {
    object <- FindClusters(object, resolution = resolutions, verbose = FALSE)
    
    # calc silhouettes for all resolutions
    calc_silhouettes <- function(resolution, assay = "integrated") {
      col <- paste0(assay, "_snn_res.", resolution)
      sil <- cluster::silhouette(
        as.numeric(object@meta.data[, col]), 
        cluster::daisy(object@reductions$pca@cell.embeddings[, dims])
      )
      return(sil[,3])
    }
    silhouettes <- map(resolutions, calc_silhouettes)
    
    # choose resolution with highest mean silhouette width
    best_width <- silhouettes %>% 
      set_names(resolutions) %>%
      map_dbl(., mean) %>%
      sort(decreasing = TRUE) %>%
      .[1] %>%
      names()
    
    if (best_width != max(resolutions)) {
      break()
    }
    resolutions <- seq(max(resolutions) + 0.2, max(resolutions) + 1, by = 0.2)
  }
  
  print(paste("Using resolution", best_width, "for clustering."))
  
  # assign clusters to active.ident slot
  ident_name <- paste0(assay, "_snn_res.", best_width)
  object@active.ident <- object@meta.data[, ident_name]
  names(object@active.ident) <- rownames(object@meta.data)
  return(object)
}

# - Cluster -------------------------------------------------------------------
cluster <- function(object, assay = "integrated") {
  if (assay == "integrated") {
    DefaultAssay(object) <- "integrated"
  } else {
    DefaultAssay(object) <- "RNA"
  }
  pcs <- object@reductions$pca@misc$sig_pcs
  # Find neighbors and cluster and different resolutions
  object <- FindNeighbors(object, dims = 1:pcs, verbose = FALSE)
  object <- choose_resolution(object, resolutions = seq(0.4, 2, by = 0.2))
  return(object)
}

# - Optimize UMAP ------------------------------------------------------------
optimize_umap <- function(object, method = "dunn") {
  pcs <- object@reductions$pca@misc$sig_pcs
  mtx <- object@reductions$pca@cell.embeddings[,1:pcs]
  neighbors <- seq(5, 50, by = 5)
  dists <- c(0.001, 0.005, 0.01, 0.05, 0.1, 0.5)
  search_grid <- expand.grid(neighbors, dists)
  # get UMAP coordinates
  get_coordinates <- function(neighbors, distance) {
    as.data.frame(uwot::umap(mtx, n_neighbors = neighbors, min_dist = distance))
  }
  results <- map(1:nrow(search_grid), 
                 ~ get_coordinates(search_grid[.x, ]$Var1, search_grid[.x, ]$Var2))
  
  if (method == "dunn") {
    # calculate Dunn index for compactness of clusters
    dunn_index <- function(result) {
      distances <- dist(result) %>% as.matrix()
      clusters <- unique(object@active.ident)
      get_min <- function(cluster) {
        min(distances[object@active.ident == cluster, object@active.ident != cluster])
      }
      get_max <- function(cluster) {
        max(distances[object@active.ident == cluster, object@active.ident == cluster])
      }
      result <- map_dbl(clusters, ~ get_min(.x) / get_max(.x))
      return(mean(result))
      }
    scores <- map_dbl(results, dunn_index)
    } else if (method == "silhouette") {
    # calculate silhouette width for each coordinate state
    calc_silhouettes <- function(result) {
      sil <- cluster::silhouette(as.numeric(object@active.ident),
                                 cluster::daisy(result)
      )
      return(mean(sil[, 3]))
    }
    scores <- map_dbl(results, calc_silhouettes)
  }
  search_grid$mean_scores <- scores
  search_grid <- rename(search_grid,
                        "n_neighbors" = Var1,
                        "min_dist" = Var2)
  search_grid <- arrange(search_grid, desc(mean_scores))
  return(search_grid)
}

# - Run UMAP ------------------------------------------------------------------
run_umap <- function(object) {
  parameters <- optimize_umap(object)
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

# - Cluster identity all cells -----------------------------------------------
find_classes <- function(object, markers) {
  classes <- read_csv("~/Programs/dropseq3/data/celltype_markers.csv")
  classes <- classes %>% mutate("score" = str_count(paper, "\\,") + 1)
  
  # grab unique clusters & classes
  clusters <- sort(unique(object@active.ident))
  unique_classes <- sort(unique(classes$cells))
  
  # find most likely class for each cluster
  find_class <- function(clstr) {
    # get list of DE genes from that cluster in descending order of enrichment
    de_genes <- 
      markers %>%
      arrange(desc(pct.1 - pct.2)) %>%
      filter(cluster == clstr) %>% 
      .$gene
    
    # find cumulative score for each gene
    class_score <- function(class) {
      this_class <- filter(classes, cells == class)
      other_classes <- filter(classes, cells != class)
      each_score <- 
        map_dbl(de_genes, 
            ~ ifelse(.x %in% this_class$gene,
                     filter(this_class, gene == .x)$score, 
                     ifelse(.x %in% other_classes$gene,
                            1 - filter(other_classes, gene == .x)$score, 0
                            )
                     )
      ) %>% unlist()
      cumulative_score <- cumsum(each_score)
      if (max(cumulative_score) == 0) {
        peak <- NA
      } else {
        peak <- min(which(cumulative_score == max(cumulative_score)))
      }
      return(peak)
    }
    
    result <- map(unique_classes, class_score) %>% unlist()
    best_hit <- which(result == max(result, na.rm = TRUE))
    if (length(best_hit) == 0) {
      return(NA)
    } else if (length(best_hit) == 1) {
      return(unique_classes[best_hit])
    } else 
      return(unique_classes[best_hit[1]])
  }
  
  # run for all clusters
  results <- map_chr(clusters, find_class)
  results <- data.frame("cluster" = clusters, "class" = results)
  
  # Find highest DE gene from hit class for each cluster
  classic_markers <- function(classes_df) {
    markers_df <- arrange(markers, desc(pct.1 - pct.2))
    clusters <- unique(classes_df$cluster)
    
    classic_marker <- function(clstr) {
      class <- filter(classes_df, cluster == clstr)$class
      if (is.na(class)) {
        hit <- NA
      } else {
        marker_genes <- filter(classes, cells == class)$gene
        hit <- filter(markers, cluster == clstr & gene %in% marker_genes)$gene[1]
      }
      return(hit)
    }
    return(map(clusters, classic_marker) %>% unlist())
  }
  results$marker <- classic_markers(results)
  
  return(results)
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
  too_few <- table(object@active.ident, object$mouse) %>%
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
  mtx <- object@assays$RNA@counts
  if (is.null(clusters)) {
    clusters <- sort(unique(object@active.ident))
  }
  if (is.null(genes)) {
    genes <- rownames(mtx)
  }
  
  cluster_expression <- function(mtx) {
    df <-
      map(clusters,
          ~ Matrix::rowMeans(mtx[genes, cluster_cells(object, .x)] > 0)
      ) %>% 
      set_names(clusters)
    return(df)
  }
  df <- cluster_expression(mtx)
  df <- bind_cols(df)
  
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
    # find 2 neighbors for every cluster based correlation of PCs
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
edgeR_test <- function(mtx, metadata, treatment) {
  library(edgeR)
  y <- DGEList(counts = mtx, group = metadata[, treatment])
  y <- calcNormFactors(y)
  design <- model.matrix(~ 0 + metadata[, treatment])
  colnames(design) <- levels(metadata[, treatment])
  y <- estimateDisp(y, design)
  
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
grab_type <- function(object, classes, type = "Neuron") {
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
    result[i-1, 1] <- i
    result[i-1, 2] <- j
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
