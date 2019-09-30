library(Seurat)
library(tidyverse)

# - Choose PCs --------------------------------------------------------------
choose_pcs <- function(object, reps = 100) {
  
}

# - Cluster -------------------------------------------------------------------
cluster <- function(object, pcs = 50) {
  object <- RunPCA(object, verbose = FALSE)
  object <- JackStraw(object, dims = pcs)
  object <- ScoreJackStraw(object, dims = 1:pcs)
  insig <- which(
    object@reductions$pca@jackstraw@overall.p.values[, "Score"] >= 0.05
  )
  while (length(insig) == 0) {
    pcs <- pcs + 50
    object <- RunPCA(object, verbose = FALSE, npcs = pcs)
    object <- JackStraw(object, dims = pcs)
    object <- ScoreJackStraw(object, dims = 1:pcs)
    insig <- which(
      object@reductions$pca@jackstraw@overall.p.values[,"Score"] >= 0.05
    )
  }
  pcs <- insig[1] - 1
  print(paste("Using", pcs, "PCs"))
  # Find neighbors and cluster
  object <- FindNeighbors(object, dims = 1:pcs)
  object <- choose_resolution(object)
  object <- RunUMAP(object, dims = 1:pcs, verbose = FALSE)
  return(object)
}

# - Choose resolution --------------------------------------------------------
choose_resolution <- function(object, 
                              resolutions = seq(0.4, 2, by = 0.2),
                              assay = "integrated") {
  object <- FindClusters(object, resolution = resolutions)
  
  # calc silhouettes for all resolutions
  calc_silhouettes <- function(resolution, assay = "integrated", pcs) {
    col <- paste0(assay, "_snn_res.", resolution)
    sil <- cluster::silhouette(
      as.numeric(object@meta.data[, col]), 
      cluster::daisy(object@reductions$pca@cell.embeddings[, 1:pcs])
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
  print(paste("Using resolution", best_width, "for clustering."))
  
  # assign clusters to active.ident slot
  ident_name <- paste0(assay, "_snn_res.", best_width)
  object@active.ident <- object@meta.data[, ident_name]
  names(object@active.ident) <- rownames(object@meta.data)
}

# find highly expressed genes -----------------------------------------------
find_expressed <- function(object, min_cells = 10, min_counts = 1) {
  calc <- rowSums(object@raw.data >= min_counts) >= min_cells
  calc <- calc[calc == TRUE]
  return(names(calc))
}

# - Cluster identity all cells -----------------------------------------------
find_classes <- function(object, markers_df) {
  classes <- read_csv("~/Programs/dropseq3/data/celltype_markers.csv")
  classes <- classes %>% mutate("score" = str_count(paper, "\\,") + 1)
  
  # grab unique clusters & classes
  clusters <- sort(unique(object@active.ident))
  unique_classes <- sort(unique(classes$cells))
  
  # find most likely class for each cluster
  find_class <- function(clstr) {
    de_genes <- arrange(markers_df, desc(pct.1 - pct.2)) %>%
      filter(cluster == clstr) %>% .$gene
    
    # find cumulative score for each gene
    class_score <- function(class) {
      this_class <- filter(classes, cells == class)
      other_classes <- filter(classes, cells != class)
      each_score <- map(de_genes, 
                        ~ ifelse(.x %in% this_class$gene, 
                                 filter(this_class, gene == .x)$score, 
                                 ifelse(.x %in% other_classes$gene,
                                        1-filter(other_classes, gene == .x)$score, 0)
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
  results <- map(clusters, find_class) %>% unlist()
  results <- data.frame("cluster" = clusters, "class" = results)
  
  # Find highest DE gene from hit class for each cluster
  classic_markers <- function(classes_df) {
    markers_df <- arrange(markers_df, desc(pct.1 - pct.2))
    clusters <- unique(classes_df$cluster)
    
    classic_marker <- function(clstr) {
      class <- filter(classes_df, cluster == clstr)$class
      if (is.na(class)) {
        hit <- NA
      } else {
        marker_genes <- filter(classes, cells == class)$gene
        hit <- filter(markers_df, cluster == clstr & gene %in% marker_genes)$gene[1]
      }
      return(hit)
    }
    return(map(clusters, classic_marker) %>% unlist())
  }
  results$marker <- classic_markers(results)
  
  return(results)
}

# - Find all conserved markers -----------------------------------------------
FindAllConservedMarkers <- function(object, ident2 = NULL, 
                                    groupby = NULL, clusters = NULL,
                                    verbose = TRUE) {
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
  markers <- select(markers, cluster, gene) %>%
    mutate("avg_logFC" = fc,
           "pct.1" = pct1,
           "pct.2" = pct2,
           "p_val_adj" = p_vals)
  
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
# merge markerless clusters with nearest neighbor by both correlation and UMAP
merge_markerless <- function(object, marker_summary) {
  markerless <- filter(marker_summary, total == 0)$cluster
  if (length(markerless) == 1) {
    cat(paste("Cluster", markerless, 
              "has no significantly enriched genes. Merging with neighbors or dropping. ")
    )
  } else if (length(markerless) > 1) {
    cat(paste("Clusters", paste(markerless, collpase = ", "), 
              "have no significantly enriched genes. Merging with neighbors or dropping. ")
    )
  } else {
    return(object)
  }
  genes <- unique(markers$gene)
  genes <- genes[genes %in% rownames(srt@assays$RNA@scale.data)]

  # find correlation between clusters and markerless
  grab_cells <- function(cluster) {
    names(object@active.ident)[object@active.ident == cluster]
  }
  correlation <- function(cluster) {
    corr <- 
      map(unique(object@active.ident),
          ~ cor(Matrix::rowMeans(object@assays$RNA@scale.data[genes, grab_cells(.x)]),
                Matrix::rowMeans(object@assays$RNA@scale.data[genes, grab_cells(cluster)]))
    ) %>% unlist()
    names(corr) <- unique(object@active.ident)
    # remove self-correlation
    corr <- corr[corr != 1]
    # only keep correlations greater than the median
    neighbor <- names(corr)[corr == max(corr)]
    neighbor <- as.integer(neighbor)
    return(neighbor)
  }
  
  # find umap neighbor
  umap_neighbor <- function(clstr) {
    cluster_centers <- data.frame(
      "cluster" = object@active.ident,
      "UMAP1" = object@reductions$umap@cell.embeddings[,1],
      "UMAP2" = object@reductions$umap@cell.embeddings[,2]
    ) %>%
      group_by(cluster) %>%
      summarize(UMAP1 = median(UMAP1), UMAP2 = median(UMAP2)) %>%
      as.data.frame() %>%
      column_to_rownames("cluster")
    dists <- 
      dist(cluster_centers) %>%
      as.matrix() %>%
      as.data.frame() %>%
      rownames_to_column("cluster") %>%
      gather(-cluster, key = "other", value = "dist")
    dists <- filter(dists, dist != 0)
    min_dist <- dists %>%
      group_by(cluster) %>%
      summarize(closest = min(dist))
    neighbor <- filter(dists, cluster == clstr & 
                         dist == filter(min_dist, cluster == clstr)$closest)
    if (neighbor$dist < median(min_dist$closest)) {
      return(neighbor$other)
      print(paste("umap", class(neighbor$other)))
    } else {
      return(NULL)
    }
  }
  
  # find neighbors for all markerless clusters
  corr <- vector()
  umap <- vector()
  for (i in seq(length(markerless))) {
    corr[i] <- correlation(markerless[i])
    umap[i] <- umap_neighbor(markerless[i])
  }
  result <- data.frame("corr" = corr, "umap" = umap)
  rownames(result) <- markerless
  
  # merge or drop cells based on neighbors
  merge_cells <- function(cluster, other) {
    cat(paste("Cluster", cluster, "and cluster", other,
              "are most similar by expression correlation and UMAP distance.",
              "Merging cluster", cluster, "into cluster", other))
    object@active.ident[object@active.ident == cluster] <- other
    object@active.ident <- factor(object@active.ident, 
                                  levels = sort(unique(object@active.ident)))
    return(object)
  }
  drop_cells <- function(cluster) {
    cat(paste("Dropping cluster", cluster, "because it has no clear larger",
              "cluster to merge into."))
    object <- SubsetData(object, ident.remove = cluster)
    object@active.ident <- factor(object@active.ident, 
                                  levels = sort(unique(object@active.ident)))
    return(object)
  }
  for (i in seq(nrow(result))) {
    # if both umap and corr are the same, merge clusters
    if (all(result[i, ] == result[i,1])) {
      object <- merge_cells(rownames(result[i,]), result[i,1])
    } else {
      object <- drop_cells(rownames(result[i,]))
    }
  }
  return(object)
}

# - Find unique genes -------------------------------------------------------
find_unique_genes <- function(object, genes = NULL, clusters = NULL,
                              top_n = 1) {
  if (is.null(clusters)) {
    clusters <- sort(unique(object@active.ident))
  }
  # grab cells from cluster and other clusters
  grab_cells <- function(cluster) {
    names(object@active.ident)[object@active.ident == cluster]
  }
  # set up data.frame
  df <- expand.grid(genes, clusters)
  df <- rename(df, "gene" = Var1, "cluster" = Var2)
  # calculate pct expressing for each gene and cluster
  pct_expressing <- function(gene, cluster) {
    sum(object@assays$RNA@counts[gene, grab_cells(cluster)] > 0) /
      length(grab_cells(cluster))
  }
  df$pct <- apply(df, 1, function(x) pct_expressing(x[1], x[2]))
  
  # find largest difference
  by_gene <- df %>%
    ungroup() %>%
    arrange(desc(pct)) %>%
    group_by(gene) %>%
    summarize("first" = max(pct), "second" = nth(pct, 2)) %>%
    mutate("diff" = first - second) %>%
    arrange(desc(diff))
  by_cluster <- df %>%
    group_by(gene) %>%
    filter(pct == max(pct)) %>%
    select(-pct)
  
  # return top n gene(s) for each cluster
  result <- inner_join(by_gene, by_cluster, by = "gene")
  result <- result %>%
    group_by(cluster) %>% 
    slice(top_n) %>%
    ungroup() %>%
    arrange(cluster)
  
  return(result)
}

# - Eigengene ---------------------------------------------------------------
eigengene <- function(object, genes) {
  genes <- genes[genes %in% rownames(object@assays$RNA@scale.data)]
  pca <- prcomp(t(object@assays$RNA@scale.data[genes, ]))
  return(pca$x[,1])
}


# - Merging clusters --------------------------------------------------------
# merge clusters according to the Tasic et al Nature 2018 criteria 
merge_clusters <- function(object) {
  # find 2 neighbors for each cluster based on Euclidean distance in UMAP
  cluster_centers <- 
    data.frame("UMAP1" = object@reductions$umap@cell.embeddings[,1],
               "UMAP2" = object@reductions$umap@cell.embeddings[,2],
               "cluster" = object@active.ident) %>%
    group_by(cluster) %>%
    summarize(UMAP1 = median(UMAP1), UMAP2 = median(UMAP2)) %>%
    as.data.frame() %>%
    column_to_rownames("cluster")
  dists <- 
    dist(cluster_centers) %>%
    as.matrix() %>%
    as.data.frame() %>%
    rownames_to_column("cluster") %>%
    gather(-cluster, key = "other", value = "dist")
  dists <- filter(dists, dist != 0)
  find_neighbors <- function(clstr) {
    dists <- filter(dists, cluster == clstr)
    dists <- arrange(dists, dist)
    dists <- slice(dists, 1:2)
    return(dists$other)
  }
  neighbors <- map(unique(dists$cluster), find_neighbors)
  neighbors <- unlist(neighbors)
  neighbors <- data.frame("cluster" = rep(unique(dists$cluster), each = 2),
                          "neighbor" = neighbors)
  
  # calculate "deScore"
  deScore <- function(ident1, ident2) {
    result <- FindMarkers(object, ident.1 = ident1, ident.2 = ident2)
    result <- mutate(result, p_val_adj = -log10(p_val_adj))
    result <- mutate(result, p_val_adj = ifelse(p_val_adj > 20, 20, p_val_adj))
    result <- sum(result$p_val_adj)
    return(result)
  }
  neighbors$deScore <- apply(neighbors, 1, function(x) deScore(x[1], x[2]))
  
  # merge clusters with deScore < 150
  to_merge <- filter(neighbors, deScore < 150)
  if (nrow(to_merge) == 0) {
    cat("No clusters to merge")
  } else {
    for (i in 1:nrow(to_merge)) {
      cat(paste("Merging clusters", to_merge[i, "cluster"], "and",
                to_merge[i, "neighbor"]))
      object@active.ident[object@active.ident == to_merge[i, "neighbor"]] <-
        to_merge[i, "cluster"]
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
