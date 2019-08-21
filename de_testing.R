# - Choose PCs --------------------------------------------------------------
choose_pcs <- function(object, reps = 100) {
  
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
FindAllConservedMarkers <- function(object, groupby = NULL, clusters = NULL,
                                    verbose = FALSE) {
  if (is.null(clusters)) {
    clusters <- sort(unique(object@active.ident))
  }
  
  # Run FindConservedMarkers on each cluster
  markers <- map(clusters, 
                 ~ FindConservedMarkers(object, .x,
                                        grouping.var = groupby,
                                        only.pos = TRUE,
                                        verbose = verbose
                                        ))
  markers <- map(markers, as.data.frame)
  markers <- map(markers, ~ rownames_to_column(.x, "gene"))
  markers <- map(clusters,
                 ~ mutate(markers[[which(clusters == .x)]], "cluster" = .x))
  markers <- bind_rows(markers)
  
  # run logitp function from metap package
  logitp <- function(markers) {
    # select only p values from data.frame
    p_df <- select(markers, ends_with("p_val_adj"))
    
    # Calcaulte C value
    calc_C <- function(p_df) {
      k <- ncol(p_df)
      C <- sqrt((k*pi^2*(5*k + 2))/(3*(5*k+4)))
      return(C)
    }
    C <- calc_C(p_df)
    
    # convert each P value to log space
    convert_p <- function(p) {
      log(p / (1-p))
    }
    p_df <- apply(p_df, c(1,2), convert_p)
    
    # Calculate t value for each gene
    calc_t <- function(p_df) {
      apply(p_df, 1, function(x) -sum(x)/C)
    }
    p_val <- 2 * pt(calc_t(p_df), ncol(p_df)-1, lower.tail = FALSE)
    return(p_val)
  }
  p_val <- logitp(markers)
  
  # calculate pct cells expressing gene and fold change
  pct1 <- rowMeans(select(markers, ends_with("pct.1")))
  pct2 <- rowMeans(select(markers, ends_with("pct.2")))
  fc <- rowMeans(select(markers, ends_with("avg_logFC")))
  
  # select only relevant columns
  markers <- select(markers, cluster, gene) %>%
    mutate("pct.1" = pct1,
           "pct.2" = pct2,
           "avg_logFC" = fc,
           "p_val_adj" = p_val)
  markers <- filter(markers, p_val_adj < 0.05)
  markers <- markers %>%
    group_by(cluster) %>%
    arrange(desc(pct.1 - pct.2))
  return(markers)
}

# - Find unique genes -------------------------------------------------------
find_unique_genes <- function(object, genes = NULL, clusters = NULL) {
  if (is.null(clusters)) {
    clusters <- sort(unique(object@active.ident))
  }
  # grab cells from cluster and other clusters
  grab_cells <- function(cluster) {
    names(object@active.ident)[object@active.ident == cluster]
  }
  # set up data.frame
  df <- expand.grid(genes, clusters)
  df <- mutate(df, "gene" = Var1, "cluster" = Var2)
  # calculate pct expressing for each gene and cluster
  pct_expressing <- function(gene, cluster) {
    sum(object@assays$RNA@counts[gene, grab_cells(cluster)] > 0) /
      length(grab_cells(cluster))
  }
  df$pct <- apply(df, 1, function(x) pct_expressing(x[1], x[2]))
  
  # find largest difference
  by_gene <- df %>%
    group_by(gene) %>%
    summarize("first" = max(pct), "second" = nth(pct, 2)) %>%
    mutate("diff" = first - second) %>%
    arrange(desc(diff))
  
  by_cluster <- df %>%
    group_by(gene) %>%
    filter(pct == max(pct)) %>%
    select(-pct)
  
  return(left_join(by_gene, by_cluster, by = "gene"))
}
