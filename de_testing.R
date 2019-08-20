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
FindAllConservedMarkers <- function(object, clusters = NULL, groupby = NULL) {
  if (is.null(clusters)) {
    clusters <- sort(unique(object@active.ident))
  }
  
  markers <- map(clusters, ~ FindConservedMarkers(object, .x, 
                                                  grouping.var = groupby))
  markers <- map(markers, as.data.frame)
  markers <- map(markers, ~ rownames_to_column(.x, "gene"))
  markers <- map(clusters,
                 ~ mutate(markers[[which(clusters == .x)]], "cluster" = .x))
  markers <- bind_rows(markers)
  return(markers)
}
