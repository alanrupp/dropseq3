library(Seurat)
library(tidyverse)

# - Summarize data ------------------------------------------------------------
summarize_data <- function(object, genes, clusters = NULL) {
  if (is.null(clusters)) {
    clusters <- unique(object@active.ident)
  } else {
    clusters <- clusters[clusters %in% object@active.ident]
  }
  
  # keep valid genes
  genes <- genes[genes %in% rownames(object@assays$RNA@data)]
  
  # grab data
  if (length(genes) == 0) {
    return(NULL)
  } else if (length(genes) == 1) {
    df <- object@assays$RNA@data[genes, ] %>% as.data.frame() %>% t() %>% as.data.frame()
    rownames(df) <- genes
  } else {
    df <- object@assays$RNA@data[genes, ] %>% as.matrix() %>% as.data.frame()
  }
  
  # make tidy
  df <- df %>%
    rownames_to_column("gene") %>%
    gather(-gene, key = "barcode", value = "counts")
  
  # add cluster info
  df <- left_join(df, data.frame("barcode" = names(object@active.ident),
                                 "cluster" = object@active.ident), by = "barcode")
  
  # filter by selected clusters
  df <- filter(df, cluster %in% clusters)
  
  # summarize
  df <- df %>%
    group_by(gene, cluster) %>%
    summarize(avg = mean(counts), prop = sum(counts > 0)/n()) %>%
    ungroup()
  
  return(df)
}


# - Heatmap plot --------------------------------------------------------------
order_clusters <- function(object, cells = NULL, genes = NULL, scale = TRUE,
                           col_names = FALSE, row_names = FALSE, 
                           heatmap_legend = FALSE) {
  if (is.null(cells)) {
    cells <- colnames(object@assays$RNA@data)
  }
  if (is.null(genes)) {
    genes <- rownames(object@assays$RNA@data)
  } 
  
  cells <- cells[cells %in% colnames(object@assays$RNA@data)]
  genes <- genes[genes %in% rownames(object@assays$RNA@data)]
  
  if (scale == TRUE) {
    matrix <- object@assays$RNA@scale.data[genes, ]
  } else {
    matrix <- object@assays$RNA@data[genes, ]
  }
  
  # functions to grab cells and calculate mean expression
  mean_expression <- function(matrix, cells) {
    rowMeans(matrix[, cells])
  }
  cells_use <- function(object, ident) {
    names(object@active.ident)[object@active.ident == ident]
  }
  
  # calculate mean expression by cluster
  mean_values <-
    map(unique(object@active.ident),
        ~mean_expression(matrix, cells = cells_use(object, .x))) %>%
    as.data.frame() %>%
    set_names(unique(object@active.ident)) %>%
    t()
  
  # auto clip high and low to center at 0
  auto_clip <- function(mtx) {
    values <- as.matrix(mtx)
    if (abs(min(values)) < max(values)) {
      clip <- abs(min(values))
      mtx <- apply(mtx, c(1,2), function(x) ifelse(x > clip, clip, x))
    } else if (abs(min(values)) < max(values)) {
      clip <- max(values)
      mtx <- apply(mtx, c(1,2), function(x) ifelse(x < -clip, -clip, x))
    }
    return(mtx)
  }
  
  if (scale == TRUE) {
    mean_values <- auto_clip(mean_values)
  }
  
  dists <- dist(mean_values)
  clusts <- hclust(dists)
  
  return(clusts$labels[clusts$order])
}

# Heatmap plot ----------------------------------------------------------
heatmap_plot <- function(object, genes = NULL, cells = NULL, scale = TRUE,
                         cluster_order = NULL,
                         label_genes = FALSE, 
                         label_clusters = TRUE, 
                         heatmap_legend = FALSE,
                         cut_clusters = 1,
                         cut_genes = 1,
                         order_genes = FALSE) {
  # filter cells
  if (is.null(cells)) {
    cells <- colnames(object@assays$RNA@data)
    idents <- unique(object@active.ident)
  } else {
    cells <- cells[cells %in% colnames(object@assays$RNA@data)]
    idents <- unique(object@active.ident[names(object@active.ident) %in% cells])
  }
  # filter genes
  if (is.null(genes)) {
    genes <- rownames(object@assays$RNA@data)
  } 
  genes <- genes[genes %in% rownames(object@assays$RNA@data)]
  
  if (scale == TRUE) {
    matrix <- object@assays$RNA@scale.data[genes, cells]
    # make red --> white --> blue color palette
    endcolors <- c("red4", "white", "royalblue4")
    color_pal <- c(colorRampPalette(c(endcolors[1], endcolors[2]))(50), 
                   colorRampPalette(c(endcolors[2], endcolors[3]))(51)[-1])
  } else {
    matrix <- object@assays$RNA@data[genes, cells]
    color_pal <- viridis::viridis(10)
  }
  
  # functions to grab cells and calculate mean expression
  mean_expression <- function(cluster) {
    cell_subset <- names(object@active.ident)[object@active.ident == cluster]
    matrix_subset <- matrix[, cell_subset]
    return(Matrix::rowMeans(matrix_subset))
  }
  
  # calculate mean expression by cluster
  mean_values <-
    map(idents, mean_expression) %>%
    as.data.frame() %>%
    set_names(idents) %>%
    t()
  
  # auto clip high and low to center at 0
  auto_clip <- function(mtx) {
    values <- as.matrix(mtx)
    if (abs(min(values)) < max(values)) {
      clip <- abs(min(values))
      mtx <- apply(mtx, c(1,2), function(x) ifelse(x > clip, clip, x))
    } else if (abs(min(values)) < max(values)) {
      clip <- max(values)
      mtx <- apply(mtx, c(1,2), function(x) ifelse(x < -clip, -clip, x))
    }
    return(mtx)
  }
  
  if (scale == TRUE) {
    mean_values <- auto_clip(mean_values)
  }
  
  # plot
  if (order_genes == TRUE) {
    plt <-
      pheatmap::pheatmap(mean_values, 
                         show_colnames = label_genes,
                         show_rownames = label_clusters,
                         treeheight_col = 0,
                         legend = heatmap_legend,
                         color = color_pal,
                         cutree_rows = cut_clusters,
                         cluster_cols = FALSE)
  } else {
    plt <-
      pheatmap::pheatmap(mean_values, 
                         show_colnames = label_genes,
                         show_rownames = label_clusters,
                         treeheight_col = 0,
                         legend = heatmap_legend,
                         color = color_pal,
                         cutree_cols = cut_genes,
                         cutree_rows = cut_clusters)
  }
  return(plt)
}


# - Heatmap block -------------------------------------------------------------
heatmap_block <- function(object,
                          genes = NULL, 
                          cells = NULL,
                          clusters = NULL, 
                          n_cells = 1000,
                          scale = TRUE,
                          label_genes = TRUE, 
                          maxmin = NULL,
                          legend = TRUE,
                          integrated = FALSE) {
  
  # 
  if (ncol(object@assays$RNA@data) < n_cells) {
    n_cells <- ncol(object@assays$RNA@data)
  }
  
  # grab relevant clusters
  if (is.null(clusters)) {
    clusters <- sort(unique(object@active.ident))
    cells <- names(object@active.ident)[object@active.ident %in% clusters]
    cells <- sample(cells, n_cells)
    idents <- object@active.ident[cells]
    idents <- sort(idents)
    cells <- names(idents)
  } else {
    clusters <- clusters[clusters %in% object@active.ident]
    cells <- names(object@active.ident)[object@active.ident %in% clusters]
    cells <- sample(cells, n_cells)
    idents <- factor(object@active.ident[cells], levels = clusters)
    idents <- sort(idents)
    cells <- names(idents)
  }
  
  # filter genes
  if (is.null(genes)) {
    genes <- object@assays$integrated@var.features
  } 
  genes <- genes[genes %in% rownames(object@assays$RNA@data)]
  
  if (scale == TRUE) {
    matrix <- object@assays$RNA@scale.data[genes, cells]
  } else {
    matrix <- object@assays$RNA@data[genes, cells]
  }
  
  # auto clip high and low to center at 0
  auto_clip <- function(mtx, clip) {
    if (is.null(clip)) {
      values <- as.matrix(mtx)
      if (abs(min(values)) < max(values)) {
        clip <- abs(min(values))
        mtx <- apply(mtx, c(1,2), function(x) ifelse(x > clip, clip, x))
      } else if (abs(min(values)) > max(values)) {
        clip <- max(values)
        mtx <- apply(mtx, c(1,2), function(x) ifelse(x < -clip, -clip, x))
      }
    } else {
      mtx <- apply(mtx, c(1,2), function(x) ifelse(x < -clip, -clip, 
                                                   ifelse(x > clip, clip, x)))
    }
    return(mtx)
  }
  
  if (scale == TRUE) {
    matrix <- auto_clip(matrix, maxmin)
  }
  
  # make matrix into tidy df
  df <- as.data.frame(matrix) %>%
    rownames_to_column("gene") %>%
    gather(-gene, key = "cell", value = "z") %>%
    mutate(cell = factor(cell, levels = cells),
           gene = factor(gene, levels = genes))
  
  
  # generate axis labels to illustrate clusters
  label_df <- data.frame("cell" = cells) %>%
    left_join(., data.frame("cell" = names(object@active.ident),
                            "cluster" = object@active.ident), by = "cell") %>%
    mutate("pos" = seq(length(cells)))
  label_bar <- label_df %>%
    group_by(cluster) %>%
    summarize(xmin = min(pos), xmax = max(pos))
  label_text <- label_df %>%
    group_by(cluster) %>%
    summarize(pos = median(pos))
  label_colors <- sample(grDevices::colors(), length(unique(label_bar$cluster)))
  
  # basic plot
  plt <- 
    ggplot(df, aes(x = cell, y = fct_rev(gene))) +
    geom_tile(aes(fill = z), color = NA, show.legend = legend) +
    scale_fill_gradient2(low = "red4", mid = "white", high = "royalblue4",
                         name = "Expression") +
    xlab(NULL) + ylab(NULL) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          panel.border = element_blank())
  
  # add axis legend to plot
  #plt <- plt + 
  #  geom_rect(data = label_bar, 
  #             aes(xmin = xmin, xmax = xmax, ymin = -2, ymax = 0, fill = cluster),
  #             inherit.aes = FALSE, show.legend = FALSE) +
  #  geom_text(data = label_text,
  #            aes(x = pos, y = -10, color = cluster, label = cluster),
  #            inherit.aes = FALSE, show.legend = FALSE) +
  #  scale_color_manual(values = label_colors)
  
  return(plt)
}


# - Violin plot ---------------------------------------------------------------
violin_plot <- function(object, genes, tx = NULL, clusters = NULL, 
                        jitter = TRUE, stacked = FALSE, order_genes = FALSE,
                        ncol = NULL, flip = FALSE, void = FALSE, 
                        order_clusters = FALSE) {
  if (is.null(clusters)) {
    clusters <- sort(unique(object@active.ident))
  }
  
  # check that genes & clusters are in object
  genes <- genes[genes %in% rownames(object@assays$RNA@data)]
  clusters <- clusters[clusters %in% object@active.ident]
  
  # function to grab cells by treatment
  grab_cells <- function(clusters) {
    names(object@active.ident)[object@active.ident %in% clusters]
  }
  
  # grab data using input cells
  if (length(genes) == 1) {
    df <- object@assays$RNA@data[genes, grab_cells(clusters)] %>%
      as.matrix() %>% t() %>% as.data.frame()
    rownames(df) <- genes
  } else {
    df <- object@assays$RNA@data[genes, grab_cells(clusters)] %>%
      as.matrix() %>% as.data.frame()
  }
  
  # make tidy
  df <- df %>%
    rownames_to_column("gene") %>%
    gather(-gene, key = "cell", value = "Counts")
  
  if (order_genes == TRUE) {
    df <- mutate(df, gene = factor(gene, levels = genes))
  }
  
  # add cluster info
  cluster_df <- data.frame("cell" = names(object@active.ident),
                           "Cluster" = object@active.ident)
  df <- left_join(df, cluster_df, by = "cell")
  
  if (order_clusters == TRUE) {
    df <- mutate(df, Cluster = factor(Cluster, levels = clusters))
  }
  
  # plot ------------
  # add treatment info
  if (!is.null(tx)) {
    tx_df <- data.frame("Treatment" = object@meta.data$treatment) %>%
      rownames_to_column("cell") %>%
      filter(Treatment %in% tx)
    df <- inner_join(df, tx_df, by = "cell") %>%
      mutate(Treatment = factor(Treatment, levels = tx))
    plt <- ggplot(df, aes(x = Treatment, y = Counts, fill = Treatment)) +
      facet_wrap(~Cluster)
  } else {
    plt <- ggplot(df, aes(x = Cluster, y = Counts, fill = Cluster))
  }
  
  # generic plot attributes
  plt <- plt + 
    geom_violin(show.legend = FALSE, scale = "width", adjust = 1) +
    theme_bw() +
    scale_y_continuous(expand = c(0,0)) +
    theme(panel.grid = element_blank())
  
  # add jitter
  if (jitter == TRUE) {
    plt <- plt + geom_jitter(show.legend = FALSE, height = 0, alpha = 0.2,
                             stroke = 0)
  }
  
  # facet wrap for multiple genes or multiple genes + treatments
  if (length(genes) > 1) {
    if (!is.null(tx)) {
      plt <- plt + facet_wrap(~gene + Cluster, scales = "free_y")
    } else {
      plt <- plt + facet_wrap(~gene, scales = "free_y")
    }
  }
  
  if (void == TRUE) {
    plt <- plt + theme_void()
  }
  if (flip == TRUE) {
    plt <- plt + coord_flip()
  }
  
  return(plt)
}


# - Stacked violin plot -----------------------------------------------------
stacked_violin <- function(object, genes, cluster_order = NULL) {
  # grab cluster order
  if (is.null(cluster_order)) {
    clusters <- sort(unique(object@active.ident))
  } else {
    clusters <- factor(
      unique(object@active.ident[object@active.ident %in% cluster_order]),
      levels = cluster_order)
  }
  
  # make plots for all but last one
  data_plots <- map(genes, 
                  ~ violin_plot(object, genes = .x, 
                                jitter = FALSE, void = TRUE) #+
                    #theme(#axis.title = element_text(color = "black"), 
                    #      axis.title.x = element_blank()) + 
                    #ylab(.x)
                  )
  
  # make plots for all gene names for left-hand side
  gene_plots <- 
    map(genes, ~ ggplot(data.frame("gene" = .x),
           aes(x = 0, y = 0, label = gene)) +
          geom_text(hjust = "inward", color = "black") +
          theme_void()
    )
  
  # combine into one list
  plots <- list()
  for (i in seq(length(genes))) {
    plots <- c(plots, gene_plots[i])
    plots <- c(plots, data_plots[i])
  }
  
  # make empty square for bottom
  void_plot <- 
    ggplot(data.frame("a" = ""), aes(x = 0, y = 0, label = a)) +
    geom_text() + 
    theme_void()
  plots[[length(plots)+1]] <- void_plot
  
  # make x axis label plot
  cluster_plot <- 
    ggplot(
      data.frame("cluster" = sort(unique(object@active.ident)),
                 "position" = seq(length(unique(object@active.ident)))),
      aes(x = position, y = 1, label = cluster)) +
    geom_text(hjust = 0.5) +
    theme_void() +
    scale_x_continuous(expand = c(0, 0), 
                       limits = c(0.5, 
                                  length(unique(object@active.ident)) + 0.5)) +
    scale_y_continuous(expand = c(0, 0)) +
    xlab(NULL) + ylab(NULL)
  plots[[length(plots)+1]] <- cluster_plot
  
  # combine into one plot_grid
  cowplot::plot_grid(plotlist = plots, 
                     ncol = 2, 
                     rel_heights = c(rep(0.95/length(genes),
                                         length(genes)), 0.05),
                     rel_widths = c(0.1, 0.9),
                     align = "h")
}


# - UMAP plot ----------------------------------------------------------------
umap_plot <- function(object, genes = NULL, cells = NULL, clusters = NULL, 
                      legend = FALSE, cluster_label = FALSE,
                      ncol = NULL, xlim = NULL, ylim = NULL) {
  # pull UMAP data
  umap <- data.frame(
    UMAP1 = object@reductions$umap@cell.embeddings[, 1],
    UMAP2 = object@reductions$umap@cell.embeddings[, 2]
  ) %>% rownames_to_column("barcode")
  
  if (is.null(clusters)) {
    clusters <- sort(unique(object@active.ident))
  } else {
    cluster_bool <- object@active.ident %in% clusters
    clusters <- clusters[cluster_bool]
    umap <- umap[cluster_bool, ]
  }
  
  pull_data <- function(genes) {
    if (length(genes) == 1) {
      df <- object@assays$RNA@counts[genes, ] %>% as.data.frame() %>% 
        set_names(genes)
    } else {
      df <- object@assays$RNA@counts %>% as.matrix() %>% t() %>% as.data.frame()
    }
    return(df)
  }
  
  results <- pull_data(genes) %>%
    rownames_to_column("barcode") %>%
    left_join(., umap, by = "barcode") %>%
    select(-barcode) %>%
    gather(-starts_with("UMAP"), key = "gene", value = "value")
  
  # plot
  plt <-
    ggplot(results, aes(x = UMAP_1, y = UMAP_2)) +
    geom_point(aes(color = value), show.legend = legend, stroke = 0) +
    scale_color_gradient(low = "gray90", high = "navyblue") +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    facet_wrap(~gene, ncol = ncol) +
    xlab("UMAP1") + ylab("UMAP2")
  if (cluster_label == TRUE) {
    cluster_center <- 
      left_join(umap, data.frame("barcode" = names(object@active.ident),
                                 "cluster" = object@active.ident), by = "barcode") %>%
      select(-barcode)
    cluster_center <- cluster_center %>%
      group_by(cluster) %>%
      summarize(UMAP_1 = mean(UMAP_1), UMAP_2 = mean(UMAP_2))
    plt <- plt + geom_text(data = cluster_center, aes(label = cluster))
  }
  if (!is.null(xlim)) {
    plt <- plt + xlim(xlim)
  }
  if (!is.null(ylim)) {
    plt <- plt + ylim(ylim)
  }
  return(plt)
}


# - Doublet UMAP plot -------------------------------------------------------
doublet_umap_plot <- function(object, doublets) {
  data.frame(
    "UMAP1" = object@reductions$umap@cell.embeddings[,1],
    "UMAP2" = object@reductions$umap@cell.embeddings[,2],
    "cell" = names(object@active.ident)
  ) %>%
    mutate("Doublet" = ifelse(cell %in% doublets, TRUE, FALSE)) %>%
    arrange(Doublet) %>%
    mutate(cell = factor(cell, levels = cell)) %>%
    ggplot(aes(x = UMAP1, y = UMAP2, color = Doublet)) +
    geom_point(show.legend = FALSE, stroke = 0) +
    scale_color_manual(values = c("gray90", "firebrick3")) +
    theme_bw() +
    theme(panel.grid = element_blank())
}


# - Dot plot ------------------------------------------------------------------
dot_plot <- function(object, genes, clusters = NULL, 
                     gene_order = FALSE,
                     cluster_order = FALSE,
                     cluster_labels = FALSE) {
  
  # get summary data
  df <- summarize_data(object, genes, clusters)
  
  if (cluster_order == TRUE) {
    df <- df %>%
      mutate(cluster = factor(cluster, levels = clusters))
  }
  
  if (gene_order == TRUE) {
    df <- df %>%
      mutate(gene = factor(gene, levels = genes))
  }
  
  # plot
  plt <- 
    ggplot(df, aes(x = gene, y = fct_rev(cluster), size = avg,
                   color = prop)) +
    geom_point(show.legend = FALSE) +
    scale_x_discrete(position = "top") +
    scale_radius(limits = c(0.01, NA)) +
    scale_color_gradient(low = "gray", high = "dodgerblue3") +
    theme_void()
  
  if (cluster_labels == TRUE) {
    plt <- plt + theme(axis.text.x = element_text(color = "black", angle = 89, 
                                                  vjust = 0.5, hjust = 0.5),
                       axis.text.y = element_text())
  } else {
    plt <- plt + theme(axis.text.x = element_text(color = "black", angle = 89, vjust = 1,
                                                  hjust = ))
  }
  
  return(plt)
}

# - Flamemap plot ------------------------------------------------------------
flamemap <- function(object, genes, cells = NULL, n_bars = 100,
                     order_genes = FALSE, cluster_labels = TRUE,
                     icicle = FALSE) {
  genes <- genes[genes %in% rownames(object@assays$RNA@data)]
  if (is.null(cells)) {
    cells <- colnames(object@assays$RNA@data)
  } else {
    cells <- cells[cells %in% colnames(object@assays$RNA@data)]
  }
  
  if (length(genes) == 0) {
    print("Invalid gene input")
    exit()
  }
  
  # grab cluster information and arrange by cluster name
  clusters <- data.frame(cell = names(object@active.ident), 
                         cluster = object@active.ident)
  clusters <- arrange(clusters, cluster, cell)
  clusters <- filter(clusters, cell %in% cells)
  
  # grab cluster info for plotting ticks on x axis
  scale_factor <- nrow(clusters)/n_bars
  cluster_stats <- 
    clusters %>% 
    mutate("pos" = seq(n())) %>% 
    group_by(cluster) %>% 
    summarize(max = max(pos), mid = mean(pos)) %>%
    mutate_if(is.numeric, ~ .x / scale_factor)
  
  # grab relevant data
  mtx <- object@assays$RNA@data[genes, filter(clusters, cell %in% cells)$cell]
  
  # reshape
  df <- 
    mtx %>%
    as.matrix() %>%
    as.data.frame() %>%
    rownames_to_column("gene") %>%
    gather(-gene, key = "cell", value = "counts") %>%
    mutate(cell = factor(cell, levels = filter(clusters, cell %in% cells)$cell)) %>%
    group_by(gene) %>%
    mutate("bar" = ntile(cell, n_bars)) %>%
    group_by(gene, bar) %>%
    arrange(desc(counts)) %>%
    ungroup() %>%
    group_by(gene, bar) %>%
    mutate("height" = as.numeric(counts > 1) / n())
  
  # order genes by user-defined order
  if (order_genes == TRUE) {
    df <- df %>%
      mutate(gene = factor(gene, levels = genes))
  }
  
  if (icicle == TRUE) {
    df <- mutate(df, counts = -counts, height = -height)
  }
  
  # plot
  plt <-
    ggplot(df, aes(x = bar, y = height, fill = counts)) +
    geom_col(show.legend = FALSE) +
    geom_hline(aes(yintercept = 0)) +
    scale_x_continuous(expand = c(0, 0), breaks = c(1, cluster_stats$max+0.5)) +
    labs(y = element_blank(), x = element_blank()) +
    facet_wrap(~gene) +
    theme_bw() +
    theme(panel.grid = element_blank(), 
          axis.text.x = element_blank(),
          axis.line.x = element_blank())
  
  if (icicle == TRUE) {
    plt <- plt + scale_fill_gradient(low = "royalblue", high = "aquamarine") +
      geom_text(data = cluster_stats,
                aes(x = mid, y = 0.015, label = cluster), inherit.aes = FALSE) +
      scale_y_continuous(expand = c(0, 0), limits = c(-1, 0.03),
                         labels = seq(1, 0, by = -0.25))
  } else {
    plt <- plt + scale_fill_gradient(low = "yellow", high = "red") +
      geom_text(data = cluster_stats,
                aes(x = mid, y = -0.015, label = cluster), inherit.aes = FALSE)  +
      scale_y_continuous(expand = c(0, 0), limits = c(-0.03, 1))
  }
  
  return(plt)
}
