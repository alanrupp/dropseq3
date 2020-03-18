library(Seurat)
library(tidyverse)
library(ggrepel)
library(ggdendro)

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
    df <- object@assays$RNA@data[genes, ] %>% as.data.frame() %>% t() %>% 
      as.data.frame()
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
  
  if (!is.null(clusters)) {  df <- filter(df, cluster %in% clusters) }
  # summarize
  df <- df %>%
    group_by(gene, cluster) %>%
    summarize(avg = mean(counts), prop = sum(counts > 0)/n()) %>%
    ungroup()
  
  return(df)
}


# - GG color hue --------------------------------------------------------------
# from John Colby on Stack Overflow
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

# - Clip matrix ---------------------------------------------------------------
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

# Heatmap plot ----------------------------------------------------------
heatmap_plot <- function(object, genes = NULL, cells = NULL, scale = TRUE,
                         label_genes = FALSE,
                         heatmap_legend = FALSE, title = NA,
                         max_expr = NULL,
                         cluster_genes = TRUE, cluster_clusters = TRUE,
                         draw_tree = TRUE, tree_scaling = 1,
                         flipped = FALSE) {
  
  # calculate mean expression by cluster
  mean_values <- cluster_means(object, genes)
  
  # cluster genes to group expression patterns
  if (cluster_genes) {
    gene_clusters <- hclust(dist(mean_values))
    gene_order <- gene_clusters$label[gene_clusters$order]
  } else {
    gene_order <- genes
  }
  
  # cluster clusters to group expression patterns
  if (cluster_clusters) {
    cluster_clusters <- hclust(dist(t(mean_values)))
    cluster_order <- cluster_clusters$label[cluster_clusters$order]
  } else {
    cluster_order <- levels(object@active.ident)
  }
  
  # tidy
  mean_values <- as.data.frame(mean_values) %>%
    rownames_to_column("gene") %>%
    gather(-gene, key = "cluster", value = "expr") %>%
    mutate(gene = factor(gene, levels = gene_order)) %>%
    mutate(cluster = factor(cluster, levels = cluster_order))
  
  if (!is.null(max_expr)) {
    mean_values <- mean_values %>%
      mutate(expr = ifelse(expr > max_expr, max_expr,
                           ifelse(expr < -max_expr, -max_expr, expr)))
  }
  
  # plot
  p <- ggplot(mean_values, aes(x = cluster, y = fct_rev(gene))) +
    geom_tile(aes(fill = expr), color = NA) +
    scale_fill_gradient2(low = "#d0587e", mid = "white", high = "#009392",
                         name = "z-score") +
    theme_void() +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    xlab(NULL) + ylab(NULL) +
    theme(axis.text.x = element_text(color = "black"))
  if (label_genes) {
    p <- p + theme(axis.text = element_text(color = "black"))
  }
  if (draw_tree) {
    dg <- dendro_data(cluster_clusters)
    p <- p + 
      geom_segment(data = dg$segments, 
                   aes(x = x, xend = xend,
                       y = (y * tree_scaling) + length(genes) + 0.5 , 
                       yend = (yend * tree_scaling) + length(genes) + 0.5))
  }
  if (flipped) {
    p <- p + coord_flip() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
  }
  return(p)
}

# - Clip matrix at upper limit ------------------------------------------------
clip_matrix <- function(mtx, limit) {
  cols <- colnames(mtx); rows <- rownames(mtx)
  clip_fun <- function(x) {
    ifelse(x > limit, limit, ifelse(x < -limit, -limit, x))
  }
  mtx <- apply(mtx, c(1, 2), clip_fun)
  colnames(mtx) <- cols; rownames(mtx) <- rows
  return(mtx)
}

# - Heatmap block -------------------------------------------------------------
heatmap_block <- function(object, genes = NULL, cells = NULL,
                          clusters = NULL, n_cells = 1000,
                          scale = TRUE, label_genes = TRUE, 
                          maxmin = NULL,
                          integrated = FALSE,
                          legend = FALSE) {
  
  # 
  if (ncol(object@assays$RNA@data) < n_cells) {
    n_cells <- ncol(object@assays$RNA@data)
  }
  # grab relevant clusters
  if (is.null(clusters)) {
    clusters <- levels(object@active.ident)
    cells <- names(object@active.ident)[object@active.ident %in% clusters]
    cells <- sample(cells, n_cells)
    idents <- sort(object@active.ident[cells])
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
  } else {
    missed <- !genes %in% rownames(object@assays$RNA@counts)
    if (sum(missed) > 0) {
      warning(paste(paste(genes[missed], collapse = ", "), 
                    "not present in dataset."))
      genes <- genes[!missed]
    }
  }
  
  if (scale) {
    matrix <- object@assays$RNA@scale.data[genes, cells]
  } else {
    matrix <- object@assays$RNA@data[genes, cells]
  }
  
  if (!is.null(maxmin) & scale) {
    matrix <- clip_matrix(matrix, maxmin)
  }
  
  # make matrix into tidy df
  df <- as.data.frame(matrix) %>%
    rownames_to_column("gene") %>%
    gather(-gene, key = "cell", value = "z") %>%
    mutate(cell = factor(cell, levels = cells),
           gene = factor(gene, levels = genes))
  
  # generate axis labels to illustrate clusters
  label_text <- data.frame("cell" = cells) %>%
    left_join(., data.frame("cell" = names(object@active.ident),
                            "cluster" = object@active.ident), 
              by = "cell") %>%
    mutate("pos" = seq(length(cells)))
  label_text_avg <- label_text %>%
    group_by(cluster) %>%
    summarize(pos = mean(pos))
  
  # - plots -----------
  plt_list <- list()
  
  # basic plot
  heatmap_plt <- 
    ggplot(df, aes(x = cell, y = fct_rev(gene))) +
    geom_tile(aes(fill = z), color = NA, show.legend = legend) +
    scale_fill_gradient2(low = "#d0587e", mid = "white", high = "#009392",
                         name = "Expression") +
    theme_void() +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    xlab(NULL) + ylab(NULL)
  plt_list[[1]] <- heatmap_plt
  
  # plot axis labels
  label_plt <- 
    ggplot(label_text_avg, aes(x = pos, y = 0)) +
    geom_text_repel(aes(label = cluster), ylim = c(-Inf, 0), direction = "y",
                    hjust = 0.5) +
    theme_void() +
    scale_x_continuous(expand = c(0, 0), limits = c(1, max(label_text$pos))) +
    scale_y_continuous(expand = c(0, 0), limits = c(-100, 2)) +
    xlab(NULL) + ylab(NULL)
  for (i in unique(label_text$cluster)) {
    label_plt <- label_plt +
      annotate("segment", 
               x = min(filter(label_text, cluster == i)$pos) + 3, 
               xend = max(filter(label_text, cluster == i)$pos) -3, 
               y = 1, 
               yend = 1, color = "black")
  }
  plt_list[[2]] <- label_plt
  
  return(cowplot::plot_grid(plotlist = plt_list, ncol = 1,
                            rel_heights = c(0.8, 0.2))
  )
}


# - Violin plot ---------------------------------------------------------------
violin_plot <- function(object, genes, x = "cluster", group = NULL,
                        fill = NULL, facet = NULL, data = "data", 
                        clusters = NULL, n_col = 1, jitter = FALSE,
                        colors = NULL, void = FALSE, flip = FALSE,
                        order_genes = TRUE, stacked = FALSE) {
  genes <- genes[genes %in% rownames(slot(object@assays$RNA, data))]
  # get data
  object$cluster <- object@active.ident
  meta <- data.frame("na" = matrix(NA, nrow = ncol(slot(object@assays$RNA, data))))
  meta[, x] <- object@meta.data[, x]
  meta <- select(meta, -na)
  if (!is.null(group)) { meta[, group] <- object@meta.data[, group] }
  if (!is.null(fill)) { meta[, fill] <- object@meta.data[, fill] }
  if (!is.null(facet)) { meta[, facet] <- object@meta.data[, facet] }
  if (is.null(colors) & is.null(fill)) { 
    colors <- gg_color_hue(length(unique(meta[, x]))) 
  } else if (is.null(colors)) {
    colors <- gg_color_hue(length(unique(meta[, fill]))) 
  }
  if (length(genes) > 1) {
    df <- data.frame(t(as.matrix(slot(object@assays$RNA, data)[genes, ])))
    df <- bind_cols(df, meta)
    df <- gather(df, 1:length(genes), key = "gene", value = "y")
  } else {
    df <- data.frame("y" = slot(object@assays$RNA, data)[genes, ], "gene" = genes)
    df <- bind_cols(meta, df)
  }
  if (order_genes) { df <- mutate(df, gene = factor(gene, levels = genes)) }
  # filter non-selected clusters
  if (!is.null(clusters)) { df <- filter(df, cluster %in% clusters) }
  # plot
  p <- ggplot(df, aes(x = !!sym(x), y = y, group = group))  +
    theme_classic() +
    scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) +
    xlab(x) + ylab("Normalized expression") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
    scale_fill_manual(values = colors)
  if (is.null(fill)) {
    p <- p + geom_violin(scale = "width", aes(fill = !!sym(x)), show.legend = FALSE)
  } else {
    p <- p + geom_violin(scale = "width", aes(fill = !!sym(fill)))
  }
  if (is.null(facet)) {
    p <- p + facet_wrap(~ gene, scales = "free_y", ncol = n_col)
  } else {
    p <- p + facet_wrap(~ vars(!!sym(facet), gene), scales = "free_y", ncol = n_col)
  }
  if (jitter) {
    p <- p + geom_jitter(width = 0.4, height = 0, alpha = 0.4, stroke = 0)
  }
  if (void) { p <- p + theme_void() }
  if (flip) { p <- p + coord_flip() }
  if (stacked) {
    p <- p + facet_wrap(~ gene, scales = "free_y", ncol = 1,
                        strip.position = "left") +
      scale_y_continuous(position = "right", limits = c(0, NA),
                         expand = c(0, 0), 
                         breaks = seq(2, max(df$y), by = 2),
                         labels = seq(2, max(df$y), by = 2)) +
      theme(strip.background = element_blank(),
            axis.ticks.y = element_blank(),
            strip.text.y = element_text(angle = 180,
                                        hjust = 1, vjust = 0.5)) +
      ylab(NULL)
  }
  return(p)
}


# - UMAP plot ----------------------------------------------------------------
umap_plot <- function(object, genes, cells = NULL, clusters = NULL, 
                      legend = TRUE, cluster_label = FALSE,
                      ncol = NULL, xlim = NULL, ylim = NULL,
                      order_genes = FALSE) {
  
  if (is.null(clusters)) {
    clusters <- sort(unique(object@active.ident))
  } else {
    clusters <- clusters[object@active.ident %in% clusters]
  }
  
  # pull expression data
  if (sum(!genes %in% rownames(object@assays$RNA@data))) {
    wrong_genes <- genes[!genes %in% rownames(object@assays$RNA@data)]
    warning(paste(
      "Warning:", paste(wrong_genes, collapse = ", "), "not in dataset."
    ))
    genes <- genes[!genes %in% wrong_genes]
  }
  pull_data <- function(genes) {
    if (length(genes) == 1) {
      df <- object@assays$RNA@data[genes, object@active.ident %in% clusters] %>% 
        as.data.frame() %>% 
        set_names(genes)
    } else {
      df <- object@assays$RNA@data[genes, object@active.ident %in% clusters] %>% 
        as.matrix() %>% t() %>% as.data.frame()
    }
    return(df)
  }
  
  # pull UMAP data
  umap <- data.frame(
    UMAP1 = object@reductions$umap@cell.embeddings[, 1],
    UMAP2 = object@reductions$umap@cell.embeddings[, 2]
  )
  umap <- umap[object@active.ident %in% clusters, ]

  
  # combine umap and expression data
  results <- pull_data(genes) %>%
    bind_cols(., umap) %>%
    gather(-starts_with("UMAP"), key = "gene", value = "value")
  
  if (order_genes) {
    results <- mutate(results, gene = factor(gene, levels = genes))
  }
  
  # standardize fill
  results <- results %>%
    group_by(gene) %>%
    mutate(value = value / max(value)) %>%
    ungroup()
  
  # plot
  plt <-
    ggplot(results, aes(x = UMAP1, y = UMAP2)) +
    geom_point(aes(color = value), show.legend = legend, stroke = 0) +
    scale_color_gradient(low = "gray90", high = "navyblue",
                         name = expression(underline("Expression"))) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          legend.text = element_blank()) +
    facet_wrap(~gene, ncol = ncol) +
    xlab("UMAP1") + ylab("UMAP2")
  if (cluster_label == TRUE) {
    cluster_center <- mutate(umap, "cluster" = object@active.ident)
    cluster_center <- cluster_center %>%
      group_by(cluster) %>%
      summarize(UMAP1 = mean(UMAP1), UMAP2 = mean(UMAP2))
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
  # grab data
  df <- data.frame(
    "UMAP1" = object@reductions$umap@cell.embeddings[,1],
    "UMAP2" = object@reductions$umap@cell.embeddings[,2],
    "cell" = names(object@active.ident)
  ) %>%
    mutate("Doublet" = ifelse(cell %in% doublets, TRUE, FALSE)) %>%
    arrange(Doublet) %>%
    mutate(cell = factor(cell, levels = cell))
  
  # add legend
  legend_pos <- c(min(df$UMAP1), max(df$UMAP2))
  
  # generate plot
  plt <- ggplot(df, aes(x = UMAP1, y = UMAP2, color = Doublet)) +
    geom_point(show.legend = FALSE, stroke = 0) +
    scale_color_manual(values = c("gray90", "firebrick3")) +
    annotate("text", x = legend_pos[1], y = legend_pos[2], 
             label = "Doublets", color = "firebrick3", hjust = 0) +
    theme_bw() +
    theme(panel.grid = element_blank())
  return(plt)
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

# - Plot GSEA scores ----------------------------------------------------------
plot_gsea_scores <- function(gsea_scores, max_val = 1, zeros = FALSE) {
  # set levels
  group_levels <- gsea_scores %>% as.data.frame() %>% 
    rownames_to_column("group") %>%
    gather(-group, key = "cluster", value = "score") %>%
    mutate_at(vars(cluster, score), as.numeric) %>%
    filter(score > 0) %>%
    mutate("place" = cluster * score) %>%
    group_by(group) %>%
    summarize("place" = median(place)) %>%
    arrange(place) %>%
    .$group
  cluster_levels <- colnames(gsea_scores)
  
  df <- gsea_scores %>% as.data.frame() %>% rownames_to_column("group") %>%
    filter(group %in% group_levels) %>%
    gather(-group, key = "cluster", value = "score") %>%
    mutate(group = factor(group, levels = group_levels),
           cluster = factor(cluster, levels = cluster_levels)) %>%
    mutate(score = ifelse(score < -max_val, -max_val,
                          ifelse(score > max_val, max_val, score))
    )
  # plot
  p <- ggplot(df, aes(x = cluster, y = fct_rev(group), fill = score)) +
    geom_tile(show.legend = FALSE) +
    scale_fill_gradient2(low = "#d0587e", mid = "white", high = "#009392",
                         name = "Score", na.value = "white") +
    theme_void() +
    theme(axis.text = element_text(color = "black"),
          axis.text.y = element_text(hjust = 1))
  return(p)
}


# - RGB plot -----------------------------------------------------------------
rgb_plot <- function(object, red = NULL, green = NULL, blue = NULL,
                     assay = "RNA", data = "data") {
  get_data <- function(color) {
    if (is.null(color)) {
      return(0)
    } else {
      d <- slot(slot(object, "assays")[[assay]], data)[color, ]
      d <- d / max(d)
    }
  }
  
  df <- data.frame(
    "UMAP1" = object@reductions$umap@cell.embeddings[,1],
    "UMAP2" = object@reductions$umap@cell.embeddings[,2],
    "red" = get_data(red),
    "green" = get_data(green),
    "blue" = get_data(blue)
  ) %>%
    mutate('cell' = factor(seq(nrow(.)))) %>%
    mutate("color" = rgb(red, green, blue))
  
  # plot
  p <- ggplot(df, aes(x = UMAP1, y = UMAP2, color = cell)) +
    scale_color_manual(values = df$color) +
    geom_point(show.legend = FALSE) +
    theme_classic() +
    theme(panel.background = element_rect(fill = "gray40", color = "black"))
  return(p)
}
