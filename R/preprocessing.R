library(reticulate)

# -- Matrix functions ---------------------------------------------------------
# downsample matrix
downsample_matrix <- function(mtx, n_counts = 1000) {
  gene_levels <- rownames(mtx); cell_levels <- colnames(mtx)
  counts <- map(1:ncol(mtx), ~ rep(gene_levels, times = mtx[, .x]))
  new_counts <- map(counts,
                    ~ sample(.x, min(n_counts, length(.x)), replace = FALSE))
  new_counts <- map(new_counts, ~ table(factor(.x, levels = gene_levels)))
  mtx <- do.call(cbind, new_counts)
  colnames(mtx) <- cell_levels
  return(Matrix::Matrix(mtx, sparse = TRUE))
}

# collapse duplicate genes
collapse_duplicate_genes <- function(mtx, method = "max") {
  if (!any(duplicated(rownames(mtx)))) { return(mtx) }
  duplicate_genes <- unique(rownames(mtx)[duplicated(rownames(mtx))])
  mtx1 <- mtx[-which(rownames(mtx) %in% duplicate_genes), ]
  if (method == "max") {
    # take row that has max value
    max_function <- function(gene) {
      indexes <- which(rownames(mtx) == gene)
      values <- Matrix::rowSums(mtx[indexes, ])
      return(indexes[values == max(values)][1])
    }
    mtx2 <- mtx[sapply(duplicate_genes, max_function, simplify = TRUE), ]
  } else if (method == "combine") {
    # combine all rows together
    combine_function <- function(gene) {
      indexes <- which(rownames(mtx) == gene)
      return(Matrix::colSums(mtx[indexes, ]))
    }
    mtx2 <- t(sapply(duplicate_genes, combine_function))
  } else if (method == "first") {
    # take first instance of gene
    mtx2 <- t(sapply(duplicate_genes, function(x) mtx[x, ]))
  }
  mtx2 <- Matrix::Matrix(mtx2, sparse = TRUE)
  gc(verbose = FALSE)
  return(rbind(mtx1, mtx2))
}

# keep genes that are in all matrices
keep_shared_genes <- function(mtx_list, duplicate_gene_method = "max") {
  duplicated_genes <- any(map_int(mtx_list, ~ sum(duplicated(rownames(.x)))))
  if (duplicated_genes) {
    mtx_list <- map(mtx_list,
                    ~ collapse_duplicate_genes(.x, method = duplicate_gene_method)
                    )
  }
  shared <- unlist(sapply(mtx_list, rownames)) %>% table() %>%
    .[. == length(mtx_list)] %>% names()
  return(map(mtx_list, ~ .x[shared, ]))
}

# remove genes that are not expressed in any cells
remove_zeros <- function(mtx) {
  if (class(mtx) == "list") {
    mtx2 <- do.call(cbind, mtx)
    expr <- Matrix::rowSums(mtx2)
    expr <- names(expr)[expr > 0]
    return(map(mtx, ~ .x[expr, ]))
  } else if (class(mtx) == "dgCMatrix") {
    expr <- Matrix::rowSums(mtx)
    expr <- names(expr)[expr > 0]
    return(mtx[expr, ])
  }
}

# remove unwanted genes
remove_unwanted_genes <- function(mtx, remove_Gm = TRUE) {
  # remove genes from mitochondrial genome
  nonmito <- read_csv("~/Programs/dropseq3/data/nonmito_genes.csv",
                      col_names = FALSE, col_types = "c") %>% .$X1
  if (class(mtx) == "list") {
    genes <- unique(unlist(map(mtx, rownames)))
    nonmito <- nonmito[nonmito %in% genes]
    mtx <- map(mtx, ~ .x[nonmito, ])
  } else if (class(mtx) == "dgCMatrix") {
    nonmito <- nonmito[nonmito %in% rownames(mtx)]
    mtx <- mtx[nonmito, ]
  }

  # remove genes starting with Gm if desired
  if (remove_Gm) {
    if (class(mtx) == "list") {
      mtx <- map(mtx, ~ .x[!str_detect(rownames(.x), "^Gm"), ])
    } else if (class(mtx) == "dgCMatrix") {
      mtx <- mtx[!str_detect(rownames(mtx), "^Gm"), ]
    }
  }
  return(mtx)
}

# remove low abundance genes
remove_low_abundance_genes <- function(mtx, min_cells = 4) {
  # calculate number of cells expressing a given gene
  if (class(mtx) == "list") {
    abund <- apply(sapply(mtx, function(x) Matrix::rowSums(x > 0)), 1, sum)
    keep <- rownames(mtx[[1]])[abund >= min_cells]
    mtx <- map(mtx, ~ .x[keep, ])
  } else if (class(mtx) == "dgCMatrix") {
    abund <- Matrix::rowSums(mtx > 0)
    keep <- names(abund)[abund >= min_cells]
    mtx <- mtx[keep, ]
  }
  return(mtx)
}

# filter
filter_cells <- function(mtx, min_genes = 500) {
  if (class(mtx) == "list") {
    mtx <- map(mtx, ~ .x[, Matrix::colSums(.x > 0) > min_genes])
  } else if (class(mtx) == "dgCMatrix") {
    mtx <- mtx[, Matrix::colSums(mtx > 0) > min_genes]
  }
  return(mtx)
}

# renaming barcodes that are duplicated across experiments
rename_duplicates <- function(mtx) {
  if (sum(duplicated(colnames(do.call(cbind, mtx)))) > 0) {
    new_names <- function(dim) {
      other <- unlist(map(mtx[-1], colnames))
      current <- colnames(mtx[[dim]])
      current <- ifelse(current %in% other, paste0(current, "-", dim), current)
      colnames(mtx[[dim]]) <- current
      return(mtx[[dim]])
    }
    mtx <- map(seq(length(mtx)), new_names)
  }
  return(mtx)
}

# calculate percent mitochondrial genes
percent_mito <- function(object) {
  mito_genes <- str_detect(rownames(object@assays$RNA@counts), "^mt")
  mito <-
    Matrix::colSums(object@assays$RNA@counts[mito_genes, ]) /
    Matrix::colSums(object@assays$RNA@counts)
  object@meta.data$percent_mito <- mito
  return(object)
}

# - Find variable genes -----------------------------------------------------
find_variable_genes <- function(object, genes = NULL, assay = "RNA",
                                data = "data", n_genes = 2000,
                                min_avg = 0.1, span = 0.3,
                                plot = FALSE) {
  # get data matrix
  mtx <- slot(object[[assay]], data)
  # keep only selected genes
  if (!is.null(genes)) {
    if (sum(!genes %in% rownames(mtx)) > 0) {
      warning(paste("Genes not in dataset:",
                    paste(genes[!genes %in% rownames(mtx)], collapse = ", ")))
    }
    genes <- genes[genes %in% rownames(mtx)]
    mtx <- mtx[genes, ]
  }
  # find mean and variance
  avg <- Matrix::rowMeans(mtx)
  above <- avg > min_avg
  df <- data.frame("avg" = avg[above],
                   "disp" = apply(mtx[above, ], 1, var) / avg[above])
  # use loess to get smoothed average
  model <- loess(disp ~ avg, data = df, span = span)
  var_genes <- model$residuals %>% sort(decreasing = TRUE)
  var_genes <- names(var_genes)[1:min(n_genes, length(var_genes))]
  slot(object[[assay]], "var.features") <- var_genes
  return(object)
}

# - Remove doublets using Scrublet ---------------------------------------------
scrublet <- function(mtx) {
  source_python("/home/alanrupp/Programs/dropseq3/python/scrublet.py")
  expected_doublet_rate <- find_expected_doublet_rate(mtx)
  doublets <- score_doublets(mtx, expected_doublet_rate)
  gc(verbose = FALSE)
  return(doublets)
}

# - Integrate data ------------------------------------------------------------
integrate_data <- function(object) {
  k <- min(200, min(sapply(object, ncol)))
  print("Finding integration anchors ...")
  anchors <- FindIntegrationAnchors(object, k.filter = k, verbose = FALSE)
  print("Integrating data ...")
  object <- IntegrateData(anchors, verbose = FALSE)
  DefaultAssay(object) <- "integrated"
  return(object)
}

# - Run scran normalize -------------------------------------------------------
run_scran_normalize <- function(object) {
  print(paste("Normalizing data with scran", packageVersion("scran"), "..."))
  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list("counts" = object@assays$RNA@counts)
    )
  print("Running dimension reduction ...")
  clusters <- scran::quickCluster(sce, min.size = 0)
  print("Computing sum factors ...")
  sce <- scran::computeSumFactors(sce, clusters = clusters)
  print("Log normalizing counts ...")
  sce <- scater::logNormCounts(sce)
  object@assays$RNA@data <- SingleCellExperiment::logcounts(sce)
  return(object)
}

# - Normalize data ----------------------------------------------------------
normalize_data <- function(object, method = "scran", batch = NULL) {
  if (method == "scran") {
    if (!is.null(batch)) {
      clusters <- levels(object@active.ident)
      object <- SplitObject(object, batch)
      object <- map(object, run_scran_normalize)
      object <- merge(object[[1]], object[2:length(object)])
      object@active.ident <- factor(object@active.ident, levels = clusters)
    } else {
      object <- run_scran_normalize(object)
    }
  } else if (method == "sctransform") {
    library(sctransform)
    metadata <- object@meta.data
    metadata$log_umi_per_gene <- log10(metadata$nCount_RNA/metadata$nFeature_RNA)
    vst_out <- vst(object@assays$RNA@counts,
                   cell_attr = metadata,
                   latent_var = 'log_umi_per_gene',
                   batch_var = batch)
    object@assays$RNA@data <- vst_out$y
  } else if (method == "TPM") {
    # median UMI counts
    mtx <- slot(object@assays[[assay]], "counts")
    med <- median(Matrix::colSums(mtx))
    mtx <- log1p(apply(mtx, 2, function(x) x / sum(x)) * med)
    object@assays$RNA@data <- Matrix::Matrix(mtx, sparse = TRUE)
  }
  gc(verbose = FALSE)
  return(object)
}

# - Scale data --------------------------------------------------------------
scale_data <- function(object, groups = NULL, assay = "RNA", data = "data",
                       genes = NULL) {
  mtx <- slot(object@assays[[assay]], data)
  if (!is.null(genes)) mtx <- mtx[genes, ]
  if (is.null(groups)) {
    scaled <- t(scale(Matrix::t(mtx)))
    rownames(scaled) <- rownames(mtx); colnames(scaled) <- colnames(mtx)
  } else {
    scaled <- map(
      unique(object@meta.data[, groups]),
      ~ t(scale(Matrix::t(mtx[, which(object@meta.data[, groups] == .x)])))
      )
    scaled <- do.call(cbind, scaled)
    scaled <- scaled[, colnames(mtx)]
  }
  gc(verbose = FALSE)
  object@assays[[assay]]@scale.data <- scaled
  return(object)
}

# - PCA knee test -----------------------------------------------------------
knee_test <- function(object) {
  n_pc = ncol(object@reductions$pca@cell.embeddings)
  total_var <- sum(object@reductions$pca@stdev^2)
  percent_var <- cumsum(object@reductions$pca@stdev^2)/total_var * n_pc
  diminishing <- which(percent_var - lag(percent_var) < 1)
  return(min(diminishing) - 1)
}

# - Quantile normalize --------------------------------------------------------
quantile_normalize <- function(df) {
  # code adapted from Dave Tang
  ranks <- apply(df, 2, rank, ties.method = "min")
  means <- apply(df, 2, sort) %>% apply(., 1, mean)
  get_values <- function(rank_values) {
    values <- means[rank_values]
    # if ties, average each value by its position
    if (any(duplicated(rank_values))) {
      dups <- unique(rank_values[duplicated(rank_values)])
      new_means <- sapply(dups, function(x) {
        missing_ranks <- seq(x, x+sum(rank_values == x)-1)
        mean(means[missing_ranks])
      })
      for (i in 1:length(dups)) {
        values[rank_values == dups[i]] <- new_means[i]
      }
    }
    return(values)
  }
  new_values <- apply(ranks, 2, get_values)
  rownames(new_values) <- rownames(df); colnames(new_values) <- colnames(df)
  return(new_values)
}
