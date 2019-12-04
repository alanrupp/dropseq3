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


# keep genes that are in all matrices
keep_shared_genes <- function(mtx_list) {
  shared <- sapply(mtx_list, rownames) %>% table() %>% 
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
find_variable_genes <- function(object) {
  df <- data.frame("avg" = Matrix::rowMeans(object@assays$RNA@data),
                   "disp" = apply(object@assays$RNA@data, 1, var) / 
                     Matrix::rowMeans(object@assays$RNA@data))
  model <- loess(disp ~ avg, data = df)
}

# - Remove doublets using Scrublet ---------------------------------------------
scrublet <- function(mtx) {
  source_python("/home/alanrupp/Programs/dropseq3/python/scrublet.py")
  expected_doublet_rate <- find_expected_doublet_rate(mtx)
  doublets <- score_doublets(mtx, expected_doublet_rate)
  return(doublets)
}


# - Integrate data ------------------------------------------------------------
integrate_data <- function(object) {
  print("Integrating data")
  k <- min(200, min(sapply(object, ncol)))
  anchors <- FindIntegrationAnchors(object, k.filter = k, verbose = FALSE)
  object <- IntegrateData(anchors, verbose = FALSE)
  DefaultAssay(object) <- "integrated"
  return(object)
}