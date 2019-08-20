# -- Matrix functions ---------------------------------------------------------
# keep genes that are in all matrices
keep_shared_genes <- function(mtx_list) {
  shared <- map(mtx, rownames) %>% 
    unlist() %>% 
    table() %>% 
    .[. == length(mtx)] %>% 
    names()
  return(map(mtx, ~ .x[shared, ]))
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
remove_unwanted_genes <- function(mtx, keep_Gm = FALSE) {
  # remove genes from mitochondrial genome
  nonmito <- read_csv("~/Programs/dropseq3/data/nonmito_genes.csv",
                      col_names = FALSE) %>% .$X1
  genes <- unique(unlist(map(mtx, rownames)))
  nonmito <- nonmito[nonmito %in% genes]
  if (class(mtx) == "list") {
    mtx <- map(mtx, ~ .x[nonmito, ])
  } else if (class(mtx) == "dgCMatrix") {
    mtx <- mtx[nonmito, ]
  }
  
  # remove genes starting with Gm if desired
  if (keep_Gm == FALSE) {
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
  abund <- sapply(mtx, function(x) Matrix::rowSums(x > 0))
  if (class(mtx) == "list") {
    abund <- rowSums(abund)
  }
  
  # select only genes above
  keep <- names(abund)[abund >= min_cells]
  if (class(mtx) == "list") {
    mtx <- map(mtx, ~ .x[keep, ])
  } else if (class(mtx) == "dgCMatrix") {
    mtx <= mtx[keep, ]
  }
  return(mtx)
}

# renaming barcodes that are duplicated across experiments
rename_duplicates <- function(mtx) {
  if (sum(duplicated(colnames(do.call(cbind, mtx)))) > 0) {
    new_names <- function(dim) {
      other <- colnames(unlist(mtx[[-dim]]))
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

# - Remove doublets ----------------------------------------------------------
remove_doublets <- function(object, markers) {
  
}
