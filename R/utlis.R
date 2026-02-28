#' @keywords internal
validate_inputs <- function(clean_tab,
                            dirty_tab,
                            clean_seqs,
                            dirty_seqs) {

  clean_tab <- as.matrix(clean_tab)
  dirty_tab <- as.matrix(dirty_tab)

  storage.mode(clean_tab) <- "double"
  storage.mode(dirty_tab) <- "double"

  if (anyNA(rownames(clean_tab)) || anyNA(rownames(dirty_tab))) {
    stop("ASV IDs (rownames) cannot contain NA.")
  }
  if (anyDuplicated(rownames(clean_tab)) || anyDuplicated(rownames(dirty_tab))) {
    stop("ASV IDs (rownames) must be unique (no duplicates).")
  }

  if (is.null(rownames(clean_tab)) || is.null(rownames(dirty_tab))) {
    stop("Both clean_tab and dirty_tab must have rownames (ASV IDs).")
  }

  if (is.null(colnames(clean_tab)) || is.null(colnames(dirty_tab))) {
    stop("Both clean_tab and dirty_tab must have colnames (sample IDs).")
  }

  if (!setequal(colnames(clean_tab), colnames(dirty_tab))) {
    stop("Sample columns must match between CLEAN and DIRTY tables.")
  }

  samp <- colnames(clean_tab)  # preserve CLEAN order
  dirty_tab <- dirty_tab[, samp, drop = FALSE]

  # ---- Sequence validation ----
  if (is.null(names(clean_seqs)) || is.null(names(dirty_seqs))) {
    stop("clean_seqs and dirty_seqs must be named character vectors.")
  }

  if (!all(rownames(clean_tab) %in% names(clean_seqs))) {
    stop("All CLEAN ASV IDs must exist in clean_seqs.")
  }

  if (!all(rownames(dirty_tab) %in% names(dirty_seqs))) {
    stop("All DIRTY ASV IDs must exist in dirty_seqs.")
  }

  # Reorder sequences to match table row order
  clean_seqs <- stats::setNames(as.character(clean_seqs[rownames(clean_tab)]),
                         rownames(clean_tab))
  dirty_seqs <- stats::setNames(as.character(dirty_seqs[rownames(dirty_tab)]),
                         rownames(dirty_tab))

  list(
    clean_tab   = clean_tab,
    dirty_tab   = dirty_tab,
    clean_seqs  = clean_seqs,
    dirty_seqs  = dirty_seqs
  )
}
