#' @keywords internal
compute_identity_matrix <- function(clean_tab,
                                    dirty_tab,
                                    clean_seqs,
                                    dirty_seqs,
                                    min_identity,
                                    par,
                                    chunk_size = 50L) {

  clean_ids <- rownames(clean_tab)
  dirty_ids <- rownames(dirty_tab)

  clean_seq_vec <- as.character(clean_seqs[clean_ids])
  dirty_seq_vec <- as.character(dirty_seqs[dirty_ids])

  if (anyNA(clean_seq_vec)) {
    bad <- clean_ids[is.na(clean_seq_vec)]
    stop("compute_identity_matrix(): NA CLEAN sequences for IDs: ",
         paste(bad, collapse = ", "))
  }
  if (anyNA(dirty_seq_vec)) {
    bad <- dirty_ids[is.na(dirty_seq_vec)]
    stop("compute_identity_matrix(): NA DIRTY sequences for IDs: ",
         paste(bad, collapse = ", "))
  }

  clean_len <- nchar(clean_seq_vec)

  # ---- Length index for CLEAN sequences ----
  clean_len_to_idx <- split(seq_along(clean_len), clean_len)

  .get_clean_candidates_by_len <- function(dl, max_len_diff) {
    lo <- dl - max_len_diff
    hi <- dl + max_len_diff
    lens <- as.integer(names(clean_len_to_idx))
    keep <- lens[lens >= lo & lens <= hi]
    if (length(keep) == 0) return(integer(0))
    unlist(clean_len_to_idx[as.character(keep)], use.names = FALSE)
  }

  message("Computing identity matrix...")

  all_len <- unique(nchar(c(clean_seq_vec, dirty_seq_vec)))
  equal_len <- length(all_len) == 1

  # -------------------------------------------------------------
  # Identity computation
  # -------------------------------------------------------------

  if (equal_len) {

    clean_split <- strsplit(clean_seq_vec, "")

    hamming_id_row <- function(i) {
      s1c <- strsplit(dirty_seq_vec[i], "")[[1]]
      vapply(clean_split, function(x) mean(x == s1c), numeric(1))
    }

    identity_list <- par$par_lapply(seq_along(dirty_seq_vec), hamming_id_row)
    identity_mat  <- do.call(rbind, identity_list)

  } else {

    lev_id_row <- function(i) {

      di <- dirty_seq_vec[i]
      dl <- nchar(di)

      maxlen_screen      <- max(dl, max(clean_len))
      max_edits_screen   <- ceiling((1 - min_identity) * maxlen_screen)

      cand_idx <- .get_clean_candidates_by_len(dl, max_edits_screen)
      out <- rep(NA_real_, length(clean_len))

      if (length(cand_idx) == 0) return(out)

      maxlen     <- pmax(dl, clean_len[cand_idx])
      max_edits  <- ceiling((1 - min_identity) * maxlen)
      feas_local <- abs(dl - clean_len[cand_idx]) <= max_edits
      if (!any(feas_local)) return(out)

      cand2   <- cand_idx[feas_local]
      maxlen2 <- maxlen[feas_local]

      dists <- stringdist::stringdist(di, clean_seq_vec[cand2], method = "lv")
      out[cand2] <- 1 - (dists / maxlen2)

      out
    }

    idx_all <- seq_along(dirty_seq_vec)

    chunk_size <- as.integer(chunk_size)
    if (is.na(chunk_size) || chunk_size < 1L) chunk_size <- 50L

    chunks <- split(idx_all, ceiling(idx_all / chunk_size))

    identity_rows <- vector("list", length(dirty_seq_vec))

    pb <- utils::txtProgressBar(min = 0, max = length(chunks), style = 3)
    on.exit(close(pb), add = TRUE)

    for (k in seq_along(chunks)) {
      idx   <- chunks[[k]]
      res_k <- par$par_lapply(idx, lev_id_row)

      for (m in seq_along(idx)) {
        identity_rows[[ idx[m] ]] <- res_k[[m]]
      }

      utils::setTxtProgressBar(pb, k)
    }

    identity_mat <- do.call(rbind, identity_rows)
  }

  # Guarantee dimnames
  identity_mat <- as.matrix(identity_mat)
  rownames(identity_mat) <- dirty_ids
  colnames(identity_mat) <- clean_ids

  message("Computing mismatch matrix...")

  # -------------------------------------------------------------
  # Mismatch computation (preserve dimnames)
  # -------------------------------------------------------------

  if (equal_len) {

    L <- all_len[1]

    mismatches_mat <- pmax(
      0L,
      as.integer(round((1 - identity_mat) * L))
    )

    mismatches_mat <- matrix(
      mismatches_mat,
      nrow = nrow(identity_mat),
      ncol = ncol(identity_mat),
      dimnames = dimnames(identity_mat)
    )

  } else {

    mismatches_mat <- matrix(
      NA_integer_,
      nrow = nrow(identity_mat),
      ncol = ncol(identity_mat),
      dimnames = dimnames(identity_mat)
    )

    for (i in seq_len(nrow(identity_mat))) {
      di_len <- nchar(dirty_seq_vec[i])
      maxlen <- pmax(di_len, clean_len)
      mismatches_mat[i, ] <- as.integer(round((1 - identity_mat[i, ]) * maxlen))
    }
  }

  list(
    identity    = identity_mat,
    mismatches  = mismatches_mat
  )
}
