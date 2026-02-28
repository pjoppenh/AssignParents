#' @keywords internal
compute_single_scores <- function(presence_mat,
                                  abundance_mat,
                                  weights = c(presence = 0.6,
                                              abundance = 0.4),
                                  scale_method = c("z","robust","rank")) {

  scale_method <- match.arg(scale_method)

  if (!identical(dim(presence_mat), dim(abundance_mat))) {
    stop("presence_mat and abundance_mat must have identical dimensions.")
  }

  n_dirty <- nrow(presence_mat)
  n_clean <- ncol(presence_mat)

  single_mat <- matrix(NA_real_,
                       nrow = n_dirty,
                       ncol = n_clean,
                       dimnames = dimnames(presence_mat))

  # normalize weights
  w <- weights / sum(weights)
  w_p <- w["presence"]
  w_a <- w["abundance"]

  for (i in seq_len(n_dirty)) {

    p_row <- presence_mat[i, ]
    a_row <- abundance_mat[i, ]

    # row-wise standardization (same logic used in original)
    zp <- .row_center_scale(p_row, method = scale_method)
    za <- .row_center_scale(a_row, method = scale_method)

    for (j in seq_len(n_clean)) {

      if (is.na(zp[j]) && is.na(za[j])) next

      wp <- if (!is.na(zp[j])) w_p else 0
      wa <- if (!is.na(za[j])) w_a else 0
      wsum <- wp + wa

      if (wsum == 0) next

      single_mat[i, j] <- (wp * zp[j] + wa * za[j]) / wsum
    }
  }

  single_mat
}

#' @keywords internal
.row_center_scale <- function(x, method = c("z","robust","rank")) {
  method <- match.arg(method)
  n <- length(x)
  if (n == 0) return(x)
  if (all(is.na(x))) return(rep(NA_real_, n))

  if (method == "z") {
    m <- mean(x, na.rm = TRUE)
    s <- stats::sd(x, na.rm = TRUE)
    if (is.na(s) || !is.finite(s) || s == 0) return(rep(0, n))
    return((x - m) / s)
  }

  if (method == "robust") {
    m <- stats::median(x, na.rm = TRUE)
    s <- stats::mad(x, na.rm = TRUE, constant = 1.4826)
    if (is.na(s) || !is.finite(s) || s == 0) return(rep(0, n))
    return((x - m) / s)
  }

  r <- rank(x, na.last = "keep", ties.method = "average")
  m <- mean(r, na.rm = TRUE)
  s <- stats::sd(r, na.rm = TRUE)
  if (is.na(s) || !is.finite(s) || s == 0) return(rep(0, n))
  (r - m) / s
}

#' @keywords internal
compute_score_matrices <- function(clean_tab,
                                   dirty_tab,
                                   identity_mat,
                                   min_identity,
                                   presence_cutoff,
                                   presence_weights,
                                   abundance_weights,
                                   min_recall,
                                   min_n,
                                   enforce_dominance,
                                   dominance_mean_ratio_max,
                                   dominance_frac_exceed_max,
                                   dominance_min_both,
                                   par) {

  clean_ids <- rownames(clean_tab)
  dirty_ids <- rownames(dirty_tab)

  lib_clean <- colSums(clean_tab)
  lib_dirty <- colSums(dirty_tab)
  eps <- .Machine$double.eps

  # ---- Preallocate ----
  n_dirty <- nrow(dirty_tab)
  n_clean <- nrow(clean_tab)

  presence_score_mat  <- matrix(NA_real_, n_dirty, n_clean, dimnames = list(dirty_ids, clean_ids))
  abundance_score_mat <- matrix(NA_real_, n_dirty, n_clean, dimnames = list(dirty_ids, clean_ids))

  recall_mat    <- matrix(NA_real_, n_dirty, n_clean, dimnames = list(dirty_ids, clean_ids))
  phi01_mat     <- matrix(NA_real_, n_dirty, n_clean, dimnames = list(dirty_ids, clean_ids))
  ppmi_tanh_mat <- matrix(NA_real_, n_dirty, n_clean, dimnames = list(dirty_ids, clean_ids))
  f1_mat        <- matrix(NA_real_, n_dirty, n_clean, dimnames = list(dirty_ids, clean_ids))

  spearman01_mat <- matrix(NA_real_, n_dirty, n_clean, dimnames = list(dirty_ids, clean_ids))
  hellinger_mat  <- matrix(NA_real_, n_dirty, n_clean, dimnames = list(dirty_ids, clean_ids))

  dominance_mean_ratio_mat  <- matrix(NA_real_, n_dirty, n_clean, dimnames = list(dirty_ids, clean_ids))
  dominance_frac_exceed_mat <- matrix(NA_real_, n_dirty, n_clean, dimnames = list(dirty_ids, clean_ids))
  dominance_n_both_mat      <- matrix(NA_integer_, n_dirty, n_clean, dimnames = list(dirty_ids, clean_ids))
  dominance_pass_mat        <- matrix(NA, n_dirty, n_clean, dimnames = list(dirty_ids, clean_ids))

  compute_dirty_row <- function(i) {

    d_counts <- dirty_tab[i, ]

    pres_score <- rep(NA_real_, n_clean)
    abund_score <- rep(NA_real_, n_clean)

    recall_v <- rep(NA_real_, n_clean)
    phi01_v  <- rep(NA_real_, n_clean)
    ppmi_v   <- rep(NA_real_, n_clean)
    f1_v     <- rep(NA_real_, n_clean)

    sp_v  <- rep(NA_real_, n_clean)
    hel_v <- rep(NA_real_, n_clean)

    dom_mean <- rep(NA_real_, n_clean)
    dom_frac <- rep(NA_real_, n_clean)
    dom_nboth <- rep(NA_integer_, n_clean)
    dom_pass <- rep(NA, n_clean)

    cand <- which(!is.na(identity_mat[i, ]) & identity_mat[i, ] >= min_identity)
    if (length(cand) == 0) {
      return(list(
        pres_score = pres_score,
        abund_score = abund_score,
        recall = recall_v,
        phi01 = phi01_v,
        ppmi_tanh = ppmi_v,
        f1 = f1_v,
        spearman01 = sp_v,
        hellinger = hel_v,
        dom_mean = dom_mean,
        dom_frac = dom_frac,
        dom_nboth = dom_nboth,
        dom_pass = dom_pass
      ))
    }

    for (j in cand) {

      p_counts <- clean_tab[j, ]

      pres <- presence_composite(
        p_counts, d_counts,
        presence_cutoff = presence_cutoff,
        weights = presence_weights,
        min_recall = min_recall,
        min_n = min_n
      )

      if (!isTRUE(pres$passed_gate)) next

      pres_score[j] <- pres$score
      recall_v[j]   <- pres$components[["recall"]]
      phi01_v[j]    <- pres$components[["phi01"]]
      ppmi_v[j]     <- pres$components[["ppmi_tanh"]]
      f1_v[j]       <- pres$components[["f1"]]

      abund <- abundance_composite(p_counts, d_counts, weights = abundance_weights)

      abund_score[j] <- abund$score
      sp_v[j]        <- abund$components[["spearman"]]
      hel_v[j]       <- abund$components[["hellinger"]]

      both <- (p_counts > presence_cutoff) & (d_counts > presence_cutoff)
      n_both <- sum(both, na.rm = TRUE)
      dom_nboth[j] <- as.integer(n_both)

      if (n_both > 0) {
        p_rel <- p_counts / pmax(lib_clean, eps)
        d_rel <- d_counts / pmax(lib_dirty, eps)
        ratio_norm  <- d_rel[both] / p_rel[both]
        mean_ratio  <- mean(ratio_norm, na.rm = TRUE)
        frac_exceed <- mean(d_rel[both] > p_rel[both], na.rm = TRUE)
      } else {
        mean_ratio <- NA_real_
        frac_exceed <- NA_real_
      }

      dom_mean[j] <- mean_ratio
      dom_frac[j] <- frac_exceed

      if (n_both >= dominance_min_both) {
        dom_pass[j] <- (is.finite(mean_ratio) && mean_ratio <= dominance_mean_ratio_max) &&
          (is.finite(frac_exceed) && frac_exceed <= dominance_frac_exceed_max)
      } else {
        dom_pass[j] <- NA
      }
    }

    list(
      pres_score = pres_score,
      abund_score = abund_score,
      recall = recall_v,
      phi01 = phi01_v,
      ppmi_tanh = ppmi_v,
      f1 = f1_v,
      spearman01 = sp_v,
      hellinger = hel_v,
      dom_mean = dom_mean,
      dom_frac = dom_frac,
      dom_nboth = dom_nboth,
      dom_pass = dom_pass
    )
  }

  idx_all <- seq_len(n_dirty)
  chunk_size <- 50L
  chunks <- split(idx_all, ceiling(idx_all / chunk_size))

  row_results <- vector("list", n_dirty)

  pb <- utils::txtProgressBar(min = 0, max = length(chunks), style = 3)

  for (k in seq_along(chunks)) {
    idx <- chunks[[k]]
    res_k <- par$par_lapply(idx, compute_dirty_row)
    for (m in seq_along(idx)) row_results[[ idx[m] ]] <- res_k[[m]]
    utils::setTxtProgressBar(pb, k)
  }
  close(pb)

  for (i in seq_len(n_dirty)) {
    r <- row_results[[i]]
    presence_score_mat[i, ]  <- r$pres_score
    abundance_score_mat[i, ] <- r$abund_score
    recall_mat[i, ]          <- r$recall
    phi01_mat[i, ]           <- r$phi01
    ppmi_tanh_mat[i, ]       <- r$ppmi_tanh
    f1_mat[i, ]              <- r$f1
    spearman01_mat[i, ]      <- r$spearman01
    hellinger_mat[i, ]       <- r$hellinger
    dominance_mean_ratio_mat[i, ]  <- r$dom_mean
    dominance_frac_exceed_mat[i, ] <- r$dom_frac
    dominance_n_both_mat[i, ]      <- r$dom_nboth
    dominance_pass_mat[i, ]        <- r$dom_pass
  }

  list(
    presence_score = presence_score_mat,
    abundance_score = abundance_score_mat,
    recall = recall_mat,
    phi01 = phi01_mat,
    ppmi_tanh = ppmi_tanh_mat,
    f1 = f1_mat,
    spearman01 = spearman01_mat,
    hellinger = hellinger_mat,
    dominance_mean_ratio = dominance_mean_ratio_mat,
    dominance_frac_exceed = dominance_frac_exceed_mat,
    dominance_n_both = dominance_n_both_mat,
    dominance_pass = dominance_pass_mat
  )
}
