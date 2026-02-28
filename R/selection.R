#' Select CLEAN parents for each DIRTY ASV (advanced neighborhood + z-margin logic)
#'
#' @keywords internal
select_parents <- function(
    identity_mat,
    mismatches_mat,
    single_score_mat,
    presence_score_mat,
    abundance_score_mat,
    recall_mat,
    phi01_mat,
    ppmi_tanh_mat,
    f1_mat,
    spearman01_mat,
    hellinger_mat,
    dominance_mean_ratio_mat,
    dominance_frac_exceed_mat,
    dominance_n_both_mat,
    dominance_pass_mat,
    min_identity,
    min_presence_score,
    enforce_dominance,
    identity_neighborhood_mismatches,
    neighborhood_z_margin,
    clean_tab
) {

  if (is.null(rownames(identity_mat)) || is.null(colnames(identity_mat))) {
    stop("identity_mat must have rownames (DIRTY IDs) and colnames (CLEAN IDs).")
  }
  if (is.null(rownames(mismatches_mat)) || is.null(colnames(mismatches_mat))) {
    stop("mismatches_mat must have rownames/colnames matching identity_mat.")
  }

  dirty_ids <- rownames(identity_mat)
  clean_ids <- colnames(identity_mat)

  # identical parent totals tie-breaker (original logic)
  clean_totals <- rowSums(clean_tab)

  # NOTE: relies on .row_center_scale() existing in your package (from scoring.R)
  if (!exists(".row_center_scale", mode = "function")) {
    stop("Internal helper .row_center_scale() not found. Define it (e.g., in scoring.R).")
  }

  pick_parent <- function(did) {

    idvec    <- identity_mat[did, ]
    misvec   <- mismatches_mat[did, ]
    sscore   <- single_score_mat[did, ]
    pscore   <- presence_score_mat[did, ]
    ascore   <- abundance_score_mat[did, ]
    dom_pass <- dominance_pass_mat[did, ]

    # dominance filter (enforce only when we have evidence)
    dom_ok <- if (enforce_dominance) {
      is.na(dom_pass) | (dom_pass == TRUE)
    } else {
      rep(TRUE, length(dom_pass))
    }

    if (all(is.na(idvec))) {
      return(c(
        parent = NA_character_, identity = NA_real_,
        single_score = NA_real_, presence_score = NA_real_, abundance_score = NA_real_,
        recall = NA_real_, phi01 = NA_real_, ppmi_tanh = NA_real_, f1 = NA_real_,
        spearman01 = NA_real_, hellinger = NA_real_,
        dom_mean_ratio = NA_real_, dom_frac_exceed = NA_real_, dom_n_both = NA_real_
      ))
    }

    # Hard identity floor
    ok_id <- !is.na(idvec) & (idvec >= min_identity)
    if (!any(ok_id)) {
      return(c(
        parent = NA_character_, identity = max(idvec, na.rm = TRUE),
        single_score = NA_real_, presence_score = NA_real_, abundance_score = NA_real_,
        recall = NA_real_, phi01 = NA_real_, ppmi_tanh = NA_real_, f1 = NA_real_,
        spearman01 = NA_real_, hellinger = NA_real_,
        dom_mean_ratio = NA_real_, dom_frac_exceed = NA_real_, dom_n_both = NA_real_
      ))
    }

    eps <- 1e-12

    # ---- perfect-identity override (bypass dominance/presence gating) ----
    perfect <- which(ok_id & (idvec >= 1 - eps))
    if (length(perfect) > 0) {

      ss_sub <- sscore[perfect]

      if (all(is.na(ss_sub))) {
        totals_sub <- clean_totals[perfect]
        parent_idx <- perfect[which(totals_sub == max(totals_sub, na.rm = TRUE))]
        if (length(parent_idx) > 1) parent_idx <- parent_idx[order(names(idvec)[parent_idx])][1]
      } else {
        parent_idx <- perfect[which(ss_sub == max(ss_sub, na.rm = TRUE))]
        if (length(parent_idx) > 1) {
          totals_sub <- clean_totals[parent_idx]
          parent_idx <- parent_idx[which(totals_sub == max(totals_sub, na.rm = TRUE))]
          if (length(parent_idx) > 1) parent_idx <- parent_idx[order(names(idvec)[parent_idx])][1]
        }
      }

      parent <- names(idvec)[parent_idx]
      chosen_id <- idvec[parent_idx]
      chosen_ps <- pscore[parent_idx]
      chosen_as <- ascore[parent_idx]
      chosen_ss <- sscore[parent_idx]

      # still respect a positive presence threshold, if set
      if (!is.na(chosen_ps) && chosen_ps < min_presence_score) {
        return(c(
          parent = NA_character_, identity = chosen_id,
          single_score = chosen_ss, presence_score = chosen_ps, abundance_score = chosen_as,
          recall = NA_real_, phi01 = NA_real_, ppmi_tanh = NA_real_, f1 = NA_real_,
          spearman01 = NA_real_, hellinger = NA_real_,
          dom_mean_ratio  = dominance_mean_ratio_mat[did, parent_idx],
          dom_frac_exceed = dominance_frac_exceed_mat[did, parent_idx],
          dom_n_both      = dominance_n_both_mat[did, parent_idx]
        ))
      }

      return(c(
        parent = parent,
        identity = chosen_id,
        single_score = chosen_ss,
        presence_score = chosen_ps,
        abundance_score = chosen_as,
        recall = recall_mat[did, parent_idx],
        phi01  = phi01_mat[did, parent_idx],
        ppmi_tanh = ppmi_tanh_mat[did, parent_idx],
        f1 = f1_mat[did, parent_idx],
        spearman01 = spearman01_mat[did, parent_idx],
        hellinger  = hellinger_mat[did, parent_idx],
        dom_mean_ratio  = dominance_mean_ratio_mat[did, parent_idx],
        dom_frac_exceed = dominance_frac_exceed_mat[did, parent_idx],
        dom_n_both      = dominance_n_both_mat[did, parent_idx]
      ))
    }

    # ---- Identity neighborhood + dominance, with presence fallback ----
    min_mis <- min(misvec[ok_id], na.rm = TRUE)
    if (!is.finite(min_mis)) {
      return(c(
        parent = NA_character_, identity = max(idvec, na.rm = TRUE),
        single_score = NA_real_, presence_score = NA_real_, abundance_score = NA_real_,
        recall = NA_real_, phi01 = NA_real_, ppmi_tanh = NA_real_, f1 = NA_real_,
        spearman01 = NA_real_, hellinger = NA_real_,
        dom_mean_ratio = NA_real_, dom_frac_exceed = NA_real_, dom_n_both = NA_real_
      ))
    }

    in_neigh <- !is.na(misvec) & (misvec <= (min_mis + identity_neighborhood_mismatches))

    cand_all <- which(in_neigh & ok_id & dom_ok)
    if (length(cand_all) == 0) {
      return(c(
        parent = NA_character_, identity = max(idvec, na.rm = TRUE),
        single_score = NA_real_, presence_score = NA_real_, abundance_score = NA_real_,
        recall = NA_real_, phi01 = NA_real_, ppmi_tanh = NA_real_, f1 = NA_real_,
        spearman01 = NA_real_, hellinger = NA_real_,
        dom_mean_ratio = NA_real_, dom_frac_exceed = NA_real_, dom_n_both = NA_real_
      ))
    }

    cand_with_presence <- cand_all[!is.na(pscore[cand_all])]
    cand <- if (length(cand_with_presence) > 0) cand_with_presence else cand_all

    # Split: best-identity vs. lower-identity within neighborhood
    base_cand <- cand[misvec[cand] == min_mis]
    alt_cand  <- setdiff(cand, base_cand)

    # If multiple at best identity, tie-break by single score -> totals -> alphabetical
    if (length(base_cand) > 1) {
      ss_sub <- sscore[base_cand]
      if (!all(is.na(ss_sub))) {
        base_cand <- base_cand[which(ss_sub == max(ss_sub, na.rm = TRUE))]
      }
      if (length(base_cand) > 1) {
        totals_sub <- clean_totals[base_cand]
        base_cand <- base_cand[which(totals_sub == max(totals_sub, na.rm = TRUE))]
        if (length(base_cand) > 1) {
          base_cand <- base_cand[order(names(idvec)[base_cand])][1]
        }
      }
    } else {
      base_cand <- base_cand[1]
    }

    ss_cand <- sscore[cand]
    choose_idx <- base_cand

    # z-margin override (only if we have single scores)
    if (!all(is.na(ss_cand))) {
      z_neigh <- .row_center_scale(ss_cand, method = "z")
      names(z_neigh) <- names(idvec)[cand]
      z_base <- z_neigh[names(idvec)[base_cand]]

      if (length(alt_cand) > 0 && isTRUE(is.finite(z_base))) {
        z_alt <- z_neigh[names(idvec)[alt_cand]]
        best_alt_local <- alt_cand[which.max(z_alt)]

        if (is.finite(z_neigh[names(idvec)[best_alt_local]]) &&
            (z_neigh[names(idvec)[best_alt_local]] - z_base) >= neighborhood_z_margin) {
          choose_idx <- best_alt_local
        }
      }
    }

    parent_idx <- choose_idx
    parent <- names(idvec)[parent_idx]

    chosen_id <- idvec[parent_idx]
    chosen_ps <- pscore[parent_idx]
    chosen_as <- ascore[parent_idx]
    chosen_ss <- sscore[parent_idx]

    # Final thresholds (presence)
    if (!is.na(chosen_ps) && chosen_ps < min_presence_score) {
      return(c(
        parent = NA_character_, identity = chosen_id,
        single_score = chosen_ss, presence_score = chosen_ps, abundance_score = chosen_as,
        recall = NA_real_, phi01 = NA_real_, ppmi_tanh = NA_real_, f1 = NA_real_,
        spearman01 = NA_real_, hellinger = NA_real_,
        dom_mean_ratio  = dominance_mean_ratio_mat[did, parent_idx],
        dom_frac_exceed = dominance_frac_exceed_mat[did, parent_idx],
        dom_n_both      = dominance_n_both_mat[did, parent_idx]
      ))
    }

    c(
      parent = parent,
      identity = chosen_id,
      single_score = chosen_ss,
      presence_score = chosen_ps,
      abundance_score = chosen_as,
      recall = recall_mat[did, parent_idx],
      phi01  = phi01_mat[did, parent_idx],
      ppmi_tanh = ppmi_tanh_mat[did, parent_idx],
      f1 = f1_mat[did, parent_idx],
      spearman01 = spearman01_mat[did, parent_idx],
      hellinger  = hellinger_mat[did, parent_idx],
      dom_mean_ratio  = dominance_mean_ratio_mat[did, parent_idx],
      dom_frac_exceed = dominance_frac_exceed_mat[did, parent_idx],
      dom_n_both      = dominance_n_both_mat[did, parent_idx]
    )
  }

  # ---- Map each DIRTY to a CLEAN parent (robust schema) ----
  expected_cols <- c(
    "parent","identity","single_score","presence_score","abundance_score",
    "recall","phi01","ppmi_tanh","f1","spearman01","hellinger",
    "dom_mean_ratio","dom_frac_exceed","dom_n_both"
  )

  rowdfs <- lapply(dirty_ids, function(did) {
    v <- pick_parent(did)
    df <- as.data.frame(as.list(v), stringsAsFactors = FALSE)

    missing <- setdiff(expected_cols, names(df))
    for (nm in missing) df[[nm]] <- if (nm == "parent") NA_character_ else NA_real_

    df[, expected_cols, drop = FALSE]
  })

  mapping <- do.call(rbind, rowdfs)
  rownames(mapping) <- dirty_ids

  mapping$dirty_id <- dirty_ids
  mapping$parent   <- as.character(mapping$parent)

  numeric_cols <- setdiff(expected_cols, "parent")
  for (nm in numeric_cols) mapping[[nm]] <- suppressWarnings(as.numeric(mapping[[nm]]))

  mapping <- mapping[, c("dirty_id","parent", numeric_cols), drop = FALSE]

  # ---- Refill mapping numerics from the source matrices (authoritative) ----
  numeric_cols <- c(
    "identity","single_score","presence_score","abundance_score",
    "recall","phi01","ppmi_tanh","f1","spearman01","hellinger",
    "dom_mean_ratio","dom_frac_exceed","dom_n_both"
  )
  for (nm in numeric_cols) if (!nm %in% colnames(mapping)) mapping[[nm]] <- NA_real_

  idx <- !is.na(mapping$parent)
  rid <- mapping$dirty_id[idx]
  cid <- mapping$parent[idx]

  get_cell <- function(M, r, c) M[r, c]

  mapping$identity[idx]        <- mapply(get_cell, MoreArgs = list(M = identity_mat), r = rid, c = cid)
  mapping$single_score[idx]    <- mapply(get_cell, MoreArgs = list(M = single_score_mat), r = rid, c = cid)
  mapping$presence_score[idx]  <- mapply(get_cell, MoreArgs = list(M = presence_score_mat), r = rid, c = cid)
  mapping$abundance_score[idx] <- mapply(get_cell, MoreArgs = list(M = abundance_score_mat), r = rid, c = cid)

  mapping$recall[idx]          <- mapply(get_cell, MoreArgs = list(M = recall_mat), r = rid, c = cid)
  mapping$phi01[idx]           <- mapply(get_cell, MoreArgs = list(M = phi01_mat), r = rid, c = cid)
  mapping$ppmi_tanh[idx]       <- mapply(get_cell, MoreArgs = list(M = ppmi_tanh_mat), r = rid, c = cid)
  mapping$f1[idx]              <- mapply(get_cell, MoreArgs = list(M = f1_mat), r = rid, c = cid)
  mapping$spearman01[idx]      <- mapply(get_cell, MoreArgs = list(M = spearman01_mat), r = rid, c = cid)
  mapping$hellinger[idx]       <- mapply(get_cell, MoreArgs = list(M = hellinger_mat), r = rid, c = cid)

  mapping$dom_mean_ratio[idx]  <- mapply(get_cell, MoreArgs = list(M = dominance_mean_ratio_mat), r = rid, c = cid)
  mapping$dom_frac_exceed[idx] <- mapply(get_cell, MoreArgs = list(M = dominance_frac_exceed_mat), r = rid, c = cid)
  mapping$dom_n_both[idx]      <- mapply(get_cell, MoreArgs = list(M = dominance_n_both_mat), r = rid, c = cid)

  for (nm in numeric_cols) mapping[[nm]] <- suppressWarnings(as.numeric(mapping[[nm]]))

  list(mapping = mapping)
}
