#' Presence composite score between CLEAN and DIRTY
#'
#' @keywords internal
presence_composite <- function(p_counts,
                               d_counts,
                               presence_cutoff = 1,
                               weights = c(phi = 0.5, ppmi = 0.3, f1 = 0.2),
                               min_recall = 0.5,
                               min_n = 6) {

  # ---- Binarize ----
  p_bin <- as.integer(p_counts > presence_cutoff)
  d_bin <- as.integer(d_counts > presence_cutoff)

  both_present <- (p_bin == 1) & (d_bin == 1)
  n_both <- sum(both_present, na.rm = TRUE)

  if (n_both < min_n) {
    return(list(
      score = NA_real_,
      components = list(
        recall = NA_real_,
        phi01 = NA_real_,
        ppmi_tanh = NA_real_,
        f1 = NA_real_
      ),
      passed_gate = FALSE
    ))
  }

  # ---- Confusion matrix ----
  tp <- sum(p_bin == 1 & d_bin == 1)
  fp <- sum(p_bin == 1 & d_bin == 0)
  fn <- sum(p_bin == 0 & d_bin == 1)
  tn <- sum(p_bin == 0 & d_bin == 0)

  # ---- Recall ----
  recall <- if ((tp + fn) > 0) tp / (tp + fn) else NA_real_

  if (is.na(recall) || recall < min_recall) {
    return(list(
      score = NA_real_,
      components = list(
        recall = recall,
        phi01 = NA_real_,
        ppmi_tanh = NA_real_,
        f1 = NA_real_
      ),
      passed_gate = FALSE
    ))
  }

  # ---- Phi coefficient ----
  denom <- sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn))
  phi <- if (denom > 0) (tp*tn - fp*fn) / denom else 0
  phi01 <- (phi + 1) / 2  # rescale to 0–1

  # ---- PPMI ----
  n <- tp + fp + fn + tn
  p_xy <- tp / n
  p_x  <- (tp + fp) / n
  p_y  <- (tp + fn) / n

  ppmi <- if (p_xy > 0 && p_x > 0 && p_y > 0) {
    max(log2(p_xy / (p_x * p_y)), 0)
  } else {
    0
  }

  ppmi_tanh <- tanh(ppmi)

  # ---- F1 ----
  precision <- if ((tp + fp) > 0) tp / (tp + fp) else 0
  f1 <- if ((precision + recall) > 0) {
    2 * precision * recall / (precision + recall)
  } else {
    0
  }

  # ---- Weighted composite ----
  w <- weights / sum(weights)

  score <- w["phi"] * phi01 +
    w["ppmi"] * ppmi_tanh +
    w["f1"] * f1

  list(
    score = score,
    components = list(
      recall = recall,
      phi01 = phi01,
      ppmi_tanh = ppmi_tanh,
      f1 = f1
    ),
    passed_gate = TRUE
  )
}
