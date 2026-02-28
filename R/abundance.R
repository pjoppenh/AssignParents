#' Abundance composite score between CLEAN and DIRTY
#'
#' @keywords internal
abundance_composite <- function(p_counts,
                                d_counts,
                                weights = c(spearman = 0.6,
                                            hellinger = 0.4)) {

  eps <- .Machine$double.eps

  # ---- Spearman correlation ----
  sp <- suppressWarnings(
    stats::cor(p_counts, d_counts,
               method = "spearman",
               use = "pairwise.complete.obs")
  )

  if (is.na(sp)) sp <- 0
  spearman01 <- (sp + 1) / 2  # rescale to 0–1

  # ---- Hellinger similarity ----
  p_rel <- p_counts / max(sum(p_counts), eps)
  d_rel <- d_counts / max(sum(d_counts), eps)

  hel_dist <- sqrt(sum((sqrt(p_rel) - sqrt(d_rel))^2))
  hellinger <- 1 - (hel_dist / sqrt(2))  # similarity form (0–1)

  if (!is.finite(hellinger)) hellinger <- 0

  # ---- Weighted composite ----
  w <- weights / sum(weights)

  score <- w["spearman"] * spearman01 +
    w["hellinger"] * hellinger

  list(
    score = score,
    components = list(
      spearman = spearman01,
      hellinger = hellinger
    )
  )
}
