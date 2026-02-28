#' Assign DIRTY ASVs to CLEAN parents
#'
#' Collapses DIRTY ASVs into CLEAN ASVs using sequence identity,
#' prevalence-aware presence similarity, abundance similarity,
#' dominance filters, and a standardized meta-score.
#'
#' @param clean_tab Numeric matrix of CLEAN ASVs (rows = ASVs, cols = samples)
#' @param dirty_tab Numeric matrix of DIRTY ASVs (rows = ASVs, cols = samples)
#' @param clean_seqs Named character vector of CLEAN sequences
#' @param dirty_seqs Named character vector of DIRTY sequences
#'
#' @param min_identity Minimum percent identity required
#' @param min_presence_score Minimum presence composite threshold
#' @param presence_cutoff Abundance threshold for presence/absence binarization
#' @param min_recall Minimum recall required in presence composite
#' @param min_n Minimum co-occurrence sample size
#' @param keep_unassigned Logical; keep DIRTY ASVs with no parent
#'
#' @param presence_weights Named numeric vector for presence composite
#' @param abundance_weights Named numeric vector for abundance composite
#' @param single_weights Named numeric vector for final meta-score
#' @param single_scale_method Scaling method for row-wise normalization
#'
#' @param enforce_dominance Logical; enforce dominance filtering
#' @param dominance_mean_ratio_max Maximum allowed normalized mean ratio
#' @param dominance_frac_exceed_max Maximum allowed exceedance fraction
#' @param dominance_min_both Minimum shared samples for dominance test
#'
#' @param identity_neighborhood_mismatches Neighborhood size in mismatches
#' @param neighborhood_z_margin Z-score margin required to override identity
#'
#' @param n_cores Number of cores for parallelization
#' @param parallel_method Parallel backend. "auto" chooses "fork" on macOS/Linux
#'   and "cluster" on Windows when n_cores > 1; otherwise "none".
#'
#' @return A list containing collapsed tables, mapping, score matrices,
#'         dominance diagnostics, and matched/unmatched sequence vectors.
#'
#' @export
assign_parents <- function(clean_tab, dirty_tab, clean_seqs, dirty_seqs,
                           min_identity = 0.98,
                           min_presence_score = 0,
                           presence_cutoff = 1,
                           min_recall = 0.5,
                           min_n = 6,
                           keep_unassigned = TRUE,
                           presence_weights = c(phi = 0.5, ppmi = 0.3, f1 = 0.2),
                           abundance_weights = c(spearman = 0.6, hellinger = 0.4),
                           single_weights = c(presence = 0.6, abundance = 0.4),
                           single_scale_method = c("z", "robust", "rank"),
                           enforce_dominance = TRUE,
                           dominance_mean_ratio_max  = 0.9,
                           dominance_frac_exceed_max = 0.35,
                           dominance_min_both = 3,
                           identity_neighborhood_mismatches = 1L,
                           neighborhood_z_margin = 1.4,
                           n_cores = 1L,
                           parallel_method = c("auto", "fork", "cluster", "none")) {

  if (!requireNamespace("stringdist", quietly = TRUE)) {
    stop("Package 'stringdist' is required but not installed.")
  }

  .notify <- function(msg) {
    message(sprintf("[%s] %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), msg))
  }

  single_scale_method <- match.arg(single_scale_method)
  parallel_method <- match.arg(parallel_method)

  # ---- Normalize n_cores + parallel_method ----
  n_cores <- as.integer(n_cores)
  if (is.na(n_cores) || n_cores <= 1L) n_cores <- 1L

  is_windows <- identical(.Platform$OS.type, "windows")

  if (parallel_method == "auto") {
    if (n_cores <= 1L) {
      parallel_method <- "none"
    } else {
      parallel_method <- if (is_windows) "cluster" else "fork"
    }
  }

  # If user explicitly requests fork on Windows, fall back
  if (parallel_method == "fork" && is_windows) {
    parallel_method <- "cluster"
  }

  # If parallel is requested but n_cores is 1, force none
  if (n_cores <= 1L) {
    parallel_method <- "none"
  }

  # ---- Input checks + harmonize ----
  x <- validate_inputs(
    clean_tab  = clean_tab,
    dirty_tab  = dirty_tab,
    clean_seqs = clean_seqs,
    dirty_seqs = dirty_seqs
  )
  clean_tab  <- x$clean_tab
  dirty_tab  <- x$dirty_tab
  clean_seqs <- x$clean_seqs
  dirty_seqs <- x$dirty_seqs

  # ---- Parallel setup ----
  # setup_parallel() should interpret:
  # - parallel_method == "none"            => sequential (par_lapply = lapply; stop = no-op)
  # - parallel_method == "fork"            => mclapply (unix only)
  # - parallel_method == "cluster"         => PSOCK cluster (Windows / explicit)
  par <- setup_parallel(
    n_cores = n_cores,
    parallel_method = parallel_method,
    export = c(
      "presence_composite",
      "abundance_composite"
    )
  )
  on.exit(try(par$stop(), silent = TRUE), add = TRUE)

  # ---- Identity + mismatches (DIRTY x CLEAN) ----
  .notify("Starting identity matrix computation...")
  id_out <- compute_identity_matrix(
    clean_tab     = clean_tab,
    dirty_tab     = dirty_tab,
    clean_seqs    = clean_seqs,
    dirty_seqs    = dirty_seqs,
    min_identity  = min_identity,
    par           = par
  )
  identity_mat   <- id_out$identity
  mismatches_mat <- id_out$mismatches
  .notify("Identity/mismatch matrices done.")

  # ---- Score matrices (presence/abundance + dominance + components) ----
  .notify("Starting DIRTYxCLEAN scoring (presence/abundance/dominance)...")
  score_out <- compute_score_matrices(
    clean_tab        = clean_tab,
    dirty_tab        = dirty_tab,
    identity_mat     = identity_mat,
    min_identity     = min_identity,
    presence_cutoff  = presence_cutoff,
    presence_weights = presence_weights,
    abundance_weights = abundance_weights,
    min_recall       = min_recall,
    min_n            = min_n,
    enforce_dominance = enforce_dominance,
    dominance_mean_ratio_max  = dominance_mean_ratio_max,
    dominance_frac_exceed_max = dominance_frac_exceed_max,
    dominance_min_both = dominance_min_both,
    par = par
  )
  .notify("Scoring complete. Building single_score matrix...")

  single_score_mat <- compute_single_scores(
    presence_mat  = score_out$presence_score,
    abundance_mat = score_out$abundance_score,
    weights       = single_weights,
    scale_method  = single_scale_method
  )

  # ---- Parent selection ----
  .notify("Selecting parents / building mapping...")
  map_out <- select_parents(
    identity_mat       = identity_mat,
    mismatches_mat     = mismatches_mat,
    single_score_mat   = single_score_mat,
    presence_score_mat = score_out$presence_score,
    abundance_score_mat = score_out$abundance_score,
    recall_mat         = score_out$recall,
    phi01_mat          = score_out$phi01,
    ppmi_tanh_mat      = score_out$ppmi_tanh,
    f1_mat             = score_out$f1,
    spearman01_mat     = score_out$spearman01,
    hellinger_mat      = score_out$hellinger,
    dominance_mean_ratio_mat = score_out$dominance_mean_ratio,
    dominance_frac_exceed_mat = score_out$dominance_frac_exceed,
    dominance_n_both_mat = score_out$dominance_n_both,
    dominance_pass_mat  = score_out$dominance_pass,
    min_identity        = min_identity,
    min_presence_score  = min_presence_score,
    enforce_dominance   = enforce_dominance,
    identity_neighborhood_mismatches = identity_neighborhood_mismatches,
    neighborhood_z_margin = neighborhood_z_margin,
    clean_tab = clean_tab
  )
  mapping <- map_out$mapping

  # ---- Collapse counts & split unmatched ----
  .notify("Collapsing tables...")
  col_out <- collapse_tables(
    clean_tab        = clean_tab,
    dirty_tab        = dirty_tab,
    mapping          = mapping,
    clean_seqs       = clean_seqs,
    dirty_seqs       = dirty_seqs,
    keep_unassigned  = keep_unassigned
  )

  # ---- Final output ----
  list(
    collapsed_table        = col_out$collapsed_table,
    unmatched_clean_table  = col_out$unmatched_clean_table,
    unassigned_dirty_table = col_out$unassigned_dirty_table,
    mapping                = mapping,
    identity               = identity_mat,
    presence_score         = score_out$presence_score,
    abundance_score        = score_out$abundance_score,
    single_score           = single_score_mat,
    dominance_mean_ratio   = score_out$dominance_mean_ratio,
    dominance_frac_exceed  = score_out$dominance_frac_exceed,
    dominance_n_both       = score_out$dominance_n_both,
    dominance_pass         = score_out$dominance_pass,
    matched_clean_ids      = col_out$matched_clean_ids,
    unmatched_clean_ids    = col_out$unmatched_clean_ids,
    matched_clean_seqs     = col_out$matched_clean_seqs,
    unmatched_clean_seqs   = col_out$unmatched_clean_seqs,
    unmatched_dirty_seqs   = col_out$unmatched_dirty_seqs
  )
}
