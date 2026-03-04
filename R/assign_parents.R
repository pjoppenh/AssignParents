#' Assign DIRTY ASVs to CLEAN parent ASVs
#'
#' Post-denoising ASV refinement: collapse likely artifact ASVs (DIRTY) into
#' parent ASVs (CLEAN) using sequence identity plus co-occurrence-based scoring,
#' with optional dominance filtering. Counts from assigned DIRTY ASVs are added
#' to the selected CLEAN parent.
#'
#' @param clean_tab Numeric matrix of CLEAN ASVs (rows = ASVs, cols = samples).
#'   Row names should be CLEAN ASV IDs.
#' @param dirty_tab Numeric matrix of DIRTY ASVs (rows = ASVs, cols = samples).
#'   Row names should be DIRTY ASV IDs.
#' @param clean_seqs Named character vector of CLEAN sequences. Names should match
#'   \code{rownames(clean_tab)}.
#' @param dirty_seqs Named character vector of DIRTY sequences. Names should match
#'   \code{rownames(dirty_tab)}.
#'
#' @param min_identity Numeric scalar in \code{[0,1]} giving the minimum sequence
#'   identity for a DIRTY->CLEAN candidate edge (e.g., \code{0.98} = 98% identity).
#'   Higher values reduce false merges but may miss true parents; lower values
#'   increase recall but risk incorrect merges. Typical values: 0.97–1.00.
#'
#' @param min_presence_score Numeric scalar threshold for accepting an assignment
#'   based on the presence composite score. If the presence composite is normalized,
#'   this is typically in \code{[0,1]}. Increase to require stronger co-occurrence
#'   evidence; decrease to allow weaker evidence.
#'
#' @param presence_cutoff Numeric scalar >= 0 (integer recommended for count tables).
#'   Values \eqn{>} \code{presence_cutoff} are treated as "present" for binarization
#'   used in co-occurrence metrics (e.g., Jaccard, recall). Common choices are
#'   \code{0} (any non-zero count is present) or \code{1} (require at least 2 reads).
#'
#' @param min_recall Numeric scalar in \code{[0,1]}. Minimum recall required in the
#'   presence composite filtering. Higher values require the DIRTY ASV's occurrences
#'   to be largely a subset of the parent's occurrences (direction depends on the
#'   exact implementation of recall in this package).
#'
#' @param min_n Integer >= 1. Minimum number of co-occurrence samples required for
#'   computing/accepting co-occurrence-based evidence (e.g., to avoid unstable scores
#'   from very small overlap).
#'
#' @param keep_unassigned Logical. If TRUE, DIRTY ASVs that do not receive an
#'   assignment (no candidate passes thresholds/filters) are retained as separate
#'   rows in the output; if FALSE they are dropped.
#'
#' @param presence_weights Named numeric vector with non-negative entries. Weights
#'   used to combine individual components of the presence composite score. Names
#'   must match the internal presence-component names computed by the package.
#'   Components not listed are treated as weight 0.
#'
#' @param abundance_weights Named numeric vector with non-negative entries. Weights
#'   used to combine individual components of the abundance composite score. Names
#'   must match the internal abundance-component names computed by the package.
#'   Components not listed are treated as weight 0.
#'
#' @param single_weights Named numeric vector with non-negative entries. Weights used
#'   to combine higher-level scores (typically aggregates of identity/presence/abundance)
#'   into a single final meta-score for ranking candidate parents. Names must match
#'   the internal final-score component names. Components not listed are treated as 0.
#'
#' @param single_scale_method Character string controlling row-wise normalization
#'   used when combining components into the final meta-score. One of:
#'   \code{"z"} (z-score), \code{"robust"} (robust z-score / median-MAD style),
#'   or \code{"rank"} (rank-based scaling).
#'
#' @param enforce_dominance Logical. If TRUE, apply dominance filtering to reject
#'   assignments where the DIRTY ASV is not consistently lower than the proposed CLEAN
#'   parent across shared samples (according to dominance parameters below).
#'
#' @param dominance_mean_ratio_max Numeric scalar in \code{(0, Inf)}. Threshold on
#'   the (normalized) mean abundance ratio used in dominance filtering when
#'   \code{enforce_dominance = TRUE}. Values < 1 enforce that DIRTY is lower than
#'   the proposed parent on average; smaller values make dominance stricter.
#'
#' @param dominance_frac_exceed_max Numeric scalar in \code{[0,1]}. Maximum allowed
#'   fraction of shared samples where the DIRTY ASV exceeds its proposed CLEAN parent
#'   (per the dominance comparison). \code{0} means DIRTY may never exceed the parent;
#'   \code{0.1} allows exceedance in up to 10% of shared samples; \code{1} effectively
#'   disables this constraint. Only used when \code{enforce_dominance = TRUE}.
#'
#' @param dominance_min_both Integer >= 1. Minimum number of shared samples where both
#'   DIRTY and proposed parent are present (per \code{presence_cutoff}) required to
#'   run the dominance test.
#'
#' @param identity_neighborhood_mismatches Integer >= 0. Neighborhood size (in mismatches)
#'   used by identity-neighborhood logic (if enabled) to consider near-ties among candidate
#'   parents around the best identity.
#'
#' @param neighborhood_z_margin Numeric scalar >= 0. Z-score margin required for the
#'   neighborhood logic to override a pure identity-based choice among near-tied parents.
#'   Larger values make overrides rarer/more conservative.
#'
#' @param n_cores Integer >= 1. Number of cores to use for parallelizable steps.
#'
#' @param parallel_method Character string specifying parallel backend: \code{"none"},
#'   \code{"fork"}, \code{"cluster"}, or \code{"auto"}. \code{"auto"} chooses \code{"fork"}
#'   on macOS/Linux and \code{"cluster"} on Windows when \code{n_cores > 1}; otherwise \code{"none"}.
#'   On Windows, \code{"fork"} is not available and will be treated as \code{"cluster"}.
#'
#' @return A named list including:
#' \describe{
#'   \item{collapsed_table}{CLEAN table with assigned DIRTY counts added to parents.}
#'   \item{unmatched_clean_table}{CLEAN rows that never receive DIRTY assignments (if computed).}
#'   \item{unassigned_dirty_table}{DIRTY rows that were not assigned (if \code{keep_unassigned = TRUE}).}
#'   \item{mapping}{Per-DIRTY mapping to selected parent plus diagnostics/scores.}
#'   \item{identity}{DIRTY x CLEAN identity matrix.}
#'   \item{presence_score, abundance_score, single_score}{Score matrices used for ranking/selection.}
#'   \item{dominance_*}{Dominance diagnostics matrices (mean ratio, exceedance fraction, etc.).}
#'   \item{matched_clean_ids, unmatched_clean_ids, matched_clean_seqs, unmatched_clean_seqs, unmatched_dirty_seqs}{Convenience vectors.}
#' }
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
