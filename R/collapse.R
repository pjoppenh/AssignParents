#' Collapse DIRTY ASVs into CLEAN parents
#'
#' @keywords internal
collapse_tables <- function(clean_tab,
                            dirty_tab,
                            mapping,
                            clean_seqs,
                            dirty_seqs,
                            keep_unassigned = TRUE) {

  if (is.null(mapping$dirty_id) || is.null(mapping$parent)) {
    stop("mapping must contain columns 'dirty_id' and 'parent'.")
  }

  collapsed_tab <- clean_tab

  dirty_ids <- mapping$dirty_id
  parents   <- mapping$parent

  assigned_idx   <- which(!is.na(parents))
  unassigned_idx <- which(is.na(parents))

  # ---- Collapse assigned DIRTY into CLEAN ----
  for (k in assigned_idx) {

    did <- dirty_ids[k]
    pid <- parents[k]

    if (!did %in% rownames(dirty_tab)) next
    if (!pid %in% rownames(collapsed_tab)) next

    collapsed_tab[pid, ] <-
      collapsed_tab[pid, ] + dirty_tab[did, ]
  }

  # ---- Build unassigned DIRTY table ----
  unassigned_dirty_table <- NULL
  if (length(unassigned_idx) > 0 && keep_unassigned) {
    unassigned_ids <- dirty_ids[unassigned_idx]
    unassigned_dirty_table <- dirty_tab[unassigned_ids, , drop = FALSE]
  }

  # ---- Identify matched CLEAN ----
  matched_clean_ids <- unique(parents[assigned_idx])
  matched_clean_ids <- matched_clean_ids[!is.na(matched_clean_ids)]

  unmatched_clean_ids <- setdiff(rownames(clean_tab), matched_clean_ids)

  unmatched_clean_table <- NULL
  if (length(unmatched_clean_ids) > 0) {
    unmatched_clean_table <-
      clean_tab[unmatched_clean_ids, , drop = FALSE]
  }

  # ---- Sequence splits ----
  matched_clean_seqs <-
    clean_seqs[matched_clean_ids]

  unmatched_clean_seqs <-
    clean_seqs[unmatched_clean_ids]

  unmatched_dirty_seqs <- NULL
  if (length(unassigned_idx) > 0) {
    unmatched_dirty_seqs <-
      dirty_seqs[dirty_ids[unassigned_idx]]
  }

  list(
    collapsed_table        = collapsed_tab,
    unmatched_clean_table  = unmatched_clean_table,
    unassigned_dirty_table = unassigned_dirty_table,
    matched_clean_ids      = matched_clean_ids,
    unmatched_clean_ids    = unmatched_clean_ids,
    matched_clean_seqs     = matched_clean_seqs,
    unmatched_clean_seqs   = unmatched_clean_seqs,
    unmatched_dirty_seqs   = unmatched_dirty_seqs
  )
}
