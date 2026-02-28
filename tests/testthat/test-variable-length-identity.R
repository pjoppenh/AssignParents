test_that("variable-length sequences use Levenshtein identity correctly", {

  clean_seqs <- c(
    C1 = "ACGTACGTACGT",
    C2 = "ACGTACGTACGTA"   # one base longer
  )

  dirty_seqs <- c(
    D1 = "ACGTACGTACGTA",  # perfect match to C2
    D2 = "ACGTACGTACGTT",  # one mismatch from C2
    D3 = "TTTTTTTTTTTTT"   # no match
  )

  clean_tab <- matrix(
    c(5, 1,
      2, 3),
    nrow = 2,
    dimnames = list(names(clean_seqs), c("S1","S2"))
  )

  dirty_tab <- matrix(
    c(1,0,
      0,1,
      3,3),
    nrow = 3,
    dimnames = list(names(dirty_seqs), c("S1","S2"))
  )

  out <- assign_parents(
    clean_tab, dirty_tab,
    clean_seqs, dirty_seqs,
    min_identity = 0.80,
    presence_cutoff = 1,
    min_n = 1,
    keep_unassigned = TRUE,
    n_cores = 1,
    parallel_method = "none"
  )

  # Structural sanity
  expect_true(is.matrix(out$collapsed_table))
  expect_true(is.data.frame(out$mapping))

  # Count conservation
  total_in  <- sum(clean_tab) + sum(dirty_tab)
  total_out <- sum(out$collapsed_table) +
    if (!is.null(out$unassigned_dirty_table))
      sum(out$unassigned_dirty_table) else 0

  expect_equal(total_in, total_out)
})
