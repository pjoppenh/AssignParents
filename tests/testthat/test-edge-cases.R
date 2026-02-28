test_that("perfect identity collapses correctly", {

  clean_seqs <- c(C1 = "ACGTACGTACGT")
  dirty_seqs <- c(D1 = "ACGTACGTACGT")

  clean_tab  <- matrix(1, nrow = 1,
                       dimnames = list("C1", "S1"))
  dirty_tab  <- matrix(2, nrow = 1,
                       dimnames = list("D1", "S1"))

  out <- assign_parents(
    clean_tab, dirty_tab,
    clean_seqs, dirty_seqs,
    min_identity = 0.95,
    presence_cutoff = 1,
    min_n = 1,
    keep_unassigned = TRUE,
    n_cores = 1,
    parallel_method = "none"
  )

  # Count conservation
  expect_equal(sum(clean_tab) + sum(dirty_tab),
               sum(out$collapsed_table) +
                 if (!is.null(out$unassigned_dirty_table))
                   sum(out$unassigned_dirty_table) else 0)

  # Perfect match should map
  expect_equal(nrow(out$mapping), 1)
})


test_that("no identity match leaves DIRTY unassigned", {

  clean_seqs <- c(C1 = "ACGTACGTACGT")
  dirty_seqs <- c(D1 = "TTTTTTTTTTTT")

  clean_tab  <- matrix(1, nrow = 1,
                       dimnames = list("C1", "S1"))
  dirty_tab  <- matrix(2, nrow = 1,
                       dimnames = list("D1", "S1"))

  out <- assign_parents(
    clean_tab, dirty_tab,
    clean_seqs, dirty_seqs,
    min_identity = 0.95,
    presence_cutoff = 1,
    min_n = 1,
    keep_unassigned = TRUE,
    n_cores = 1,
    parallel_method = "none"
  )

  # DIRTY should remain unassigned
  expect_true(!is.null(out$unassigned_dirty_table))
  expect_equal(sum(out$unassigned_dirty_table), 2)
})
