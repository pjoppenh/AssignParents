test_that("parallel runs are deterministic", {

  case <- make_random_case(
    seed = 456,
    n_clean = 40,
    n_dirty = 80,
    n_samp = 6,
    L = 90
  )

  run_once <- function() {
    assign_parents(
      case$clean_tab,
      case$dirty_tab,
      case$clean_seqs,
      case$dirty_seqs,
      min_identity = 0.95,
      presence_cutoff = 1,
      min_n = 2,
      keep_unassigned = TRUE,
      n_cores = 2,
      parallel_method = "cluster"
    )$mapping
  }

  m1 <- run_once()
  m2 <- run_once()
  m3 <- run_once()

  expect_equal(m1, m2)
  expect_equal(m2, m3)
})
