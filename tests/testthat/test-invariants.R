test_that("assign_parents preserves structural invariants", {

  case <- make_random_case(
    seed = 789,
    n_clean = 30,
    n_dirty = 60,
    n_samp = 4,
    L = 70
  )

  out <- assign_parents(
    case$clean_tab,
    case$dirty_tab,
    case$clean_seqs,
    case$dirty_seqs,
    min_identity = 0.95,
    presence_cutoff = 1,
    min_n = 2,
    keep_unassigned = TRUE,
    n_cores = 1,
    parallel_method = "none"
  )

  # 1) collapsed table must be matrix
  expect_true(is.matrix(out$collapsed_table))

  # 2) mapping must be data.frame
  expect_true(is.data.frame(out$mapping))

  # 3) column names preserved
  expect_identical(
    colnames(out$collapsed_table),
    colnames(case$clean_tab)
  )

  # 4) count conservation
  total_in  <- sum(case$clean_tab) + sum(case$dirty_tab)
  total_out <- sum(out$collapsed_table) +
    if (!is.null(out$unassigned_dirty_table))
      sum(out$unassigned_dirty_table) else 0

  expect_equal(total_in, total_out)

  # 5) mapped DIRTY IDs must exist
  dirty_ids <- rownames(case$dirty_tab)
  dirty_col <- intersect(
    c("dirty_id","DIRTY","dirty","dirty_asv"),
    names(out$mapping)
  )

  expect_length(dirty_col, 1)
  expect_true(all(out$mapping[[dirty_col]] %in% dirty_ids))
})
