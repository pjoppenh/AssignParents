test_that("parallel and sequential produce identical results", {

  case <- make_random_case(
    seed = 123,
    n_clean = 50,
    n_dirty = 100,
    n_samp = 5,
    L = 80
  )

  out_seq <- assign_parents(
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

  out_cluster <- assign_parents(
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
  )

  # On mac/linux you can also test fork
  if (.Platform$OS.type != "windows") {
    out_fork <- assign_parents(
      case$clean_tab,
      case$dirty_tab,
      case$clean_seqs,
      case$dirty_seqs,
      min_identity = 0.95,
      presence_cutoff = 1,
      min_n = 2,
      keep_unassigned = TRUE,
      n_cores = 2,
      parallel_method = "fork"
    )

    expect_equal(out_seq$collapsed_table, out_fork$collapsed_table)
    expect_equal(out_seq$mapping, out_fork$mapping)
  }

  expect_equal(out_seq$collapsed_table, out_cluster$collapsed_table)
  expect_equal(out_seq$mapping, out_cluster$mapping)
})
