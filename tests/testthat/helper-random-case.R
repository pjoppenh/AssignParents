make_random_case <- function(seed,
                             n_clean = 50,
                             n_dirty = 100,
                             n_samp = 5,
                             L = 80) {

  set.seed(seed)

  bases <- c("A","C","G","T")

  rand_seq <- function(n, L) {
    replicate(n, paste0(sample(bases, L, replace = TRUE), collapse = ""))
  }

  clean_seqs <- rand_seq(n_clean, L)
  dirty_seqs <- rand_seq(n_dirty, L)

  names(clean_seqs) <- paste0("C", seq_len(n_clean))
  names(dirty_seqs) <- paste0("D", seq_len(n_dirty))

  clean_tab <- matrix(
    rpois(n_clean * n_samp, lambda = 5),
    nrow = n_clean,
    dimnames = list(names(clean_seqs),
                    paste0("S", seq_len(n_samp)))
  )

  dirty_tab <- matrix(
    rpois(n_dirty * n_samp, lambda = 3),
    nrow = n_dirty,
    dimnames = list(names(dirty_seqs),
                    paste0("S", seq_len(n_samp)))
  )

  list(
    clean_tab = clean_tab,
    dirty_tab = dirty_tab,
    clean_seqs = clean_seqs,
    dirty_seqs = dirty_seqs
  )
}
