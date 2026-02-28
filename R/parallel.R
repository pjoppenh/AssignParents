#' @keywords internal
setup_parallel <- function(n_cores = 1L,
                           parallel_method = c("auto", "fork", "cluster", "none"),
                           export = NULL) {

  parallel_method <- match.arg(parallel_method)

  n_cores <- as.integer(n_cores)
  if (is.na(n_cores) || n_cores <= 1L) n_cores <- 1L

  if (parallel_method == "auto") {
    parallel_method <- if (n_cores > 1L) "cluster" else "none"
  }

  # Sequential
  if (parallel_method == "none" || n_cores <= 1L) {
    return(list(
      par_lapply = function(X, FUN, ...) base::lapply(X, FUN, ...),
      stop = function() invisible(NULL)
    ))
  }

  # Fork (optional; not Windows)
  if (parallel_method == "fork") {
    if (identical(.Platform$OS.type, "windows")) {
      parallel_method <- "cluster"
    } else {
      return(list(
        par_lapply = function(X, FUN, ...) parallel::mclapply(X, FUN, mc.cores = n_cores, ...),
        stop = function() invisible(NULL)
      ))
    }
  }

  # Cluster (PSOCK)
  if (parallel_method == "cluster") {
    cl <- parallel::makeCluster(n_cores)

    # Export functions/objects so workers can run code during devtools::load_all()
    # 'export' should be a character vector of function names.
    if (!is.null(export) && length(export) > 0) {
      parallel::clusterExport(cl, varlist = export, envir = parent.frame())
    }

    return(list(
      par_lapply = function(X, FUN, ...) parallel::parLapplyLB(cl, X, FUN, ...),
      stop = function() parallel::stopCluster(cl)
    ))
  }

  stop("setup_parallel(): unsupported parallel_method = ", parallel_method)
}
