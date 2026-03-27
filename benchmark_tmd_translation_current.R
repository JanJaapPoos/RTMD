benchmark_tmd_translation <- function(
  filename = "Subset_NS_EC_forTMD30_21_2_2025.nc",
  date = as.Date("2025-01-01"),
  n_iter = 3
) {
  stopifnot(n_iter >= 1)

  root <- "M:/TMD3.0 translation"
  loader_r <- file.path(root, "load_tmd_translation_current_pure.R")
  loader_rcpp <- file.path(root, "load_tmd_translation_current.R")

  if (!file.exists(loader_r)) stop("Missing pure-R loader: ", loader_r, call. = FALSE)
  if (!file.exists(loader_rcpp)) stop("Missing current loader: ", loader_rcpp, call. = FALSE)
  if (!file.exists(file.path(root, filename))) stop("Missing model file: ", file.path(root, filename), call. = FALSE)

  paired_lat <- seq(52, 55, 0.1)
  paired_lon <- rep(2.95, length(paired_lat))
  grid_lat <- seq(52, 55, 0.25)
  grid_lon <- seq(2, 4, 0.25)

  cases <- list(
    list(name = "single_point", args = list(filename, 52.5, 2.95, date, input_mode = "paired")),
    list(name = "paired_points", args = list(filename, paired_lat, paired_lon, date, input_mode = "paired")),
    list(name = "grid_points", args = list(filename, grid_lat, grid_lon, date, input_mode = "grid"))
  )

  time_call <- function(expr, n) {
    times <- numeric(n)
    value <- NULL
    for (i in seq_len(n)) {
      gc()
      t0 <- proc.time()[["elapsed"]]
      value <- force(expr)
      t1 <- proc.time()[["elapsed"]]
      times[i] <- t1 - t0
    }
    list(times = times, value = value)
  }

  summarize_times <- function(x) {
    c(min = min(x), median = stats::median(x), mean = mean(x), max = max(x))
  }

  compare_values <- function(a, b) {
    if (!identical(dim(a), dim(b))) {
      return(list(equal = FALSE, reason = "different dimensions", max_abs_diff = NA_real_))
    }
    a_num <- as.vector(a)
    b_num <- as.vector(b)
    if (length(a_num) != length(b_num)) {
      return(list(equal = FALSE, reason = "different lengths", max_abs_diff = NA_real_))
    }
    diff <- suppressWarnings(max(Mod(a_num - b_num), na.rm = TRUE))
    if (!is.finite(diff)) diff <- 0
    list(equal = isTRUE(all.equal(a, b, tolerance = 1e-10)), reason = NULL, max_abs_diff = diff)
  }

  run_suite <- function(loader) {
    env <- new.env(parent = baseenv())
    env$source <- function(file, ...) sys.source(file, envir = env)
    sys.source(loader, envir = env)
    if (!exists("tmd_predict", envir = env, inherits = FALSE) || !is.function(env$tmd_predict)) {
      stop("Loader did not define tmd_predict in the benchmark environment.", call. = FALSE)
    }
    lapply(cases, function(case) {
      res <- time_call(do.call(env$tmd_predict, case$args), n_iter)
      list(name = case$name, times = res$times, value = res$value)
    })
  }

  cat("Benchmarking pure-R loader...\n")
  res_r <- run_suite(loader_r)

  cat("Benchmarking current loader...\n")
  res_rcpp <- run_suite(loader_rcpp)

  rows <- vector("list", length(cases))
  checks <- vector("list", length(cases))

  for (i in seq_along(cases)) {
    s_r <- summarize_times(res_r[[i]]$times)
    s_cpp <- summarize_times(res_rcpp[[i]]$times)
    cmp <- compare_values(res_r[[i]]$value, res_rcpp[[i]]$value)

    denom <- unname(s_cpp["median"])
    speedup <- if (denom == 0) NA_real_ else unname(s_r["median"] / denom)

    rows[[i]] <- data.frame(
      case = cases[[i]]$name,
      r_min = s_r["min"],
      r_median = s_r["median"],
      r_mean = s_r["mean"],
      current_min = s_cpp["min"],
      current_median = s_cpp["median"],
      current_mean = s_cpp["mean"],
      speedup_median = speedup,
      stringsAsFactors = FALSE
    )

    checks[[i]] <- data.frame(
      case = cases[[i]]$name,
      equal = cmp$equal,
      max_abs_diff = cmp$max_abs_diff,
      reason = if (is.null(cmp$reason)) "" else cmp$reason,
      stringsAsFactors = FALSE
    )
  }

  timing_table <- do.call(rbind, rows)
  check_table <- do.call(rbind, checks)

  cat("\nTiming summary (seconds):\n")
  print(timing_table, row.names = FALSE)

  cat("\nOutput comparison:\n")
  print(check_table, row.names = FALSE)

  invisible(list(
    timing = timing_table,
    checks = check_table,
    raw = list(r = res_r, current = res_rcpp)
  ))
}

if (sys.nframe() == 0) {
  benchmark_tmd_translation()
}
