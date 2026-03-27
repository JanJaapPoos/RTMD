# Optional Rcpp acceleration for the hottest numerical kernels.

.tmd_rcpp_env <- new.env(parent = emptyenv())
.tmd_rcpp_env$attempted <- FALSE
.tmd_rcpp_env$available <- FALSE

.tmd_init_rcpp <- function() {
  if (.tmd_rcpp_env$attempted) {
    return(.tmd_rcpp_env$available)
  }
  .tmd_rcpp_env$attempted <- TRUE

  if (!requireNamespace("Rcpp", quietly = TRUE)) {
    return(FALSE)
  }

  cpp_path <- file.path(get0(".tmd_translation_dir", ifnotfound = getwd()), "tmd_rcpp_backend.cpp")
  if (!file.exists(cpp_path)) {
    return(FALSE)
  }

  ok <- tryCatch({
    Rcpp::sourceCpp(cpp_path, rebuild = FALSE, verbose = FALSE)
    exists("interp_linear_multi_cpp", mode = "function", inherits = TRUE) &&
      exists("tmd_harp_cpp", mode = "function", inherits = TRUE)
  }, error = function(e) FALSE)

  .tmd_rcpp_env$available <- isTRUE(ok)
  .tmd_rcpp_env$available
}

.interp_linear_multi_r <- .interp_linear_multi
tmd_harp_r <- tmd_harp

.interp_linear_multi <- function(x, y, z, xi, yi, extrap = NA_real_) {
  if (!.tmd_init_rcpp()) {
    return(.interp_linear_multi_r(x, y, z, xi, yi, extrap = extrap))
  }

  g <- .prepare_interp_grid(x, y, z)
  z3 <- .as_3d_array(g$z)
  shape <- .query_shape(xi)
  xi_v <- shape$vector
  yi_v <- as.vector(yi)

  out <- interp_linear_multi_cpp(
    x = g$x,
    y = g$y,
    z = as.vector(z3),
    zdim = as.integer(dim(z3)),
    xi = xi_v,
    yi = yi_v,
    extrap = extrap
  )

  array(out, dim = c(shape$dim, dim(z3)[3]))
}

tmd_harp <- function(t, hc, constituents, p, N, ph, omega) {
  if (!.tmd_init_rcpp()) {
    return(tmd_harp_r(t, hc, constituents, p, N, ph, omega))
  }

  t <- .matlab_datenum(t)
  constituents <- as.character(constituents)
  hc <- as.matrix(hc)
  if (nrow(hc) != length(constituents)) {
    stop("Rows of hc must match the number of constituents.", call. = FALSE)
  }

  time_s <- (t - 727564) * 86400
  nodal <- tmd_nodal(constituents, p, N)

  tmd_harp_cpp(
    time_s = as.numeric(time_s),
    hc_re = unname(Re(hc)),
    hc_im = unname(Im(hc)),
    pu = unname(nodal$pu),
    pf = unname(nodal$pf),
    ph = as.numeric(ph),
    omega = as.numeric(omega)
  )
}
