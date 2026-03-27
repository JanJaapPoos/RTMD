# Fix output reshaping in tmd_predict() to match MATLAB behavior.

tmd_predict <- function(filename, lat, lon, time, ptype = "h", ..., constituents = NULL, coasts = "nan", InferMinor = NULL, input_mode = c("paired", "grid")) {
  grids <- .gridify_lat_lon(lat, lon, input_mode = input_mode)
  lat <- grids$lat
  lon <- grids$lon

  .assert_tmd_file(filename)
  if (!identical(dim(lat), dim(lon)) && length(lat) != length(lon)) {
    stop("Dimensions of lat and lon must agree.", call. = FALSE)
  }
  if (!ptype %in% c("h", "z", "u", "U", "v", "V")) {
    stop("ptype must be one of 'h', 'z', 'u', 'U', 'v', 'V'.", call. = FALSE)
  }
  if (ptype == "z") {
    ptype <- "h"
  }

  con_list <- .get_model_constituents(filename)
  infer_minor_constituents <- TRUE
  dots <- list(...)
  if (!is.null(dots$constituents) && is.null(constituents)) constituents <- dots$constituents
  if (!is.null(dots$coasts)) coasts <- dots$coasts
  if (!is.null(dots$InferMinor) && is.null(InferMinor)) InferMinor <- dots$InferMinor
  if (!is.null(constituents)) {
    con_list <- as.character(constituents)
    infer_minor_constituents <- FALSE
  }
  if (!is.null(InferMinor)) infer_minor_constituents <- isTRUE(InferMinor)

  time_num <- .matlab_datenum(time)
  input_grid_size <- if (is.null(dim(lat))) c(length(lat), 1) else dim(lat)
  input_time_size <- if (is.null(dim(time_num))) c(length(time_num), 1) else dim(time_num)
  lat_vec <- as.vector(lat)
  time_vec <- as.vector(time_num)

  map_solution <- length(lat_vec) > 1 && length(lat_vec) != length(time_vec)

  hc <- tmd_interp(filename, ptype, lat, lon, constituents = con_list, coasts = coasts, input_mode = "paired")
  hc_mat <- .reshape_hc_to_matrix(hc, length(con_list))
  astro <- tmd_astrol(time_vec)
  const_info <- tmd_constit(con_list)

  if (map_solution) {
    z <- matrix(NA_real_, nrow = length(lat_vec), ncol = length(time_vec))
    isf <- apply(hc_mat, 2, function(v) all(is.finite(v)))
    for (k in seq_along(time_vec)) {
      hhat <- tmd_harp(time_vec[k], hc_mat[, isf, drop = FALSE], con_list, astro$p[k], astro$N[k], const_info$ph, const_info$omega)
      d_minor <- 0
      if (infer_minor_constituents) {
        d_minor <- tmd_infer_minor(hc_mat[, isf, drop = FALSE], con_list, time_vec[k], astro$s[k], astro$h[k], astro$p[k], astro$N[k])
      }
      z[isf, k] <- d_minor + hhat
    }
  } else {
    hhat <- tmd_harp(time_vec, hc_mat, con_list, astro$p, astro$N, const_info$ph, const_info$omega)
    d_minor <- 0
    if (infer_minor_constituents) {
      d_minor <- tmd_infer_minor(hc_mat, con_list, time_vec, astro$s, astro$h, astro$p, astro$N)
    }
    z <- d_minor + hhat
  }

  if (identical(input_grid_size, c(1, 1))) {
    return(array(z, dim = input_time_size))
  }

  spatial_n <- prod(input_grid_size)
  trailing_n <- length(z) / spatial_n
  if (abs(trailing_n - round(trailing_n)) > .Machine$double.eps^0.5) {
    stop("Cannot reshape prediction output to the requested spatial/time dimensions.", call. = FALSE)
  }
  trailing_n <- as.integer(round(trailing_n))

  if (trailing_n <= 1) {
    return(array(z, dim = input_grid_size))
  }

  array(z, dim = c(input_grid_size[1], input_grid_size[2], trailing_n))
}
