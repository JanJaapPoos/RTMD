tmd_predict <- function(filename, lat, lon, time, ptype = "h", ..., constituents = NULL, coasts = "nan", InferMinor = NULL) {
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
  if (!is.null(dots$constituents) && is.null(constituents)) {
    constituents <- dots$constituents
  }
  if (!is.null(dots$coasts)) {
    coasts <- dots$coasts
  }
  if (!is.null(dots$InferMinor) && is.null(InferMinor)) {
    InferMinor <- dots$InferMinor
  }
  if (!is.null(constituents)) {
    con_list <- as.character(constituents)
    infer_minor_constituents <- FALSE
  }
  if (!is.null(InferMinor)) {
    infer_minor_constituents <- isTRUE(InferMinor)
  }

  time_num <- .matlab_datenum(time)
  input_grid_size <- if (is.null(dim(lat))) c(length(lat), 1) else dim(lat)
  input_time_size <- if (is.null(dim(time_num))) c(length(time_num), 1) else dim(time_num)
  lat_vec <- as.vector(lat)
  time_vec <- as.vector(time_num)
  map_solution <- length(lat_vec) > 1 && length(lat_vec) != length(time_vec)

  hc <- tmd_interp(filename, ptype, lat, lon, constituents = con_list, coasts = coasts)
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
  array(z, dim = c(input_grid_size[1], input_grid_size[2], length(time_vec)))
}

tmd_ellipse <- function(filename, constituent, lati, loni) {
  if (!is.character(constituent) || length(constituent) != 1) {
    stop("Constituent must be a single string.", call. = FALSE)
  }
  u <- tmd_interp(filename, "u", lati, loni, constituents = constituent)
  v <- tmd_interp(filename, "v", lati, loni, constituents = constituent)
  u <- .drop_last_dim_if_scalar(u)
  v <- .drop_last_dim_if_scalar(v)

  t1p <- Re(u) + Im(v)
  t2p <- Re(v) - Im(u)
  t1m <- Re(u) - Im(v)
  t2m <- Re(v) + Im(u)
  ap <- sqrt(t1p^2 + t2p^2) / 2
  am <- sqrt(t1m^2 + t2m^2) / 2

  ep <- atan2(t2p, t1p)
  ep <- ep + 2 * pi * (ep < 0)
  ep <- 180 * ep / pi
  em <- atan2(t2m, t1m)
  em <- em + 2 * pi * (em < 0)
  em <- 180 * em / pi

  umajor <- ap + am
  uminor <- ap - am
  uincl <- 0.5 * (em + ep)
  uincl <- uincl - 180 * (uincl > 180)
  uphase <- 0.5 * (em - ep)
  uphase <- uphase + 360 * (uphase < 0)
  uphase <- uphase - 360 * (uphase >= 360)
  list(umajor = umajor, uminor = uminor, uphase = uphase, uincl = uincl)
}

cube2rect <- function(A3, mask = NULL) {
  d <- dim(A3)
  if (length(d) != 3) {
    stop("A3 must be a 3D array.", call. = FALSE)
  }
  if (!is.null(mask)) {
    if (!all(dim(mask) == d[1:2])) {
      stop("Dimensions of mask must match the first two dimensions of A3.", call. = FALSE)
    }
    if (!is.logical(mask)) {
      stop("mask must be logical.", call. = FALSE)
    }
  }
  A2 <- aperm(A3, c(3, 1, 2))
  dim(A2) <- c(d[3], d[1] * d[2])
  if (!is.null(mask)) {
    A2 <- A2[, mask, drop = FALSE]
  }
  A2
}

rect2cube <- function(A2, gridsize_or_mask) {
  if (is.logical(gridsize_or_mask) && length(dim(gridsize_or_mask)) == 2) {
    mask <- gridsize_or_mask
    gridsize <- dim(mask)
    full <- if (is.logical(A2)) matrix(FALSE, nrow = nrow(A2), ncol = length(mask)) else matrix(NA, nrow = nrow(A2), ncol = length(mask))
    full[, mask] <- A2
    return(aperm(array(full, dim = c(nrow(A2), gridsize[1], gridsize[2])), c(2, 3, 1)))
  }
  gridsize <- as.integer(gridsize_or_mask)
  if (length(gridsize) == 2) {
    if (gridsize[1] * gridsize[2] != ncol(A2)) {
      stop("gridsize must match the number of columns in A2.", call. = FALSE)
    }
    return(aperm(array(A2, dim = c(nrow(A2), gridsize[1], gridsize[2])), c(2, 3, 1)))
  }
  if (length(gridsize) == 3) {
    return(aperm(array(A2, dim = c(gridsize[3], gridsize[1], gridsize[2])), c(2, 3, 1)))
  }
  stop("gridsize_or_mask must be either a 2D logical mask or a 2/3 element size vector.", call. = FALSE)
}

tidal_range <- function(filename_or_hc, conList = NULL, mask = NULL) {
  t <- seq(730486, 730852, by = 1 / 48)
  if (is.character(filename_or_hc) && length(filename_or_hc) == 1) {
    .assert_tmd_file(filename_or_hc)
    conList <- .get_model_constituents(filename_or_hc)
    mask <- tmd_data(filename_or_hc, "mask")$Z
    hc <- tmd_data(filename_or_hc, "h")$Z
  } else {
    hc <- filename_or_hc
    if (is.null(conList) || is.null(mask)) {
      stop("If the first input is numeric, conList and mask must also be supplied.", call. = FALSE)
    }
  }

  hc_rect <- cube2rect(hc, mask)
  astro <- tmd_astrol(t)
  const_info <- tmd_constit(conList)
  z_min <- rep(0, ncol(hc_rect))
  z_max <- rep(0, ncol(hc_rect))

  for (k in seq_along(t)) {
    hhat <- tmd_harp(t[k], t(hc_rect), conList, astro$p[k], astro$N[k], const_info$ph, const_info$omega)
    d_minor <- tmd_infer_minor(t(hc_rect), conList, t[k], astro$s[k], astro$h[k], astro$p[k], astro$N[k])
    z_now <- d_minor + hhat
    z_max <- pmax(z_max, z_now)
    z_min <- pmin(z_min, z_now)
  }

  R <- rect2cube(matrix(z_max - z_min, nrow = 1), mask)
  R[is.na(R)] <- 0
  .drop_last_dim_if_scalar(R)
}
