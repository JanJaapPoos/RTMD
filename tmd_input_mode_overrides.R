# Explicit paired/grid input handling for public entry points.

.gridify_lat_lon <- function(lat, lon, input_mode = c("paired", "grid")) {
  input_mode <- match.arg(input_mode)
  if (input_mode == "grid" && is.null(dim(lat)) && is.null(dim(lon)) && length(lat) > 1 && length(lon) > 1) {
    lat_vals <- as.vector(lat)
    lon_vals <- as.vector(lon)
    lat <- matrix(rep(lat_vals, times = length(lon_vals)), nrow = length(lat_vals), ncol = length(lon_vals))
    lon <- matrix(rep(lon_vals, each = length(lat_vals)), nrow = length(lat_vals), ncol = length(lon_vals))
  }
  list(lat = lat, lon = lon)
}

tmd_interp <- function(filename, variable, lati, loni, ..., constituents = NULL, coasts = "nan", input_mode = c("paired", "grid")) {
  grids <- .gridify_lat_lon(lati, loni, input_mode = input_mode)
  lati <- grids$lat
  loni <- grids$lon

  .assert_tmd_file(filename)
  if (!is.character(variable) || length(variable) != 1) stop("Input variable must be a string.", call. = FALSE)
  if (!.islatlon(lati, loni)) stop("Inputs lati and loni must be valid geographic coordinates.", call. = FALSE)
  dots <- list(...)
  if (!is.null(dots$constituents) && is.null(constituents)) constituents <- dots$constituents
  if (!is.null(dots$coasts)) coasts <- dots$coasts
  if (is.numeric(coasts) && length(coasts) == 1 && is.na(coasts)) coasts <- "nan"
  coasts <- tolower(as.character(coasts))
  proj4 <- .proj4_string(filename)
  if (identical(proj4, "+proj=stere +lat_0=90 +lat_ts=70 +lon_0=-45 +x_0=0 +y_0=0 +datum=WGS84 +units=km +no_defs +type=crs")) {
    xy <- tmd_ll2ps(lati, loni, 70, -45, "N"); xi <- xy$x; yi <- xy$y; xi[lati < 0] <- NA; yi[lati < 0] <- NA
  } else if (identical(proj4, "+proj=stere +lat_0=-90 +lat_ts=-71 +lon_0=-70 +x_0=0 +y_0=0 +datum=WGS84 +units=km +no_defs +type=crs")) {
    xy <- tmd_ll2ps(lati, loni, -71, -70, "S"); xi <- xy$x; yi <- xy$y; xi[lati > 0] <- NA; yi[lati > 0] <- NA
  } else if (identical(proj4, "+proj=stere +lat_0=90 +lat_ts=70 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=km +no_defs +type=crs")) {
    xy <- tmd_ll2ps(lati, loni, 70, 0, "N"); xi <- xy$x; yi <- xy$y
  } else {
    xi <- loni; yi <- lati; xi[xi < 0] <- xi[xi < 0] + 360
  }
  subset_bounds <- cbind(as.vector(xi), as.vector(yi))
  vl <- tolower(variable)
  data_obj <- switch(vl,
    hph = tmd_data(filename, "h", bounds = subset_bounds, constituents = constituents),
    uph = tmd_data(filename, "U", bounds = subset_bounds, constituents = constituents),
    vph = tmd_data(filename, "V", bounds = subset_bounds, constituents = constituents),
    wct = tmd_data(filename, variable, bounds = subset_bounds),
    flexure = tmd_data(filename, variable, bounds = subset_bounds),
    mask = tmd_data(filename, variable, bounds = subset_bounds),
    tmd_data(filename, variable, bounds = subset_bounds, constituents = constituents)
  )
  z <- data_obj$Z; x <- data_obj$x_or_lon; y <- data_obj$y_or_lat
  zi <- if (variable == "mask") .interp_nearest(x, y, z, xi, yi, extrap = 2) else .interp_linear_multi(x, y, z, xi, yi, extrap = NA_real_)
  if (vl %in% c("ham", "uam", "vam")) zi <- Mod(zi) else if (vl %in% c("hph", "uph", "vph")) zi <- Arg(zi)
  if (!variable %in% c("mask", "flexure", "wct")) {
    if (coasts == "nan") {
      mask <- tmd_data(filename, "mask", bounds = subset_bounds)$Z
      maski <- .interp_nearest(x, y, as.numeric(mask != 1), xi, yi, extrap = 2) == 0
      if (length(dim(zi)) == 2) zi[!as.vector(maski), ] <- NA else zi[!maski] <- NA
    } else if (coasts %in% c("flex", "flexure")) {
      flex <- tmd_data(filename, "flexure", bounds = subset_bounds)$Z
      flexi <- .interp_linear(x, y, flex, xi, yi, extrap = NA_real_)
      if (length(dim(zi)) == 2) zi <- zi * as.vector(flexi) else zi <- zi * flexi
    }
  }
  .drop_last_dim_if_scalar(zi)
}

tmd_predict <- function(filename, lat, lon, time, ptype = "h", ..., constituents = NULL, coasts = "nan", InferMinor = NULL, input_mode = c("paired", "grid")) {
  grids <- .gridify_lat_lon(lat, lon, input_mode = input_mode)
  lat <- grids$lat
  lon <- grids$lon

  .assert_tmd_file(filename)
  if (!identical(dim(lat), dim(lon)) && length(lat) != length(lon)) stop("Dimensions of lat and lon must agree.", call. = FALSE)
  if (!ptype %in% c("h", "z", "u", "U", "v", "V")) stop("ptype must be one of 'h', 'z', 'u', 'U', 'v', 'V'.", call. = FALSE)
  if (ptype == "z") ptype <- "h"
  con_list <- .get_model_constituents(filename)
  infer_minor_constituents <- TRUE
  dots <- list(...)
  if (!is.null(dots$constituents) && is.null(constituents)) constituents <- dots$constituents
  if (!is.null(dots$coasts)) coasts <- dots$coasts
  if (!is.null(dots$InferMinor) && is.null(InferMinor)) InferMinor <- dots$InferMinor
  if (!is.null(constituents)) { con_list <- as.character(constituents); infer_minor_constituents <- FALSE }
  if (!is.null(InferMinor)) infer_minor_constituents <- isTRUE(InferMinor)
  time_num <- .matlab_datenum(time)
  input_grid_size <- if (is.null(dim(lat))) c(length(lat), 1) else dim(lat)
  input_time_size <- if (is.null(dim(time_num))) c(length(time_num), 1) else dim(time_num)
  lat_vec <- as.vector(lat); time_vec <- as.vector(time_num)
  map_solution <- length(lat_vec) > 1 && length(lat_vec) != length(time_vec)
  hc <- tmd_interp(filename, ptype, lat, lon, constituents = con_list, coasts = coasts, input_mode = "paired")
  hc_mat <- .reshape_hc_to_matrix(hc, length(con_list))
  astro <- tmd_astrol(time_vec); const_info <- tmd_constit(con_list)
  if (map_solution) {
    z <- matrix(NA_real_, nrow = length(lat_vec), ncol = length(time_vec))
    isf <- apply(hc_mat, 2, function(v) all(is.finite(v)))
    for (k in seq_along(time_vec)) {
      hhat <- tmd_harp(time_vec[k], hc_mat[, isf, drop = FALSE], con_list, astro$p[k], astro$N[k], const_info$ph, const_info$omega)
      d_minor <- 0
      if (infer_minor_constituents) d_minor <- tmd_infer_minor(hc_mat[, isf, drop = FALSE], con_list, time_vec[k], astro$s[k], astro$h[k], astro$p[k], astro$N[k])
      z[isf, k] <- d_minor + hhat
    }
  } else {
    hhat <- tmd_harp(time_vec, hc_mat, con_list, astro$p, astro$N, const_info$ph, const_info$omega)
    d_minor <- 0
    if (infer_minor_constituents) d_minor <- tmd_infer_minor(hc_mat, con_list, time_vec, astro$s, astro$h, astro$p, astro$N)
    z <- d_minor + hhat
  }
  if (identical(input_grid_size, c(1, 1))) return(array(z, dim = input_time_size))
  array(z, dim = c(input_grid_size[1], input_grid_size[2], length(time_vec)))
}

tmd_ellipse <- function(filename, constituent, lati, loni, input_mode = c("paired", "grid")) {
  grids <- .gridify_lat_lon(lati, loni, input_mode = input_mode)
  lati <- grids$lat
  loni <- grids$lon

  if (!is.character(constituent) || length(constituent) != 1) stop("Constituent must be a single string.", call. = FALSE)
  u <- .drop_last_dim_if_scalar(tmd_interp(filename, "u", lati, loni, constituents = constituent, input_mode = "paired"))
  v <- .drop_last_dim_if_scalar(tmd_interp(filename, "v", lati, loni, constituents = constituent, input_mode = "paired"))
  t1p <- Re(u) + Im(v); t2p <- Re(v) - Im(u); t1m <- Re(u) - Im(v); t2m <- Re(v) + Im(u)
  ap <- sqrt(t1p^2 + t2p^2) / 2; am <- sqrt(t1m^2 + t2m^2) / 2
  ep <- atan2(t2p, t1p); ep <- 180 * (ep + 2 * pi * (ep < 0)) / pi
  em <- atan2(t2m, t1m); em <- 180 * (em + 2 * pi * (em < 0)) / pi
  umajor <- ap + am; uminor <- ap - am; uincl <- 0.5 * (em + ep); uincl <- uincl - 180 * (uincl > 180)
  uphase <- 0.5 * (em - ep); uphase <- uphase + 360 * (uphase < 0); uphase <- uphase - 360 * (uphase >= 360)
  list(umajor = umajor, uminor = uminor, uphase = uphase, uincl = uincl)
}
