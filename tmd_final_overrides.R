# Final consolidated overrides for the translated TMD R implementation.

.assert_tmd_file <- function(filename) {
  .require_ncdf4()
  if (!grepl("\\.nc$", filename, ignore.case = TRUE)) {
    stop("Input filename must end in .nc.", call. = FALSE)
  }
  if (!file.exists(filename)) {
    stop(sprintf("Cannot find %s. Check the path and try again.", filename), call. = FALSE)
  }

  version_ok <- FALSE
  structure_ok <- FALSE
  try({
    nc <- ncdf4::nc_open(filename)
    on.exit(ncdf4::nc_close(nc), add = TRUE)
    version_val <- tryCatch(ncdf4::ncatt_get(nc, 0, "tmd_version")$value, error = function(e) NULL)
    version_num <- suppressWarnings(as.numeric(version_val))
    version_ok <- !is.null(version_num) && !is.na(version_num) && version_num >= 3.0
    proj4_val <- tryCatch(ncdf4::ncatt_get(nc, "mapping", "spatial_proj4")$value, error = function(e) NULL)
    cons_val <- tryCatch(ncdf4::ncatt_get(nc, "constituents", "constituent_order")$value, error = function(e) NULL)
    structure_ok <- !is.null(proj4_val) && !is.null(cons_val) && nzchar(trimws(cons_val))
  }, silent = TRUE)

  if (!(version_ok || structure_ok)) {
    stop(sprintf("%s is not compatible with TMD3.0+.", filename), call. = FALSE)
  }
}

.get_att <- function(filename, varid, attname) {
  nc <- ncdf4::nc_open(filename)
  on.exit(ncdf4::nc_close(nc), add = TRUE)
  val <- ncdf4::ncatt_get(nc, varid, attname)$value
  if (is.null(val)) {
    stop(sprintf("Missing NetCDF attribute '%s' for '%s'.", attname, varid), call. = FALSE)
  }
  val
}

.ncvar_get_subset <- function(filename, varname, start = NULL, count = NULL) {
  nc <- ncdf4::nc_open(filename)
  on.exit(ncdf4::nc_close(nc), add = TRUE)
  if (is.null(start) || is.null(count)) {
    return(ncdf4::ncvar_get(nc, varname, collapse_degen = FALSE))
  }
  ncdf4::ncvar_get(nc, varname, start = start, count = count, collapse_degen = FALSE)
}

.nc_get_att_open <- function(nc, varid, attname) {
  val <- ncdf4::ncatt_get(nc, varid, attname)$value
  if (is.null(val)) {
    stop(sprintf("Missing NetCDF attribute '%s' for '%s'.", attname, varid), call. = FALSE)
  }
  val
}

.read_runs <- function(idx) {
  idx <- sort(as.integer(idx))
  if (length(idx) == 0) {
    return(list())
  }
  split(idx, cumsum(c(TRUE, diff(idx) != 1L)))
}

.complex_with_dim <- function(re, im) {
  re3 <- .as_3d_array(re)
  im3 <- .as_3d_array(im)
  array(complex(real = as.vector(re3), imaginary = as.vector(im3)), dim = dim(re3))
}

.gridify_lat_lon <- function(lat, lon) {
  if (is.null(dim(lat)) && is.null(dim(lon)) && length(lat) > 1 && length(lon) > 1) {
    lat_vals <- as.vector(lat)
    lon_vals <- as.vector(lon)
    lat <- matrix(rep(lat_vals, times = length(lon_vals)), nrow = length(lat_vals), ncol = length(lon_vals))
    lon <- matrix(rep(lon_vals, each = length(lat_vals)), nrow = length(lat_vals), ncol = length(lon_vals))
  }
  list(lat = lat, lon = lon)
}

.query_shape <- function(x) {
  d <- dim(x)
  if (is.null(d)) {
    list(dim = c(length(x)), vector = as.vector(x), is_scalar = length(x) == 1)
  } else {
    list(dim = d, vector = as.vector(x), is_scalar = FALSE)
  }
}

.restore_query_shape <- function(values, shape) {
  if (length(shape$dim) == 1) {
    if (shape$is_scalar) {
      return(values[1])
    }
    return(array(values, dim = shape$dim))
  }
  array(values, dim = shape$dim)
}

.normalize_grid_data <- function(z, x, y) {
  d <- dim(z)
  nx <- length(x)
  ny <- length(y)
  if (is.null(d)) {
    if (length(z) == nx * ny) {
      return(array(z, dim = c(ny, nx)))
    }
    stop("Grid data dimensions do not match x/y coordinates.", call. = FALSE)
  }
  if (length(d) == 2) {
    if (all(d == c(ny, nx))) return(z)
    if (all(d == c(nx, ny))) return(t(z))
    if (prod(d) == nx * ny) return(array(as.vector(z), dim = c(ny, nx)))
    stop("Grid data dimensions do not match x/y coordinates.", call. = FALSE)
  }
  if (length(d) >= 3) {
    pairs <- utils::combn(seq_along(d), 2)
    for (j in seq_len(ncol(pairs))) {
      a <- pairs[1, j]
      b <- pairs[2, j]
      if (all(d[c(a, b)] == c(ny, nx))) return(aperm(z, c(a, b, setdiff(seq_along(d), c(a, b)))))
      if (all(d[c(a, b)] == c(nx, ny))) return(aperm(z, c(b, a, setdiff(seq_along(d), c(a, b)))))
    }
    n_other <- prod(d) / (nx * ny)
    if (is.finite(n_other) && abs(n_other - round(n_other)) < .Machine$double.eps^0.5) {
      return(array(as.vector(z), dim = c(ny, nx, as.integer(round(n_other)))))
    }
    stop("Grid data dimensions do not match x/y coordinates.", call. = FALSE)
  }
  stop("Unsupported grid dimensionality.", call. = FALSE)
}

.prepare_interp_grid <- function(x, y, z) {
  if (is.matrix(x) || is.matrix(y)) {
    stop("Interpolation requires vector x/y grids.", call. = FALSE)
  }
  x <- as.numeric(x)
  y <- as.numeric(y)
  if (anyNA(x) || anyNA(y)) {
    stop("Interpolation grids cannot contain NA values.", call. = FALSE)
  }
  z <- .normalize_grid_data(z, x, y)
  if (length(x) > 1 && is.unsorted(x, strictly = FALSE)) {
    ox <- order(x)
    x <- x[ox]
    if (length(dim(z)) == 2) z <- z[, ox, drop = FALSE] else z <- z[, ox, , drop = FALSE]
  }
  if (length(y) > 1 && is.unsorted(y, strictly = FALSE)) {
    oy <- order(y)
    y <- y[oy]
    if (length(dim(z)) == 2) z <- z[oy, , drop = FALSE] else z <- z[oy, , , drop = FALSE]
  }
  list(x = x, y = y, z = z)
}

.nearest_indices <- function(grid, pts) {
  idx <- findInterval(pts, grid)
  idx[idx < 1] <- 1L
  idx[idx >= length(grid)] <- length(grid) - 1L
  left <- grid[idx]
  right <- grid[idx + 1L]
  idx[abs(pts - right) < abs(pts - left)] <- idx[abs(pts - right) < abs(pts - left)] + 1L
  idx[is.na(pts)] <- NA_integer_
  idx
}

.interp_nearest <- function(x, y, z, xi, yi, extrap = NA_real_) {
  g <- .prepare_interp_grid(x, y, z)
  shape <- .query_shape(xi)
  xi_v <- shape$vector
  yi_v <- as.vector(yi)
  out <- rep(extrap, length(xi_v))
  inside <- !is.na(xi_v) & !is.na(yi_v) & xi_v >= min(g$x) & xi_v <= max(g$x) & yi_v >= min(g$y) & yi_v <= max(g$y)
  if (any(inside)) {
    ix <- .nearest_indices(g$x, xi_v[inside])
    iy <- .nearest_indices(g$y, yi_v[inside])
    out[inside] <- g$z[cbind(iy, ix)]
  }
  .restore_query_shape(out, shape)
}

.interp_linear <- function(x, y, z, xi, yi, extrap = NA_real_) {
  g <- .prepare_interp_grid(x, y, z)
  shape <- .query_shape(xi)
  xi_v <- shape$vector
  yi_v <- as.vector(yi)
  out <- rep(extrap, length(xi_v))
  inside <- !is.na(xi_v) & !is.na(yi_v) & xi_v >= min(g$x) & xi_v <= max(g$x) & yi_v >= min(g$y) & yi_v <= max(g$y)
  if (!any(inside)) return(.restore_query_shape(out, shape))
  ix <- findInterval(xi_v[inside], g$x)
  iy <- findInterval(yi_v[inside], g$y)
  ix[ix < 1] <- 1L; iy[iy < 1] <- 1L
  ix[ix >= length(g$x)] <- length(g$x) - 1L
  iy[iy >= length(g$y)] <- length(g$y) - 1L
  x1 <- g$x[ix]; x2 <- g$x[ix + 1L]; y1 <- g$y[iy]; y2 <- g$y[iy + 1L]
  tx <- ifelse(x2 == x1, 0, (xi_v[inside] - x1) / (x2 - x1))
  ty <- ifelse(y2 == y1, 0, (yi_v[inside] - y1) / (y2 - y1))
  z11 <- g$z[cbind(iy, ix)]
  z21 <- g$z[cbind(iy, ix + 1L)]
  z12 <- g$z[cbind(iy + 1L, ix)]
  z22 <- g$z[cbind(iy + 1L, ix + 1L)]
  out[inside] <- (1 - tx) * (1 - ty) * z11 + tx * (1 - ty) * z21 + (1 - tx) * ty * z12 + tx * ty * z22
  .restore_query_shape(out, shape)
}

.interp_linear_multi <- function(x, y, z, xi, yi, extrap = NA_real_) {
  g <- .prepare_interp_grid(x, y, z)
  z3 <- .as_3d_array(g$z)
  shape <- .query_shape(xi)
  xi_v <- shape$vector
  yi_v <- as.vector(yi)
  ns <- dim(z3)[3]
  out <- matrix(extrap, nrow = length(xi_v), ncol = ns)
  inside <- !is.na(xi_v) & !is.na(yi_v) & xi_v >= min(g$x) & xi_v <= max(g$x) & yi_v >= min(g$y) & yi_v <= max(g$y)
  if (any(inside)) {
    ix <- findInterval(xi_v[inside], g$x)
    iy <- findInterval(yi_v[inside], g$y)
    ix[ix < 1] <- 1L; iy[iy < 1] <- 1L
    ix[ix >= length(g$x)] <- length(g$x) - 1L
    iy[iy >= length(g$y)] <- length(g$y) - 1L
    x1 <- g$x[ix]; x2 <- g$x[ix + 1L]; y1 <- g$y[iy]; y2 <- g$y[iy + 1L]
    tx <- ifelse(x2 == x1, 0, (xi_v[inside] - x1) / (x2 - x1))
    ty <- ifelse(y2 == y1, 0, (yi_v[inside] - y1) / (y2 - y1))
    w11 <- (1 - tx) * (1 - ty); w21 <- tx * (1 - ty); w12 <- (1 - tx) * ty; w22 <- tx * ty
    nrowz <- dim(z3)[1]
    idx11 <- iy + (ix - 1L) * nrowz
    idx21 <- iy + ix * nrowz
    idx12 <- (iy + 1L) + (ix - 1L) * nrowz
    idx22 <- (iy + 1L) + ix * nrowz
    zmat <- matrix(z3, nrow = nrowz * dim(z3)[2], ncol = ns)
    out[inside, ] <- w11 * zmat[idx11, , drop = FALSE] + w21 * zmat[idx21, , drop = FALSE] +
      w12 * zmat[idx12, , drop = FALSE] + w22 * zmat[idx22, , drop = FALSE]
  }
  array(out, dim = c(shape$dim, ns))
}

tmd_conlist <- function(filename) {
  .assert_tmd_file(filename)
  .get_model_constituents(filename)
}

tmd_data <- function(filename, variable, ..., constituents = NULL, bounds = NULL, geo = FALSE) {
  .assert_tmd_file(filename)
  dots <- list(...)
  if (!is.null(dots$constituents) && is.null(constituents)) constituents <- dots$constituents
  if (!is.null(dots$bounds) && is.null(bounds)) bounds <- dots$bounds
  if (!is.null(dots$geo)) geo <- isTRUE(dots$geo)
  nc <- ncdf4::nc_open(filename)
  on.exit(ncdf4::nc_close(nc), add = TRUE)
  con_list <- .split_constituents(.nc_get_att_open(nc, "constituents", "constituent_order"))
  global_model <- identical(.nc_get_att_open(nc, "mapping", "spatial_proj4"), "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
  if (!is.null(bounds)) {
    bounds <- as.matrix(bounds)
    if (ncol(bounds) != 2) stop("Spatial bounds must be Mx2 in the form [xi yi] or [lon lat].", call. = FALSE)
    if (global_model) bounds[, 1] <- ifelse(bounds[, 1] < 0, bounds[, 1] + 360, bounds[, 1])
  }
  if (is.null(constituents)) {
    constituents <- con_list
  } else {
    constituents <- as.character(constituents)
    if (!all(constituents %in% con_list)) stop("Some requested constituents do not exist in the model file.", call. = FALSE)
  }
  if (global_model) {
    x_or_lon <- as.numeric(ncdf4::ncvar_get(nc, "lon", collapse_degen = FALSE))
    y_or_lat <- as.numeric(ncdf4::ncvar_get(nc, "lat", collapse_degen = FALSE))
  } else {
    x_or_lon <- as.numeric(ncdf4::ncvar_get(nc, "x", collapse_degen = FALSE))
    y_or_lat <- as.numeric(ncdf4::ncvar_get(nc, "y", collapse_degen = FALSE))
  }
  dx <- abs(diff(x_or_lon[1:2])); dy <- abs(diff(y_or_lat[1:2]))
  if (!is.null(bounds)) {
    xl <- c(min(bounds[, 1]) - 2 * dx, max(bounds[, 1]) + 2 * dx)
    yl <- c(min(bounds[, 2]) - 2 * dy, max(bounds[, 2]) + 2 * dy)
    ri <- which(y_or_lat >= yl[1] & y_or_lat <= yl[2]); ci <- which(x_or_lon >= xl[1] & x_or_lon <= xl[2])
    x_or_lon <- x_or_lon[ci]; y_or_lat <- y_or_lat[ri]
  } else {
    ri <- seq_along(y_or_lat); ci <- seq_along(x_or_lon)
  }
  if (length(ci) == 0 || length(ri) == 0) stop("No tide data available in the region of interest.", call. = FALSE)
  con_index <- match(constituents, con_list)
  read_2d <- function(varname) .permute_xy(ncdf4::ncvar_get(nc, varname, start = c(ci[1], ri[1]), count = c(length(ci), length(ri)), collapse_degen = FALSE))
  read_3d <- function(varname, cind) {
    runs <- .read_runs(cind)
    out <- array(NA_real_, dim = c(length(ri), length(ci), length(cind)))
    pos <- 1L
    for (run in runs) {
      raw <- ncdf4::ncvar_get(nc, varname, start = c(ci[1], ri[1], run[1]), count = c(length(ci), length(ri), length(run)), collapse_degen = FALSE)
      block <- .as_3d_array(.permute_xy(raw))
      out[, , pos:(pos + length(run) - 1L)] <- block
      pos <- pos + length(run)
    }
    out
  }
  h_parts <- NULL; u_parts <- NULL; v_parts <- NULL
  get_h <- function() { if (is.null(h_parts)) h_parts <<- list(re = read_3d("hRe", con_index), im = read_3d("hIm", con_index)); h_parts }
  get_u <- function() { if (is.null(u_parts)) u_parts <<- list(re = read_3d("URe", con_index), im = read_3d("UIm", con_index)); u_parts }
  get_v <- function() { if (is.null(v_parts)) v_parts <<- list(re = read_3d("VRe", con_index), im = read_3d("VIm", con_index)); v_parts }
  Z <- switch(variable,
    mask = as.logical(read_2d("mask")),
    flexure = read_2d("flexure") / 100,
    wct = read_2d("wct"),
    hRe = get_h()$re, hIm = get_h()$im,
    URe = get_u()$re, UIm = get_u()$im,
    VRe = get_v()$re, VIm = get_v()$im,
    uRe = get_u()$re, uIm = get_u()$im,
    vRe = get_v()$re, vIm = get_v()$im,
    h = .complex_with_dim(get_h()$re, get_h()$im),
    hAm = Mod(.complex_with_dim(get_h()$re, get_h()$im)),
    hPh = Arg(.complex_with_dim(get_h()$re, get_h()$im)),
    u = .complex_with_dim(get_u()$re, get_u()$im),
    U = .complex_with_dim(get_u()$re, get_u()$im),
    uAm = Mod(.complex_with_dim(get_u()$re, get_u()$im)),
    UAm = Mod(.complex_with_dim(get_u()$re, get_u()$im)),
    uPh = Arg(.complex_with_dim(get_u()$re, get_u()$im)),
    UPh = Arg(.complex_with_dim(get_u()$re, get_u()$im)),
    v = .complex_with_dim(get_v()$re, get_v()$im),
    V = .complex_with_dim(get_v()$re, get_v()$im),
    vAm = Mod(.complex_with_dim(get_v()$re, get_v()$im)),
    VAm = Mod(.complex_with_dim(get_v()$re, get_v()$im)),
    vPh = Arg(.complex_with_dim(get_v()$re, get_v()$im)),
    VPh = Arg(.complex_with_dim(get_v()$re, get_v()$im)),
    stop(sprintf("The requested variable %s is not available in %s.", variable, filename), call. = FALSE)
  )
  if (variable %in% c("uRe", "uIm", "u", "uAm", "vRe", "vIm", "v", "vAm")) Z <- .array_divide_by_matrix(Z, pmax(read_2d("wct"), 10))
  if (geo && !global_model) {
    x_or_lon <- .permute_xy(ncdf4::ncvar_get(nc, "lon", start = c(ci[1], ri[1]), count = c(length(ci), length(ri)), collapse_degen = FALSE))
    y_or_lat <- .permute_xy(ncdf4::ncvar_get(nc, "lat", start = c(ci[1], ri[1]), count = c(length(ci), length(ri)), collapse_degen = FALSE))
  }
  list(Z = Z, x_or_lon = x_or_lon, y_or_lat = y_or_lat, conList = constituents)
}

tmd_interp <- function(filename, variable, lati, loni, ..., constituents = NULL, coasts = "nan") {
  grids <- .gridify_lat_lon(lati, loni); lati <- grids$lat; loni <- grids$lon
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

tmd_predict <- function(filename, lat, lon, time, ptype = "h", ..., constituents = NULL, coasts = "nan", InferMinor = NULL) {
  grids <- .gridify_lat_lon(lat, lon); lat <- grids$lat; lon <- grids$lon
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
  hc <- tmd_interp(filename, ptype, lat, lon, constituents = con_list, coasts = coasts)
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

tmd_ellipse <- function(filename, constituent, lati, loni) {
  grids <- .gridify_lat_lon(lati, loni); lati <- grids$lat; loni <- grids$lon
  if (!is.character(constituent) || length(constituent) != 1) stop("Constituent must be a single string.", call. = FALSE)
  u <- .drop_last_dim_if_scalar(tmd_interp(filename, "u", lati, loni, constituents = constituent))
  v <- .drop_last_dim_if_scalar(tmd_interp(filename, "v", lati, loni, constituents = constituent))
  t1p <- Re(u) + Im(v); t2p <- Re(v) - Im(u); t1m <- Re(u) - Im(v); t2m <- Re(v) + Im(u)
  ap <- sqrt(t1p^2 + t2p^2) / 2; am <- sqrt(t1m^2 + t2m^2) / 2
  ep <- atan2(t2p, t1p); ep <- 180 * (ep + 2 * pi * (ep < 0)) / pi
  em <- atan2(t2m, t1m); em <- 180 * (em + 2 * pi * (em < 0)) / pi
  umajor <- ap + am; uminor <- ap - am; uincl <- 0.5 * (em + ep); uincl <- uincl - 180 * (uincl > 180)
  uphase <- 0.5 * (em - ep); uphase <- uphase + 360 * (uphase < 0); uphase <- uphase - 360 * (uphase >= 360)
  list(umajor = umajor, uminor = uminor, uphase = uphase, uincl = uincl)
}
