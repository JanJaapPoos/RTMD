# Shared helpers for the TMD MATLAB-to-R translation.

.require_ncdf4 <- function() {
  if (!requireNamespace("ncdf4", quietly = TRUE)) {
    stop("Package 'ncdf4' is required. Install it with install.packages('ncdf4').", call. = FALSE)
  }
}

.matlab_datenum <- function(x) {
  if (inherits(x, "POSIXt")) {
    return(as.numeric(as.POSIXct(x, tz = "UTC")) / 86400 + 719529)
  }
  if (inherits(x, "Date")) {
    return(as.numeric(x) + 719529)
  }
  if (is.numeric(x)) {
    return(x)
  }
  stop("Time input must be numeric, Date, or POSIXt.", call. = FALSE)
}

.split_constituents <- function(x) {
  if (length(x) == 0 || is.na(x) || !nzchar(x)) {
    return(character())
  }
  strsplit(trimws(x), "\\s+")[[1]]
}

.assert_tmd_file <- function(filename) {
  .require_ncdf4()
  if (!grepl("\\.nc$", filename, ignore.case = TRUE)) {
    stop("Input filename must end in .nc.", call. = FALSE)
  }
  if (!file.exists(filename)) {
    stop(sprintf("Cannot find %s. Check the path and try again.", filename), call. = FALSE)
  }
  tmd_version <- tryCatch(
    ncdf4::ncatt_get(filename, 0, "tmd_version")$value,
    error = function(e) NULL
  )
  if (is.null(tmd_version) || is.na(as.numeric(tmd_version)) || as.numeric(tmd_version) < 3.0) {
    stop(sprintf("%s is not compatible with TMD3.0+.", filename), call. = FALSE)
  }
}

.get_att <- function(filename, varid, attname) {
  val <- ncdf4::ncatt_get(filename, varid, attname)$value
  if (is.null(val)) {
    stop(sprintf("Missing NetCDF attribute '%s' for '%s'.", attname, varid), call. = FALSE)
  }
  val
}

.get_model_constituents <- function(filename) {
  .split_constituents(.get_att(filename, "constituents", "constituent_order"))
}

.proj4_string <- function(filename) {
  .get_att(filename, "mapping", "spatial_proj4")
}

.is_global_model <- function(filename) {
  identical(.proj4_string(filename), "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
}

.ncvar_get_subset <- function(filename, varname, start = NULL, count = NULL) {
  if (is.null(start) || is.null(count)) {
    return(ncdf4::ncvar_get(filename, varname))
  }
  ncdf4::ncvar_get(filename, varname, start = start, count = count)
}

.permute_xy <- function(x) {
  d <- dim(x)
  if (is.null(d)) {
    return(x)
  }
  if (length(d) == 2) {
    return(t(x))
  }
  if (length(d) == 3) {
    return(aperm(x, c(2, 1, 3)))
  }
  stop("Unsupported NetCDF variable dimensionality.", call. = FALSE)
}

.as_3d_array <- function(x) {
  d <- dim(x)
  if (is.null(d)) {
    return(array(x, dim = c(length(x), 1, 1)))
  }
  if (length(d) == 2) {
    return(array(x, dim = c(d[1], d[2], 1)))
  }
  if (length(d) == 3) {
    return(x)
  }
  stop("Unsupported array dimensionality.", call. = FALSE)
}

.array_divide_by_matrix <- function(arr, mat) {
  arr3 <- .as_3d_array(arr)
  out <- arr3
  for (k in seq_len(dim(arr3)[3])) {
    out[, , k] <- arr3[, , k] / mat
  }
  out
}

.drop_last_dim_if_scalar <- function(x) {
  d <- dim(x)
  if (!is.null(d) && length(d) >= 1 && d[length(d)] == 1) {
    return(array(x, dim = d[-length(d)]))
  }
  x
}

.reshape_hc_to_matrix <- function(hc, n_constituents) {
  d <- dim(hc)
  if (is.null(d) || length(d) == 1) {
    return(matrix(hc, nrow = n_constituents, ncol = 1))
  }
  perm <- c(length(d), seq_len(length(d) - 1))
  tmp <- aperm(hc, perm)
  matrix(tmp, nrow = n_constituents, ncol = length(tmp) / n_constituents)
}

.islatlon <- function(lat, lon) {
  is.numeric(lat) &&
    is.numeric(lon) &&
    !any(abs(lat) > 90, na.rm = TRUE) &&
    !any(lon > 360, na.rm = TRUE) &&
    !any(lon < -180, na.rm = TRUE)
}

.prepare_interp_grid <- function(x, y, z) {
  if (is.matrix(x) || is.matrix(y)) {
    stop("Interpolation requires vector x/y grids.", call. = FALSE)
  }
  list(x = as.numeric(x), y = as.numeric(y), z = z)
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
  out <- array(extrap, dim = dim(xi))
  inside <- !is.na(xi) & !is.na(yi) &
    xi >= min(g$x) & xi <= max(g$x) &
    yi >= min(g$y) & yi <= max(g$y)
  if (!any(inside)) {
    return(out)
  }
  ix <- .nearest_indices(g$x, xi[inside])
  iy <- .nearest_indices(g$y, yi[inside])
  out[inside] <- g$z[cbind(iy, ix)]
  out
}

.interp_linear <- function(x, y, z, xi, yi, extrap = NA_real_) {
  g <- .prepare_interp_grid(x, y, z)
  out <- array(extrap, dim = dim(xi))
  inside <- !is.na(xi) & !is.na(yi) &
    xi >= min(g$x) & xi <= max(g$x) &
    yi >= min(g$y) & yi <= max(g$y)
  if (!any(inside)) {
    return(out)
  }

  xi_v <- xi[inside]
  yi_v <- yi[inside]
  ix <- findInterval(xi_v, g$x)
  iy <- findInterval(yi_v, g$y)
  ix[ix < 1] <- 1L
  iy[iy < 1] <- 1L
  ix[ix >= length(g$x)] <- length(g$x) - 1L
  iy[iy >= length(g$y)] <- length(g$y) - 1L

  x1 <- g$x[ix]
  x2 <- g$x[ix + 1L]
  y1 <- g$y[iy]
  y2 <- g$y[iy + 1L]

  tx <- ifelse(x2 == x1, 0, (xi_v - x1) / (x2 - x1))
  ty <- ifelse(y2 == y1, 0, (yi_v - y1) / (y2 - y1))

  z11 <- g$z[cbind(iy, ix)]
  z21 <- g$z[cbind(iy, ix + 1L)]
  z12 <- g$z[cbind(iy + 1L, ix)]
  z22 <- g$z[cbind(iy + 1L, ix + 1L)]

  out[inside] <- (1 - tx) * (1 - ty) * z11 +
    tx * (1 - ty) * z21 +
    (1 - tx) * ty * z12 +
    tx * ty * z22
  out
}

.interpolate_slices <- function(x, y, z, xi, yi, method = c("linear", "nearest"), extrap = NA_real_) {
  method <- match.arg(method)
  z3 <- .as_3d_array(z)
  out <- array(NA, dim = c(dim(xi), dim(z3)[3]))
  for (k in seq_len(dim(z3)[3])) {
    if (method == "nearest") {
      out[, , k] <- .interp_nearest(x, y, z3[, , k], xi, yi, extrap = extrap)
    } else {
      out[, , k] <- .interp_linear(x, y, z3[, , k], xi, yi, extrap = extrap)
    }
  }
  out
}
