tmd <- function(TMDFunctionName = "tmd") {
  docs <- c(
    tmd = "Available translated R functions: tmd_data, tmd_interp, tmd_predict, tmd_conlist, tmd_ellipse, tidal_range, tmd_astrol, tmd_constit, tmd_nodal, tmd_harp, tmd_infer_minor, tmd_ll2ps, tmd_ps2ll.",
    predict = "Use tmd_predict(filename, lat, lon, time, ptype = 'h', ...).",
    tmd_predict = "Use tmd_predict(filename, lat, lon, time, ptype = 'h', ...).",
    interp = "Use tmd_interp(filename, variable, lati, loni, ...).",
    tmd_interp = "Use tmd_interp(filename, variable, lati, loni, ...).",
    data = "Use tmd_data(filename, variable, ...).",
    tmd_data = "Use tmd_data(filename, variable, ...).",
    conlist = "Use tmd_conlist(filename).",
    tmd_conlist = "Use tmd_conlist(filename).",
    ellipse = "Use tmd_ellipse(filename, constituent, lati, loni).",
    tmd_ellipse = "Use tmd_ellipse(filename, constituent, lati, loni)."
  )
  key <- tolower(TMDFunctionName)
  msg <- docs[[key]]
  if (is.null(msg)) {
    warning(sprintf("No translated R help entry for '%s'.", TMDFunctionName), call. = FALSE)
    msg <- docs[["tmd"]]
  }
  message(msg)
  invisible(msg)
}

tmd_conlist <- function(filename) {
  .assert_tmd_file(filename)
  .get_model_constituents(filename)
}

tmd_data <- function(filename, variable, ..., constituents = NULL, bounds = NULL, geo = FALSE) {
  .assert_tmd_file(filename)
  dots <- list(...)
  if (!is.null(dots$constituents) && is.null(constituents)) {
    constituents <- dots$constituents
  }
  if (!is.null(dots$bounds) && is.null(bounds)) {
    bounds <- dots$bounds
  }
  if (!is.null(dots$geo)) {
    geo <- isTRUE(dots$geo)
  }

  con_list <- .get_model_constituents(filename)
  global_model <- .is_global_model(filename)

  if (!is.null(bounds)) {
    bounds <- as.matrix(bounds)
    if (ncol(bounds) != 2) {
      stop("Spatial bounds must be Mx2 in the form [xi yi] or [lon lat].", call. = FALSE)
    }
    if (global_model) {
      bounds[, 1] <- ifelse(bounds[, 1] < 0, bounds[, 1] + 360, bounds[, 1])
    }
  }

  if (is.null(constituents)) {
    constituents <- con_list
  } else {
    constituents <- as.character(constituents)
    if (!all(constituents %in% con_list)) {
      stop("Some requested constituents do not exist in the model file.", call. = FALSE)
    }
  }

  if (global_model) {
    x_or_lon <- as.numeric(ncdf4::ncvar_get(filename, "lon"))
    y_or_lat <- as.numeric(ncdf4::ncvar_get(filename, "lat"))
  } else {
    x_or_lon <- as.numeric(ncdf4::ncvar_get(filename, "x"))
    y_or_lat <- as.numeric(ncdf4::ncvar_get(filename, "y"))
  }

  dx <- abs(diff(x_or_lon[1:2]))
  dy <- abs(diff(y_or_lat[1:2]))
  if (!is.null(bounds)) {
    xl <- c(min(bounds[, 1]) - 2 * dx, max(bounds[, 1]) + 2 * dx)
    yl <- c(min(bounds[, 2]) - 2 * dy, max(bounds[, 2]) + 2 * dy)
    ri <- which(y_or_lat >= yl[1] & y_or_lat <= yl[2])
    ci <- which(x_or_lon >= xl[1] & x_or_lon <= xl[2])
    x_or_lon <- x_or_lon[ci]
    y_or_lat <- y_or_lat[ri]
  } else {
    ri <- seq_along(y_or_lat)
    ci <- seq_along(x_or_lon)
  }
  if (length(ci) == 0 || length(ri) == 0) {
    stop("No tide data available in the region of interest.", call. = FALSE)
  }

  con_index <- match(constituents, con_list)
  read_2d <- function(varname) {
    .permute_xy(.ncvar_get_subset(filename, varname, start = c(ci[1], ri[1]), count = c(length(ci), length(ri))))
  }
  read_3d <- function(varname, cind) {
    .permute_xy(.ncvar_get_subset(filename, varname, start = c(ci[1], ri[1], cind), count = c(length(ci), length(ri), length(cind))))
  }

  Z <- switch(
    variable,
    mask = as.logical(read_2d("mask")),
    flexure = read_2d("flexure") / 100,
    wct = read_2d("wct"),
    hRe = read_3d("hRe", con_index),
    hIm = read_3d("hIm", con_index),
    URe = read_3d("URe", con_index),
    UIm = read_3d("UIm", con_index),
    VRe = read_3d("VRe", con_index),
    VIm = read_3d("VIm", con_index),
    uRe = read_3d("URe", con_index),
    uIm = read_3d("UIm", con_index),
    vRe = read_3d("VRe", con_index),
    vIm = read_3d("VIm", con_index),
    h = complex(real = read_3d("hRe", con_index), imaginary = read_3d("hIm", con_index)),
    hAm = Mod(complex(real = read_3d("hRe", con_index), imaginary = read_3d("hIm", con_index))),
    hPh = Arg(complex(real = read_3d("hRe", con_index), imaginary = read_3d("hIm", con_index))),
    u = complex(real = read_3d("URe", con_index), imaginary = read_3d("UIm", con_index)),
    U = complex(real = read_3d("URe", con_index), imaginary = read_3d("UIm", con_index)),
    uAm = Mod(complex(real = read_3d("URe", con_index), imaginary = read_3d("UIm", con_index))),
    UAm = Mod(complex(real = read_3d("URe", con_index), imaginary = read_3d("UIm", con_index))),
    uPh = Arg(complex(real = read_3d("URe", con_index), imaginary = read_3d("UIm", con_index))),
    UPh = Arg(complex(real = read_3d("URe", con_index), imaginary = read_3d("UIm", con_index))),
    v = complex(real = read_3d("VRe", con_index), imaginary = read_3d("VIm", con_index)),
    V = complex(real = read_3d("VRe", con_index), imaginary = read_3d("VIm", con_index)),
    vAm = Mod(complex(real = read_3d("VRe", con_index), imaginary = read_3d("VIm", con_index))),
    VAm = Mod(complex(real = read_3d("VRe", con_index), imaginary = read_3d("VIm", con_index))),
    vPh = Arg(complex(real = read_3d("VRe", con_index), imaginary = read_3d("VIm", con_index))),
    VPh = Arg(complex(real = read_3d("VRe", con_index), imaginary = read_3d("VIm", con_index))),
    stop(sprintf("The requested variable %s is not available in %s.", variable, filename), call. = FALSE)
  )

  if (variable %in% c("uRe", "uIm", "u", "uAm", "vRe", "vIm", "v", "vAm")) {
    wct <- pmax(read_2d("wct"), 10)
    Z <- .array_divide_by_matrix(Z, wct)
  }

  if (geo && !global_model) {
    x_or_lon <- .permute_xy(.ncvar_get_subset(filename, "lon", start = c(ci[1], ri[1]), count = c(length(ci), length(ri))))
    y_or_lat <- .permute_xy(.ncvar_get_subset(filename, "lat", start = c(ci[1], ri[1]), count = c(length(ci), length(ri))))
  }

  list(Z = Z, x_or_lon = x_or_lon, y_or_lat = y_or_lat, conList = constituents)
}

tmd_ll2ps <- function(lat, lon, SLAT, SLON, HEMI) {
  CDR <- 57.29577951
  E2 <- 6.694379852e-3
  E <- sqrt(E2)
  pi_single <- 3.141592654
  RE <- 6378.1370

  if (abs(SLAT) == 90) {
    RHO <- 2 * RE / (((1 + E)^(1 + E) * (1 - E)^(1 - E))^(E / 2))
  } else {
    SL <- abs(SLAT) / CDR
    TC <- tan(pi_single / 4 - SL / 2) / (((1 - E * sin(SL)) / (1 + E * sin(SL)))^(E / 2))
    MC <- cos(SL) / sqrt(1 - E2 * (sin(SL)^2))
    RHO <- RE * MC / TC
  }
  lat_abs <- abs(lat) / CDR
  T <- tan(pi_single / 4 - lat_abs / 2) / (((1 - E * sin(lat_abs)) / (1 + E * sin(lat_abs)))^(E / 2))
  lon_r <- -(lon - SLON) / CDR
  x <- -RHO * T * sin(lon_r)
  y <- RHO * T * cos(lon_r)
  if (toupper(HEMI) == "N") {
    y <- -y
  }
  list(x = x, y = y)
}

tmd_ps2ll <- function(X, Y, SLAT, SLON, HEMI) {
  CDR <- 57.29577951
  E2 <- 6.694379852e-3
  E <- sqrt(E2)
  pi_single <- 3.141592654
  RE <- 6378.1370
  SGN <- if (toupper(HEMI) == "S") -1 else 1
  if (toupper(HEMI) == "N") {
    Y <- -Y
  }
  SLAT <- abs(SLAT)
  SL <- SLAT / CDR
  RHO <- sqrt(X^2 + Y^2)
  lat <- rep(NA_real_, length(RHO))
  lon <- rep(NA_real_, length(RHO))
  near_pole <- RHO < 0.1
  lat[near_pole] <- 90 * SGN
  lon[near_pole] <- 0
  use <- !near_pole
  if (any(use)) {
    CM <- cos(SL) / sqrt(1 - E2 * (sin(SL)^2))
    T <- tan((pi_single / 4) - (SL / 2)) / (((1 - E * sin(SL)) / (1 + E * sin(SL)))^(E / 2))
    if (abs(SLAT - 90) < 1e-5) {
      T <- RHO[use] * sqrt((1 + E)^(1 + E) * (1 - E)^(1 - E)) / (2 * RE)
    } else {
      T <- RHO[use] * T / (RE * CM)
    }
    a1 <- 5 * E2^2 / 24
    a2 <- E2^3 / 12
    a3 <- 7 * E2^2 / 48
    a4 <- 29 * E2^3 / 240
    a5 <- 7 * E2^3 / 120
    CHI <- (pi_single / 2) - 2 * atan(T)
    lat[use] <- CHI + ((E2 / 2) + a1 + a2) * sin(2 * CHI) + (a3 + a4) * sin(4 * CHI) + a5 * sin(6 * CHI)
    lat[use] <- SGN * lat[use] * CDR
    lon[use] <- -(atan2(-X[use], Y[use]) * CDR) + SLON
    lon[lon < -180] <- lon[lon < -180] + 360
    lon[lon > 180] <- lon[lon > 180] - 360
  }
  list(lat = lat, lon = lon)
}

tmd_interp <- function(filename, variable, lati, loni, ..., constituents = NULL, coasts = "nan") {
  .assert_tmd_file(filename)
  if (!is.character(variable) || length(variable) != 1) {
    stop("Input variable must be a string.", call. = FALSE)
  }
  if (!.islatlon(lati, loni)) {
    stop("Inputs lati and loni must be valid geographic coordinates.", call. = FALSE)
  }

  dots <- list(...)
  if (!is.null(dots$constituents) && is.null(constituents)) {
    constituents <- dots$constituents
  }
  if (!is.null(dots$coasts)) {
    coasts <- dots$coasts
  }
  if (is.numeric(coasts) && length(coasts) == 1 && is.na(coasts)) {
    coasts <- "nan"
  }
  coasts <- tolower(as.character(coasts))

  proj4 <- .proj4_string(filename)
  if (identical(proj4, "+proj=stere +lat_0=90 +lat_ts=70 +lon_0=-45 +x_0=0 +y_0=0 +datum=WGS84 +units=km +no_defs +type=crs")) {
    xy <- tmd_ll2ps(lati, loni, 70, -45, "N")
    xi <- xy$x
    yi <- xy$y
    xi[lati < 0] <- NA
    yi[lati < 0] <- NA
  } else if (identical(proj4, "+proj=stere +lat_0=-90 +lat_ts=-71 +lon_0=-70 +x_0=0 +y_0=0 +datum=WGS84 +units=km +no_defs +type=crs")) {
    xy <- tmd_ll2ps(lati, loni, -71, -70, "S")
    xi <- xy$x
    yi <- xy$y
    xi[lati > 0] <- NA
    yi[lati > 0] <- NA
  } else if (identical(proj4, "+proj=stere +lat_0=90 +lat_ts=70 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=km +no_defs +type=crs")) {
    xy <- tmd_ll2ps(lati, loni, 70, 0, "N")
    xi <- xy$x
    yi <- xy$y
  } else {
    xi <- loni
    yi <- lati
    xi[xi < 0] <- xi[xi < 0] + 360
  }

  variable_lower <- tolower(variable)
  subset_bounds <- cbind(as.vector(xi), as.vector(yi))
  data_obj <- switch(
    variable_lower,
    hph = tmd_data(filename, "h", bounds = subset_bounds, constituents = constituents),
    uph = tmd_data(filename, "U", bounds = subset_bounds, constituents = constituents),
    vph = tmd_data(filename, "V", bounds = subset_bounds, constituents = constituents),
    wct = tmd_data(filename, variable, bounds = subset_bounds),
    flexure = tmd_data(filename, variable, bounds = subset_bounds),
    mask = tmd_data(filename, variable, bounds = subset_bounds),
    tmd_data(filename, variable, bounds = subset_bounds, constituents = constituents)
  )
  z <- data_obj$Z
  x <- data_obj$x_or_lon
  y <- data_obj$y_or_lat

  if (variable == "mask") {
    zi <- .interp_nearest(x, y, z, xi, yi, extrap = 2)
  } else {
    zi <- .interpolate_slices(x, y, z, xi, yi, method = "linear", extrap = NA_real_)
  }
  if (variable_lower %in% c("ham", "uam", "vam")) {
    zi <- Mod(zi)
  } else if (variable_lower %in% c("hph", "uph", "vph")) {
    zi <- Arg(zi)
  }

  if (!variable %in% c("mask", "flexure", "wct")) {
    if (coasts == "nan") {
      mask <- tmd_data(filename, "mask", bounds = subset_bounds)$Z
      maski <- .interp_nearest(x, y, as.numeric(mask != 1), xi, yi, extrap = 2) == 0
      for (k in seq_len(dim(zi)[3])) {
        tmp <- zi[, , k]
        tmp[!maski] <- NA
        zi[, , k] <- tmp
      }
    } else if (coasts %in% c("flex", "flexure")) {
      flex <- tmd_data(filename, "flexure", bounds = subset_bounds)$Z
      flexi <- .interp_linear(x, y, flex, xi, yi, extrap = NA_real_)
      for (k in seq_len(dim(zi)[3])) {
        zi[, , k] <- zi[, , k] * flexi
      }
    }
  }

  .drop_last_dim_if_scalar(zi)
}
