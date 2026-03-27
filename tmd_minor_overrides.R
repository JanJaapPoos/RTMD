# Fix dimension handling in minor-constituent inference for vector times.

tmd_infer_minor <- function(hc, constituents, t, s, h, p, N) {
  constituents <- tolower(as.character(constituents))
  t <- .matlab_datenum(t)
  hc <- as.matrix(hc)
  M <- ncol(hc)
  Tn <- length(t)
  if (Tn > 1 && !(M == 1 || M == Tn)) {
    stop("If neither times nor location are scalar, they must be vectors of the same length.", call. = FALSE)
  }

  t_1992 <- t - 727564
  t_hr <- (t_1992 - floor(t_1992)) * 24
  t1 <- 15 * t_hr
  t2 <- 30 * t_hr
  cid8 <- c("q1","o1","p1","k1","n2","m2","s2","k2")
  Lib <- match(cid8, constituents)
  Lia <- !is.na(Lib)
  z8 <- hc[Lib, , drop = FALSE]
  if (sum(Lia) <= 5) {
    stop("Not enough constituents for inference.", call. = FALSE)
  }

  minor_cons <- c("2q1","sigma1","rho1","m1","m1","chi1","pi1","phi1","theta1","j1","oo1","2n2","mu2","nu2","lambda2","l2","l2","t2")
  ind <- !(minor_cons %in% constituents)
  rad <- pi / 180
  PP <- 282.8

  zmin <- matrix(0 + 0i, nrow = 18, ncol = M)
  zmin[1,] <- 0.263*z8[1,] - 0.0252*z8[2,]
  zmin[2,] <- 0.297*z8[1,] - 0.0264*z8[2,]
  zmin[3,] <- 0.164*z8[1,] + 0.0048*z8[2,]
  zmin[4,] <- 0.0140*z8[2,] + 0.0101*z8[4,]
  zmin[5,] <- 0.0389*z8[2,] + 0.0282*z8[4,]
  zmin[6,] <- 0.0064*z8[2,] + 0.0060*z8[4,]
  zmin[7,] <- 0.0030*z8[2,] + 0.0171*z8[4,]
  zmin[8,] <- -0.0015*z8[2,] + 0.0152*z8[4,]
  zmin[9,] <- -0.0065*z8[2,] + 0.0155*z8[4,]
  zmin[10,] <- -0.0389*z8[2,] + 0.0836*z8[4,]
  zmin[11,] <- -0.0431*z8[2,] + 0.0613*z8[4,]
  zmin[12,] <- 0.264*z8[5,] - 0.0253*z8[6,]
  zmin[13,] <- 0.298*z8[5,] - 0.0264*z8[6,]
  zmin[14,] <- 0.165*z8[5,] + 0.00487*z8[6,]
  zmin[15,] <- 0.0040*z8[6,] + 0.0074*z8[7,]
  zmin[16,] <- 0.0131*z8[6,] + 0.0326*z8[7,]
  zmin[17,] <- 0.0033*z8[6,] + 0.0082*z8[7,]
  zmin[18,] <- 0.0585*z8[7,]

  arg <- matrix(0, nrow = 18, ncol = Tn)
  arg[1,] <- t1 - 4*s + h + 2*p - 90
  arg[2,] <- t1 - 4*s + 3*h - 90
  arg[3,] <- t1 - 3*s + 3*h - p - 90
  arg[4,] <- t1 - s + h - p + 90
  arg[5,] <- t1 - s + h + p + 90
  arg[6,] <- t1 - s + 3*h - p + 90
  arg[7,] <- t1 - 2*h + PP - 90
  arg[8,] <- t1 + 3*h + 90
  arg[9,] <- t1 + s - h + p + 90
  arg[10,] <- t1 + s + h - p + 90
  arg[11,] <- t1 + 2*s + h + 90
  arg[12,] <- t2 - 4*s + 2*h + 2*p
  arg[13,] <- t2 - 4*s + 4*h
  arg[14,] <- t2 - 3*s + 4*h - p
  arg[15,] <- t2 - s + p + 180
  arg[16,] <- t2 - s + 2*h - p + 180
  arg[17,] <- t2 - s + 2*h + p
  arg[18,] <- t2 - h + PP

  sinn <- sin(N * rad)
  cosn <- cos(N * rad)
  sin2n <- sin(2 * N * rad)
  cos2n <- cos(2 * N * rad)

  f <- matrix(1, nrow = 18, ncol = Tn)
  f[1,] <- sqrt((1 + 0.189*cosn - 0.0058*cos2n)^2 + (0.189*sinn - 0.0058*sin2n)^2)
  f[2,] <- f[1,]; f[3,] <- f[1,]
  f[4,] <- sqrt((1 + 0.185*cosn)^2 + (0.185*sinn)^2)
  f[5,] <- sqrt((1 + 0.201*cosn)^2 + (0.201*sinn)^2)
  f[6,] <- sqrt((1 + 0.221*cosn)^2 + (0.221*sinn)^2)
  f[10,] <- sqrt((1 + 0.198*cosn)^2 + (0.198*sinn)^2)
  f[11,] <- sqrt((1 + 0.640*cosn + 0.134*cos2n)^2 + (0.640*sinn + 0.134*sin2n)^2)
  f[12,] <- sqrt((1 - 0.0373*cosn)^2 + (0.0373*sinn)^2)
  f[13,] <- f[12,]; f[14,] <- f[12,]; f[16,] <- f[12,]
  f[17,] <- sqrt((1 + 0.441*cosn)^2 + (0.441*sinn)^2)

  u <- matrix(0, nrow = 18, ncol = Tn)
  u[1,] <- atan2(0.189*sinn - 0.0058*sin2n, 1 + 0.189*cosn - 0.0058*sin2n) / rad
  u[2,] <- u[1,]; u[3,] <- u[1,]
  u[4,] <- atan2(0.185*sinn, 1 + 0.185*cosn) / rad
  u[5,] <- atan2(-0.201*sinn, 1 + 0.201*cosn) / rad
  u[6,] <- atan2(-0.221*sinn, 1 + 0.221*cosn) / rad
  u[10,] <- atan2(-0.198*sinn, 1 + 0.198*cosn) / rad
  u[11,] <- atan2(-0.640*sinn - 0.134*sin2n, 1 + 0.640*cosn + 0.134*cos2n) / rad
  u[12,] <- atan2(-0.0373*sinn, 1 - 0.0373*cosn) / rad
  u[13,] <- u[12,]; u[14,] <- u[12,]; u[16,] <- u[12,]
  u[17,] <- atan2(-0.441*sinn, 1 + 0.441*cosn) / rad

  if (Tn == 1) {
    tmp <- (arg[ind, 1] + u[ind, 1]) * rad
    return(as.vector(colSums(
      Re(zmin[ind, , drop = FALSE]) * matrix(f[ind, 1] * cos(tmp), nrow = sum(ind), ncol = M) +
      Im(zmin[ind, , drop = FALSE]) * matrix(f[ind, 1] * sin(tmp), nrow = sum(ind), ncol = M)
    )))
  }

  if (M == 1) {
    tmp <- (arg[ind, , drop = FALSE] + u[ind, , drop = FALSE]) * rad
    re_mat <- matrix(Re(zmin[ind, 1]), nrow = sum(ind), ncol = Tn)
    im_mat <- matrix(Im(zmin[ind, 1]), nrow = sum(ind), ncol = Tn)
    return(colSums(
      re_mat * (f[ind, , drop = FALSE] * cos(tmp)) +
      im_mat * (f[ind, , drop = FALSE] * sin(tmp))
    ))
  }

  tmp <- (arg[ind, , drop = FALSE] + u[ind, , drop = FALSE]) * rad
  rowSums(
    t(Re(zmin[ind, , drop = FALSE])) * t(f[ind, , drop = FALSE] * cos(tmp)) +
    t(Im(zmin[ind, , drop = FALSE])) * t(f[ind, , drop = FALSE] * sin(tmp))
  )
}
