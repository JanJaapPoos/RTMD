tmd_astrol <- function(t) {
  t <- .matlab_datenum(t)
  t_mjd <- t - 727564 + 48622
  T <- t_mjd - 51544.4993
  list(
    s = (218.3164 + 13.17639648 * T) %% 360,
    h = (280.4661 + 0.98564736 * T) %% 360,
    p = (83.3535 + 0.11140353 * T) %% 360,
    N = (125.0445 - 0.05295377 * T) %% 360
  )
}

tmd_constit <- function(constituents) {
  c_data <- c("m2","s2","k1","o1","n2","p1","k2","q1","2n2","mu2","nu2","l2","t2","j1","m1","oo1","rho1","mf","mm","ssa","m4","ms4","mn4","m6","m8","mk3","s6","2sm2","2mk3","s1","2q1","m3","sa")
  ispec_data <- c(2,2,1,1,2,1,2,1,2,2,2,2,2,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,3,0)
  amp_data <- c(0.2441,0.112743,0.141565,0.100661,0.046397,0.046848,0.030684,0.019273,0.006141,0.007408,0.008811,0.006931,0.006608,0.007915,0.007915,0.004338,0.003661,0.042041,0.022191,0.019567,0,0,0,0,0,0,0,0,0,7.6464e-04,0.002565,0.003192,0.0031)
  phase_data <- c(1.731557546,0,0.173003674,1.558553872,6.050721243,6.110181633,3.487600001,5.877717569,4.086699633,3.463115091,5.427136701,0.553986502,0.052841931,2.137025284,2.436575100,1.929046130,5.254133027,1.756042456,1.964021610,3.487600001,3.463115091,1.731557546,1.499093481,5.194672637,6.926230184,1.904561220,0,4.551627762,3.809122439,0,3.913707,5.738991,6.2303435066085)
  omega_data <- c(1.405189e-04,1.454441e-04,7.292117e-05,6.759774e-05,1.378797e-04,7.252295e-05,1.458423e-04,6.495854e-05,1.352405e-04,1.355937e-04,1.382329e-04,1.431581e-04,1.452450e-04,7.556036e-05,7.028195e-05,7.824458e-05,6.531174e-05,0.053234e-04,0.026392e-04,0.003982e-04,2.810377e-04,2.859630e-04,2.783984e-04,4.215566e-04,5.620755e-04,2.134402e-04,4.363323e-04,1.503693e-04,2.081166e-04,7.2722e-05,0.6231934e-04,2.107783523e-04,1.990968403205456e-07)
  alpha_data <- c(0.693,0.693,0.736,0.695,0.693,0.706,0.693,0.695,0.693,0.693,0.693,0.693,0.693,0.695,0.695,0.695,0.695,0.693,0.693,0.693,0.693,0.693,0.693,0.693,0.693,0.693,0.693,0.693,0.693,0.693,0.693,0.802,0.693)
  kk <- match(tolower(as.character(constituents)), c_data)
  if (all(is.na(kk))) {
    return(list(ispec = -1, amp = 0, ph = 0, omega = 0, alpha = 0))
  }
  list(ispec = ispec_data[kk], amp = amp_data[kk], ph = phase_data[kk], omega = omega_data[kk], alpha = alpha_data[kk])
}

tmd_nodal <- function(constituents, p, N) {
  cid0 <- c("sa","ssa","mm","msf","mf","mt","alpha1","2q1","sigma1","q1","rho1","o1","tau1","m1","chi1","pi1","p1","s1","k1","psi1","phi1","theta1","j1","oo1","2n2","mu2","n2","nu2","m2a","m2","m2b","lambda2","l2","t2","s2","r2","k2","eta2","mns2","2sm2","m3","mk3","s3","mn4","m4","ms4","mk4","s4","s5","m6","s6","s7","s8")
  rad <- pi / 180
  nT <- length(p)
  sinn <- sin(N * rad)
  cosn <- cos(N * rad)
  sin2n <- sin(2 * N * rad)
  cos2n <- cos(2 * N * rad)
  sin3n <- sin(3 * N * rad)

  f <- matrix(1, nrow = nT, ncol = 53)
  f[,3] <- 1 - 0.130 * cosn
  f[,5] <- 1.043 + 0.414 * cosn
  f[,6] <- sqrt((1 + 0.203 * cosn + 0.040 * cos2n)^2 + (0.203 * sinn + 0.040 * sin2n)^2)
  f[,8] <- sqrt((1 + 0.188 * cosn)^2 + (0.188 * sinn)^2)
  f[,9] <- f[,8]; f[,10] <- f[,8]; f[,11] <- f[,8]
  f[,12] <- sqrt((1 + 0.189 * cosn - 0.0058 * cos2n)^2 + (0.189 * sinn - 0.0058 * sin2n)^2)
  tmp1 <- 1.36 * cos(p * rad) + 0.267 * cos((p - N) * rad)
  tmp2 <- 0.64 * sin(p * rad) + 0.135 * sin((p - N) * rad)
  f[,14] <- sqrt(tmp1^2 + tmp2^2)
  f[,15] <- sqrt((1 + 0.221 * cosn)^2 + (0.221 * sinn)^2)
  f[,19] <- sqrt((1 + 0.1158 * cosn - 0.0029 * cos2n)^2 + (0.1554 * sinn - 0.0029 * sin2n)^2)
  f[,23] <- sqrt((1 + 0.169 * cosn)^2 + (0.227 * sinn)^2)
  f[,24] <- sqrt((1 + 0.640 * cosn + 0.134 * cos2n)^2 + (0.640 * sinn + 0.134 * sin2n)^2)
  f[,25] <- sqrt((1 - 0.03731 * cosn + 0.00052 * cos2n)^2 + (0.03731 * sinn - 0.00052 * sin2n)^2)
  f[,26] <- f[,25]; f[,27] <- f[,25]; f[,28] <- f[,25]; f[,30] <- f[,25]
  temp1 <- 1 - 0.25 * cos(2 * p * rad) - 0.11 * cos((2 * p - N) * rad) - 0.04 * cosn
  temp2 <- 0.25 * sin(2 * p * rad) + 0.11 * sin((2 * p - N) * rad) + 0.04 * sinn
  f[,33] <- sqrt(temp1^2 + temp2^2)
  f[,37] <- sqrt((1 + 0.2852 * cosn + 0.0324 * cos2n)^2 + (0.3108 * sinn + 0.0324 * sin2n)^2)
  f[,38] <- sqrt((1 + 0.436 * cosn)^2 + (0.436 * sinn)^2)
  f[,39] <- f[,30]^2; f[,40] <- f[,30]; f[,42] <- f[,19] * f[,30]
  f[,44] <- f[,30]^2; f[,45] <- f[,44]; f[,46] <- f[,30]; f[,47] <- f[,30] * f[,37]; f[,50] <- f[,30]^3

  u <- matrix(0, nrow = nT, ncol = 53)
  u[,5] <- (-23.7 * sinn + 2.7 * sin2n - 0.4 * sin3n) * rad
  u[,6] <- atan(-(0.203 * sinn + 0.040 * sin2n) / (1 + 0.203 * cosn + 0.040 * cos2n))
  u[,8] <- atan(0.189 * sinn / (1 + 0.189 * cosn))
  u[,9] <- u[,8]; u[,10] <- u[,8]; u[,11] <- u[,8]
  u[,12] <- (10.8 * sinn - 1.3 * sin2n + 0.2 * sin3n) * rad
  u[,14] <- atan2(tmp2, tmp1)
  u[,15] <- atan(-0.221 * sinn / (1 + 0.221 * cosn))
  u[,19] <- atan((-0.1554 * sinn + 0.0029 * sin2n) / (1 + 0.1158 * cosn - 0.0029 * cos2n))
  u[,23] <- atan(-0.227 * sinn / (1 + 0.169 * cosn))
  u[,24] <- atan(-(0.640 * sinn + 0.134 * sin2n) / (1 + 0.640 * cosn + 0.134 * cos2n))
  u[,25] <- atan((-0.03731 * sinn + 0.00052 * sin2n) / (1 - 0.03731 * cosn + 0.00052 * cos2n))
  u[,26] <- u[,25]; u[,27] <- u[,25]; u[,28] <- u[,25]; u[,30] <- u[,25]
  u[,33] <- atan(-temp2 / temp1)
  u[,37] <- atan(-(0.3108 * sinn + 0.0324 * sin2n) / (1 + 0.2852 * cosn + 0.0324 * cos2n))
  u[,38] <- atan(-0.436 * sinn / (1 + 0.436 * cosn))
  u[,39] <- u[,30] * 2; u[,40] <- u[,30]; u[,41] <- 1.5 * u[,30]; u[,42] <- u[,30] + u[,19]
  u[,44] <- u[,30] * 2; u[,45] <- u[,44]; u[,46] <- u[,30]; u[,47] <- u[,30] + u[,37]; u[,50] <- u[,30] * 3

  Lib <- match(tolower(as.character(constituents)), cid0)
  list(pu = u[, Lib, drop = FALSE], pf = f[, Lib, drop = FALSE])
}

tmd_harp <- function(t, hc, constituents, p, N, ph, omega) {
  t <- .matlab_datenum(t)
  hc <- as.matrix(hc)
  time_s <- (t - 727564) * 86400
  nodal <- tmd_nodal(constituents, p, N)
  tmp <- outer(time_s, omega) + matrix(ph, nrow = length(t), ncol = length(constituents), byrow = TRUE) + nodal$pu
  w_re <- nodal$pf * cos(tmp)
  w_im <- nodal$pf * sin(tmp)

  if (length(t) == 1) {
    return(as.vector(colSums(Re(hc) * matrix(w_re[1, ], nrow = nrow(hc), ncol = ncol(hc)) +
                                Im(hc) * matrix(w_im[1, ], nrow = nrow(hc), ncol = ncol(hc)))))
  }
  if (ncol(hc) == 1) {
    return(drop(w_re %*% Re(hc[, 1, drop = FALSE]) + w_im %*% Im(hc[, 1, drop = FALSE])))
  }
  if (ncol(hc) == length(t)) {
    return(rowSums(t(Re(hc)) * w_re + t(Im(hc)) * w_im))
  }
  stop("hc must have one column, or one column per input time.", call. = FALSE)
}

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
  zmin[1,] <- 0.263*z8[1,] - 0.0252*z8[2,]; zmin[2,] <- 0.297*z8[1,] - 0.0264*z8[2,]; zmin[3,] <- 0.164*z8[1,] + 0.0048*z8[2,]
  zmin[4,] <- 0.0140*z8[2,] + 0.0101*z8[4,]; zmin[5,] <- 0.0389*z8[2,] + 0.0282*z8[4,]; zmin[6,] <- 0.0064*z8[2,] + 0.0060*z8[4,]
  zmin[7,] <- 0.0030*z8[2,] + 0.0171*z8[4,]; zmin[8,] <- -0.0015*z8[2,] + 0.0152*z8[4,]; zmin[9,] <- -0.0065*z8[2,] + 0.0155*z8[4,]
  zmin[10,] <- -0.0389*z8[2,] + 0.0836*z8[4,]; zmin[11,] <- -0.0431*z8[2,] + 0.0613*z8[4,]; zmin[12,] <- 0.264*z8[5,] - 0.0253*z8[6,]
  zmin[13,] <- 0.298*z8[5,] - 0.0264*z8[6,]; zmin[14,] <- 0.165*z8[5,] + 0.00487*z8[6,]; zmin[15,] <- 0.0040*z8[6,] + 0.0074*z8[7,]
  zmin[16,] <- 0.0131*z8[6,] + 0.0326*z8[7,]; zmin[17,] <- 0.0033*z8[6,] + 0.0082*z8[7,]; zmin[18,] <- 0.0585*z8[7,]

  arg <- matrix(0, nrow = 18, ncol = Tn)
  arg[1,] <- t1 - 4*s + h + 2*p - 90; arg[2,] <- t1 - 4*s + 3*h - 90; arg[3,] <- t1 - 3*s + 3*h - p - 90
  arg[4,] <- t1 - s + h - p + 90; arg[5,] <- t1 - s + h + p + 90; arg[6,] <- t1 - s + 3*h - p + 90
  arg[7,] <- t1 - 2*h + PP - 90; arg[8,] <- t1 + 3*h + 90; arg[9,] <- t1 + s - h + p + 90
  arg[10,] <- t1 + s + h - p + 90; arg[11,] <- t1 + 2*s + h + 90; arg[12,] <- t2 - 4*s + 2*h + 2*p
  arg[13,] <- t2 - 4*s + 4*h; arg[14,] <- t2 - 3*s + 4*h - p; arg[15,] <- t2 - s + p + 180
  arg[16,] <- t2 - s + 2*h - p + 180; arg[17,] <- t2 - s + 2*h + p; arg[18,] <- t2 - h + PP

  sinn <- sin(N * rad); cosn <- cos(N * rad); sin2n <- sin(2 * N * rad); cos2n <- cos(2 * N * rad)
  f <- matrix(1, nrow = 18, ncol = Tn)
  f[1,] <- sqrt((1 + 0.189*cosn - 0.0058*cos2n)^2 + (0.189*sinn - 0.0058*sin2n)^2)
  f[2,] <- f[1,]; f[3,] <- f[1,]; f[4,] <- sqrt((1 + 0.185*cosn)^2 + (0.185*sinn)^2); f[5,] <- sqrt((1 + 0.201*cosn)^2 + (0.201*sinn)^2)
  f[6,] <- sqrt((1 + 0.221*cosn)^2 + (0.221*sinn)^2); f[10,] <- sqrt((1 + 0.198*cosn)^2 + (0.198*sinn)^2)
  f[11,] <- sqrt((1 + 0.640*cosn + 0.134*cos2n)^2 + (0.640*sinn + 0.134*sin2n)^2)
  f[12,] <- sqrt((1 - 0.0373*cosn)^2 + (0.0373*sinn)^2); f[13,] <- f[12,]; f[14,] <- f[12,]; f[16,] <- f[12,]
  f[17,] <- sqrt((1 + 0.441*cosn)^2 + (0.441*sinn)^2)

  u <- matrix(0, nrow = 18, ncol = Tn)
  u[1,] <- atan2(0.189*sinn - 0.0058*sin2n, 1 + 0.189*cosn - 0.0058*sin2n) / rad
  u[2,] <- u[1,]; u[3,] <- u[1,]; u[4,] <- atan2(0.185*sinn, 1 + 0.185*cosn) / rad; u[5,] <- atan2(-0.201*sinn, 1 + 0.201*cosn) / rad
  u[6,] <- atan2(-0.221*sinn, 1 + 0.221*cosn) / rad; u[10,] <- atan2(-0.198*sinn, 1 + 0.198*cosn) / rad
  u[11,] <- atan2(-0.640*sinn - 0.134*sin2n, 1 + 0.640*cosn + 0.134*cos2n) / rad
  u[12,] <- atan2(-0.0373*sinn, 1 - 0.0373*cosn) / rad; u[13,] <- u[12,]; u[14,] <- u[12,]; u[16,] <- u[12,]; u[17,] <- atan2(-0.441*sinn, 1 + 0.441*cosn) / rad

  if (Tn == 1) {
    tmp <- (arg[ind, 1] + u[ind, 1]) * rad
    return(as.vector(colSums(Re(zmin[ind, , drop = FALSE]) * matrix(f[ind, 1] * cos(tmp), nrow = sum(ind), ncol = M) +
                                Im(zmin[ind, , drop = FALSE]) * matrix(f[ind, 1] * sin(tmp), nrow = sum(ind), ncol = M))))
  }
  if (M == 1) {
    tmp <- (arg[ind, , drop = FALSE] + u[ind, , drop = FALSE]) * rad
    return(colSums(Re(zmin[ind, 1, drop = FALSE]) * (f[ind, , drop = FALSE] * cos(tmp)) +
                     Im(zmin[ind, 1, drop = FALSE]) * (f[ind, , drop = FALSE] * sin(tmp))))
  }
  tmp <- (arg[ind, , drop = FALSE] + u[ind, , drop = FALSE]) * rad
  rowSums(t(Re(zmin[ind, , drop = FALSE])) * t(f[ind, , drop = FALSE] * cos(tmp)) +
            t(Im(zmin[ind, , drop = FALSE])) * t(f[ind, , drop = FALSE] * sin(tmp)))
}
