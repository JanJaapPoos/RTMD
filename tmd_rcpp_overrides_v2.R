# Avoid coercion warnings by using the Rcpp backend only for real-valued interpolation.

.interp_linear_multi <- function(x, y, z, xi, yi, extrap = NA_real_) {
  if (is.complex(z) || !.tmd_init_rcpp()) {
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
