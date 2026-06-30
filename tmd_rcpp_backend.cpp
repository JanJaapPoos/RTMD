#include <Rcpp.h>
#include <algorithm>
#include <cmath>

using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix interp_linear_multi_cpp(
    NumericVector x,
    NumericVector y,
    NumericVector z,
    IntegerVector zdim,
    NumericVector xi,
    NumericVector yi,
    double extrap = NA_REAL) {

  const int ny = zdim[0];
  const int nx = zdim[1];
  const int ns = zdim[2];
  const int nq = xi.size();

  NumericMatrix out(nq, ns);
  std::fill(out.begin(), out.end(), extrap);

  const double xmin = x[0];
  const double xmax = x[x.size() - 1];
  const double ymin = y[0];
  const double ymax = y[y.size() - 1];
  const int plane = ny * nx;

  for (int q = 0; q < nq; ++q) {
    const double xq = xi[q];
    const double yq = yi[q];

    if (NumericVector::is_na(xq) || NumericVector::is_na(yq) ||
        xq < xmin || xq > xmax || yq < ymin || yq > ymax) {
      continue;
    }

    int ix = std::upper_bound(x.begin(), x.end(), xq) - x.begin();
    int iy = std::upper_bound(y.begin(), y.end(), yq) - y.begin();
    if (ix <= 0) ix = 1;
    if (iy <= 0) iy = 1;
    if (ix >= x.size()) ix = x.size() - 1;
    if (iy >= y.size()) iy = y.size() - 1;

    const int ix0 = ix - 1;
    const int iy0 = iy - 1;

    const double x1 = x[ix0];
    const double x2 = x[ix];
    const double y1 = y[iy0];
    const double y2 = y[iy];

    const double tx = (x2 == x1) ? 0.0 : (xq - x1) / (x2 - x1);
    const double ty = (y2 == y1) ? 0.0 : (yq - y1) / (y2 - y1);

    const double w11 = (1.0 - tx) * (1.0 - ty);
    const double w21 = tx * (1.0 - ty);
    const double w12 = (1.0 - tx) * ty;
    const double w22 = tx * ty;

    const int idx11 = iy0 + ny * ix0;
    const int idx21 = iy0 + ny * (ix0 + 1);
    const int idx12 = (iy0 + 1) + ny * ix0;
    const int idx22 = (iy0 + 1) + ny * (ix0 + 1);

    for (int s = 0; s < ns; ++s) {
      const int offset = s * plane;
      out(q, s) =
        w11 * z[offset + idx11] +
        w21 * z[offset + idx21] +
        w12 * z[offset + idx12] +
        w22 * z[offset + idx22];
    }
  }

  return out;
}

// [[Rcpp::export]]
NumericVector tmd_harp_cpp(
    NumericVector time_s,
    NumericMatrix hc_re,
    NumericMatrix hc_im,
    NumericMatrix pu,
    NumericMatrix pf,
    NumericVector ph,
    NumericVector omega) {

  const int T = time_s.size();
  const int C = hc_re.nrow();
  const int M = hc_re.ncol();

  if (hc_im.nrow() != C || hc_im.ncol() != M) {
    stop("hc_re and hc_im dimensions must match.");
  }
  if (pu.ncol() != C || pf.ncol() != C || pu.nrow() != T || pf.nrow() != T) {
    stop("pu/pf dimensions must match time and constituent counts.");
  }
  if (ph.size() != C || omega.size() != C) {
    stop("ph/omega lengths must match the number of constituents.");
  }

  NumericVector out;

  if (T == 1) {
    out = NumericVector(M);
    for (int m = 0; m < M; ++m) {
      double acc = 0.0;
      for (int c = 0; c < C; ++c) {
        const double tmp = omega[c] * time_s[0] + ph[c] + pu(0, c);
        acc += pf(0, c) * (hc_re(c, m) * std::cos(tmp) + hc_im(c, m) * std::sin(tmp));
      }
      out[m] = acc;
    }
    return out;
  }

  if (M == 1) {
    out = NumericVector(T);
    for (int t = 0; t < T; ++t) {
      double acc = 0.0;
      for (int c = 0; c < C; ++c) {
        const double tmp = omega[c] * time_s[t] + ph[c] + pu(t, c);
        acc += pf(t, c) * (hc_re(c, 0) * std::cos(tmp) + hc_im(c, 0) * std::sin(tmp));
      }
      out[t] = acc;
    }
    return out;
  }

  if (M == T) {
    out = NumericVector(T);
    for (int t = 0; t < T; ++t) {
      double acc = 0.0;
      for (int c = 0; c < C; ++c) {
        const double tmp = omega[c] * time_s[t] + ph[c] + pu(t, c);
        acc += pf(t, c) * (hc_re(c, t) * std::cos(tmp) + hc_im(c, t) * std::sin(tmp));
      }
      out[t] = acc;
    }
    return out;
  }

  stop("hc must have one column, or one column per input time.");
  return NumericVector::create();
}
