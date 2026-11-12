/*******************************************************************************
*
*  McStas, neutron ray-tracing package
*
* Library: share/special-functions-lib.c
*
* %Identification
* Written by: Victor Laliena, Daniel Lomholt Christensen
* Date: 12/11/2026
* Origin: University of Zaragoza, Niels Bohr Institute @ UCPH
*
* Associated Legendre functions, spherical harmonics and Gauss-Legendre
* quadrature. See special-functions-lib.h.
*
*******************************************************************************/

#ifndef SPECIAL_FUNCTIONS_LIB_H
#include "special-functions-lib.h"
#endif

#ifndef SPECIAL_FUNCTIONS_LIB_C
#define SPECIAL_FUNCTIONS_LIB_C

double
plgndr (int l, int m, double x) {
  double fact, pll, pmm, pmmp1, somx2;
  int i, ll;

  if (m < 0 || m > l || fabs (x) > 1.0) {
    printf ("Bad arguments in routine plgndr\n");
    printf ("m: %d l: %d x: %f\n", m, l, x);
    printf ("Exiting...\n");
    exit (-1);
  }
  // Compute Pmm
  pmm = 1.0;
  if (m > 0) {
    somx2 = sqrt ((1.0 - x) * (1.0 + x));
    fact = 1.0;
    for (i = 1; i <= m; i++) {
      pmm *= -fact * somx2;
      fact += 2.0;
    }
  }
  if (l == m) {
    return pmm;
  } else {
    // Compute Pmm+1
    pmmp1 = x * (2 * m + 1) * pmm;
    if (l == (m + 1))
      return pmmp1;
    else {
      // Compute Plm, l > m + 1.
      for (ll = m + 2; ll <= l; ll++) {
        pll = (x * (2 * ll - 1) * pmmp1 - (ll + m - 1) * pmm) / (ll - m);
        pmm = pmmp1;
        pmmp1 = pll;
      }
      return pll;
    }
  }
}

// legendre associated functions normalized for computations (see notes)
double
plgndr_mod (int l, int m, double x) {
  int k;
  int ma = abs (m);
  double plm, norm = 1.0;
  for (k = l + ma; k > l - ma; k--) {
    norm /= k;
  }
  plm = sqrt (norm) * plgndr (l, ma, x);
  if (m < 0 && ma % 2) {
    plm = -plm;
  }
  return plm;
}

// Gauss-Legendre nodes x and weights w of order n on [-1, 1]
void
gaussLegendre (int n, double* x, double* w) {
  int i, j, it;
  double z, z1, p1, p2, p3, pp = 1;

  for (i = 0; i < (n + 1) / 2; i++) {
    z = cos (PI * (i + 0.75) / (n + 0.5));
    for (it = 0; it < 100; it++) {
      p1 = 1;
      p2 = 0;
      for (j = 1; j <= n; j++) {
        p3 = p2;
        p2 = p1;
        p1 = ((2 * j - 1) * z * p2 - (j - 1) * p3) / j;
      }
      pp = n * (z * p1 - p2) / (z * z - 1);
      z1 = z;
      z = z1 - p1 / pp;
      if (fabs (z - z1) < 1e-15) {
        break;
      }
    }
    x[i] = -z;
    x[n - 1 - i] = z;
    w[i] = w[n - 1 - i] = 2 / ((1 - z * z) * pp * pp);
  }
}

// sqn[l * (lmax + 1) + m] = sqrt ((l - m)! / (l + m)!) for 0 <= m <= l <= lmax
void
legendreNorms (int lmax, double* sqn) {
  int l, m, k;
  double norm;

  for (l = 0; l <= lmax; l++) {
    for (m = 0; m <= l; m++) {
      norm = 1.0;
      for (k = l + m; k > l - m; k--) {
        norm /= k;
      }
      sqn[l * (lmax + 1) + m] = sqrt (norm);
    }
  }
}

// plgndr_mod (l, m, x) for 0 <= m <= l <= lmax at p[l * (lmax + 1) + m], with sqn = sqrt ((l - m)! / (l + m)!) at the same index
void
allPlgndrMod (int lmax, double x, double* sqn, double* p) {
  int l, m, n = lmax + 1;
  double pmm = 1.0, somx2 = sqrt ((1.0 - x) * (1.0 + x));

  for (m = 0; m <= lmax; m++) {
    if (m > 0) {
      pmm *= -(2 * m - 1) * somx2;
    }
    p[m * n + m] = pmm;
    if (m < lmax) {
      p[(m + 1) * n + m] = x * (2 * m + 1) * pmm;
    }
    for (l = m + 2; l <= lmax; l++) {
      p[l * n + m] = (x * (2 * l - 1) * p[(l - 1) * n + m] - (l + m - 1) * p[(l - 2) * n + m]) / (l - m);
    }
  }
  for (l = 0; l <= lmax; l++) {
    for (m = 0; m <= l; m++) {
      p[l * n + m] *= sqn[l * n + m];
    }
  }
}

// Adds N_lm (R h_q) = plgndr_mod (l, m, cos theta) exp (i m phi) for 0 <= m <= l <= lmax to (Fr, Fi)[(l * (lmax + 1) + m) * nq + q]
void
addHarmonics (Rotation R, int lmax, int nq, Coords* hq, double* sqn, double* p, double* Fr, double* Fi) {
  int q, l, m, k;
  double phi, cm, sm;
  Coords v;

  for (q = 0; q < nq; q++) {
    v = rot_apply (R, hq[q]);
    allPlgndrMod (lmax, fmin (1.0, fmax (-1.0, v.z)), sqn, p);
    phi = atan2 (v.y, v.x);
    for (m = 0; m <= lmax; m++) {
      cm = cos (m * phi);
      sm = sin (m * phi);
      for (l = m; l <= lmax; l++) {
        k = l * (lmax + 1) + m;
        Fr[k * nq + q] += p[k] * cm;
        Fi[k * nq + q] += p[k] * sm;
      }
    }
  }
}

// Value at quadrature point q of a harmonic stored by addHarmonics, with N_l,-m = (-1)^m conj (N_lm)
void
harmonicValue (int l, int m, int lmax, int nq, int q, double* Fr, double* Fi, double* re, double* im) {
  int k = (l * (lmax + 1) + abs (m)) * nq + q;
  double s = (m < 0 && abs (m) % 2) ? -1 : 1;

  *re = s * Fr[k];
  *im = (m < 0 ? -s : s) * Fi[k];
}

// W_mn = (2l + 1) / (4 pi) sum_q w_q F_lm (h_q) conj (G_ln (h_q)) at W[(l + m) * (2l + 1) + l + n]: the matrix with
// F_lm (h) = sum_n W_mn G_ln (h), for G_ln the unrotated harmonics
void
projectHarmonics (int l, int lmax, int nq, double* wq, double* Fr, double* Fi, double* Gr, double* Gi, double* Wr, double* Wi) {
  int m, n, q, d = 2 * l + 1;
  double fr, fi, gr, gi, sr, si;

  for (m = -l; m <= l; m++) {
    for (n = -l; n <= l; n++) {
      sr = 0;
      si = 0;
      for (q = 0; q < nq; q++) {
        harmonicValue (l, m, lmax, nq, q, Fr, Fi, &fr, &fi);
        harmonicValue (l, n, lmax, nq, q, Gr, Gi, &gr, &gi);
        sr += wq[q] * (fr * gr + fi * gi);
        si += wq[q] * (fi * gr - fr * gi);
      }
      Wr[(l + m) * d + l + n] = d / (4 * PI) * sr;
      Wi[(l + m) * d + l + n] = d / (4 * PI) * si;
    }
  }
}

#endif
