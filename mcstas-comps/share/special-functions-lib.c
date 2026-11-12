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
* Associated Legendre functions. See special-functions-lib.h.
*
*******************************************************************************/

#ifndef SPECIAL_FUNCTIONS_LIB_H
#include "special-functions-lib.h"
#endif

#ifndef SPECIAL_FUNCTIONS_LIB_C
#define SPECIAL_FUNCTIONS_LIB_C

double
associated_legendre (int l, int m, double x) {
  double odd_factor, p_lm, p_mm, p_m1m, sin_theta;
  int i, degree;

  if (m < 0 || m > l || fabs (x) > 1.0) {
    printf ("Bad arguments in routine associated_legendre\n");
    printf ("m: %d l: %d x: %f\n", m, l, x);
    printf ("Exiting...\n");
    exit (-1);
  }
  // P_m^m = (-1)^m (2m - 1)!! sin(theta)^m
  p_mm = 1.0;
  if (m > 0) {
    sin_theta = sqrt ((1.0 - x) * (1.0 + x));
    odd_factor = 1.0;
    for (i = 1; i <= m; i++) {
      p_mm *= -odd_factor * sin_theta;
      odd_factor += 2.0;
    }
  }
  if (l == m) {
    return p_mm;
  } else {
    // P_(m+1)^m
    p_m1m = x * (2 * m + 1) * p_mm;
    if (l == (m + 1))
      return p_m1m;
    else {
      // P_l^m for l > m + 1 by recurrence in the degree
      for (degree = m + 2; degree <= l; degree++) {
        p_lm = (x * (2 * degree - 1) * p_m1m - (degree + m - 1) * p_mm) / (degree - m);
        p_mm = p_m1m;
        p_m1m = p_lm;
      }
      return p_lm;
    }
  }
}

double
scaled_associated_legendre (int l, int m, double x) {
  int k;
  int abs_m = abs (m);
  double value, factorial_ratio = 1.0;
  for (k = l + abs_m; k > l - abs_m; k--) {
    factorial_ratio /= k;
  }
  value = sqrt (factorial_ratio) * associated_legendre (l, abs_m, x);
  if (m < 0 && abs_m % 2) {
    value = -value;
  }
  return value;
}

void
legendre_scale_factors (int lmax, double* scale) {
  int l, m, k;
  double factorial_ratio;

  for (l = 0; l <= lmax; l++) {
    for (m = 0; m <= l; m++) {
      factorial_ratio = 1.0;
      for (k = l + m; k > l - m; k--) {
        factorial_ratio /= k;
      }
      scale[l * (lmax + 1) + m] = sqrt (factorial_ratio);
    }
  }
}

void
all_scaled_associated_legendre (int lmax, double x, double* scale, double* values) {
  int l, m, n = lmax + 1;
  double p_mm = 1.0, sin_theta = sqrt ((1.0 - x) * (1.0 + x));

  for (m = 0; m <= lmax; m++) {
    if (m > 0) {
      p_mm *= -(2 * m - 1) * sin_theta;
    }
    values[m * n + m] = p_mm;
    if (m < lmax) {
      values[(m + 1) * n + m] = x * (2 * m + 1) * p_mm;
    }
    for (l = m + 2; l <= lmax; l++) {
      values[l * n + m] = (x * (2 * l - 1) * values[(l - 1) * n + m] - (l + m - 1) * values[(l - 2) * n + m]) / (l - m);
    }
  }
  for (l = 0; l <= lmax; l++) {
    for (m = 0; m <= l; m++) {
      values[l * n + m] *= scale[l * n + m];
    }
  }
}

#endif
