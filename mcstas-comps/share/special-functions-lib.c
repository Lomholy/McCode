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
gauss_legendre_quadrature (int n, double* nodes, double* weights) {
  int i, degree, iteration;
  double z, z_previous, p_n, p_n1, p_n2, derivative = 1;

  for (i = 0; i < (n + 1) / 2; i++) {
    // Newton iteration for the i-th root of P_n, from an asymptotic first guess
    z = cos (PI * (i + 0.75) / (n + 0.5));
    for (iteration = 0; iteration < 100; iteration++) {
      p_n = 1;
      p_n1 = 0;
      for (degree = 1; degree <= n; degree++) {
        p_n2 = p_n1;
        p_n1 = p_n;
        p_n = ((2 * degree - 1) * z * p_n1 - (degree - 1) * p_n2) / degree;
      }
      derivative = n * (z * p_n - p_n1) / (z * z - 1);
      z_previous = z;
      z = z_previous - p_n / derivative;
      if (fabs (z - z_previous) < 1e-15) {
        break;
      }
    }
    nodes[i] = -z;
    nodes[n - 1 - i] = z;
    weights[i] = weights[n - 1 - i] = 2 / ((1 - z * z) * derivative * derivative);
  }
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

void
add_rotated_harmonics (Rotation rotation, int lmax, int n_points, Coords* points, double* scale, double* legendre_values, double* harmonic_re,
                       double* harmonic_im) {
  int point, l, m, lm;
  double phi, cos_m_phi, sin_m_phi;
  Coords direction;

  for (point = 0; point < n_points; point++) {
    direction = rot_apply (rotation, points[point]);
    all_scaled_associated_legendre (lmax, fmin (1.0, fmax (-1.0, direction.z)), scale, legendre_values);
    phi = atan2 (direction.y, direction.x);
    for (m = 0; m <= lmax; m++) {
      cos_m_phi = cos (m * phi);
      sin_m_phi = sin (m * phi);
      for (l = m; l <= lmax; l++) {
        lm = l * (lmax + 1) + m;
        harmonic_re[lm * n_points + point] += legendre_values[lm] * cos_m_phi;
        harmonic_im[lm * n_points + point] += legendre_values[lm] * sin_m_phi;
      }
    }
  }
}

void
stored_harmonic (int l, int m, int lmax, int n_points, int point, double* harmonic_re, double* harmonic_im, double* re, double* im) {
  int index = (l * (lmax + 1) + abs (m)) * n_points + point;
  double sign = (m < 0 && abs (m) % 2) ? -1 : 1;

  *re = sign * harmonic_re[index];
  *im = (m < 0 ? -sign : sign) * harmonic_im[index];
}

void
harmonic_projection_matrix (int l, int lmax, int n_points, double* weights, double* rotated_re, double* rotated_im, double* reference_re,
                            double* reference_im, double* matrix_re, double* matrix_im) {
  int m, n, point, size = 2 * l + 1;
  double rotated_value_re, rotated_value_im, reference_value_re, reference_value_im, sum_re, sum_im;

  for (m = -l; m <= l; m++) {
    for (n = -l; n <= l; n++) {
      sum_re = 0;
      sum_im = 0;
      for (point = 0; point < n_points; point++) {
        stored_harmonic (l, m, lmax, n_points, point, rotated_re, rotated_im, &rotated_value_re, &rotated_value_im);
        stored_harmonic (l, n, lmax, n_points, point, reference_re, reference_im, &reference_value_re, &reference_value_im);
        sum_re += weights[point] * (rotated_value_re * reference_value_re + rotated_value_im * reference_value_im);
        sum_im += weights[point] * (rotated_value_im * reference_value_re - rotated_value_re * reference_value_im);
      }
      matrix_re[(l + m) * size + l + n] = size / (4 * PI) * sum_re;
      matrix_im[(l + m) * size + l + n] = size / (4 * PI) * sum_im;
    }
  }
}

#endif
