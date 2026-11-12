/*******************************************************************************
*
*  McStas, neutron ray-tracing package
*
* Library: share/special-functions-lib.h
*
* %Identification
* Written by: Victor Laliena, Daniel Lomholt Christensen
* Date: 12/11/2026
* Origin: University of Zaragoza, Niels Bohr Institute @ UCPH
*
* Associated Legendre functions, spherical harmonics and Gauss-Legendre
* quadrature. Used by Texture_process.
*
* scaled_associated_legendre is the associated Legendre function scaled to
*
*   N_lm (x) = sqrt ((l - |m|)! / (l + |m|)!) P_l^|m| (x),  N_l,-m = (-1)^m N_lm,
*
* so that N_lm (cos theta) exp (i m phi) = sqrt (4 pi / (2l + 1)) Y_lm (theta, phi).
* Arrays over 0 <= m <= l <= lmax store (l, m) at l * (lmax + 1) + m. Harmonic
* tables over points store (l, m) at point p at (l * (lmax + 1) + m) * n_points + p.
*
* Usage: %include "special-functions-lib" in the SHARE section of the component.
*
*******************************************************************************/

#ifndef SPECIAL_FUNCTIONS_LIB_H
#define SPECIAL_FUNCTIONS_LIB_H

// Associated Legendre function P_l^m (x), 0 <= m <= l, |x| <= 1
double associated_legendre (int l, int m, double x);
// N_lm (x) for -l <= m <= l
double scaled_associated_legendre (int l, int m, double x);
// scale = sqrt ((l - m)! / (l + m)!) for 0 <= m <= l <= lmax
void legendre_scale_factors (int lmax, double* scale);
// N_lm (x) for 0 <= m <= l <= lmax into values, from the scale of legendre_scale_factors
void all_scaled_associated_legendre (int lmax, double x, double* scale, double* values);
// Gauss-Legendre nodes and weights of order n on [-1, 1]
void gauss_legendre_quadrature (int n, double* nodes, double* weights);
// Adds N_lm (cos theta) exp (i m phi), 0 <= m <= l <= lmax, at the directions rotation * points to (harmonic_re, harmonic_im);
// legendre_values is a work array of (lmax + 1)^2 values
void add_rotated_harmonics (Rotation rotation, int lmax, int n_points, Coords* points, double* scale, double* legendre_values, double* harmonic_re,
                            double* harmonic_im);
// Harmonic (l, m), -l <= m <= l, at one point of a table filled by add_rotated_harmonics
void stored_harmonic (int l, int m, int lmax, int n_points, int point, double* harmonic_re, double* harmonic_im, double* re, double* im);
// Matrix W of degree l, at W[(l + m) * (2l + 1) + l + n], with rotated_lm = sum_n W_mn reference_ln:
// W_mn = (2l + 1) / (4 pi) sum_points weight * rotated_lm * conj (reference_ln)
void harmonic_projection_matrix (int l, int lmax, int n_points, double* weights, double* rotated_re, double* rotated_im, double* reference_re,
                                 double* reference_im, double* matrix_re, double* matrix_im);

#endif
