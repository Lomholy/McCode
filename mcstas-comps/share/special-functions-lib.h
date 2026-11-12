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
* plgndr_mod is the associated Legendre function scaled to
*
*   N_lm (x) = sqrt ((l - |m|)! / (l + |m|)!) P_l^|m| (x),  N_l,-m = (-1)^m N_lm,
*
* so that N_lm (cos theta) exp (i m phi) = sqrt (4 pi / (2l + 1)) Y_lm (theta, phi).
* Arrays over 0 <= m <= l <= lmax store (l, m) at l * (lmax + 1) + m.
*
* Usage: %include "special-functions-lib" in the SHARE section of the component.
*
*******************************************************************************/

#ifndef SPECIAL_FUNCTIONS_LIB_H
#define SPECIAL_FUNCTIONS_LIB_H

// Associated Legendre function P_l^m (x), 0 <= m <= l, |x| <= 1
double plgndr (int l, int m, double x);
// N_lm (x) for -l <= m <= l
double plgndr_mod (int l, int m, double x);
// sqn = sqrt ((l - m)! / (l + m)!) for 0 <= m <= l <= lmax
void legendreNorms (int lmax, double* sqn);
// N_lm (x) for 0 <= m <= l <= lmax into p, from the sqn of legendreNorms
void allPlgndrMod (int lmax, double x, double* sqn, double* p);
// Gauss-Legendre nodes x and weights w of order n on [-1, 1]
void gaussLegendre (int n, double* x, double* w);
// Adds the harmonics N_lm (R h_q) exp (i m phi), 0 <= m <= l <= lmax, at the nq points R h_q to (Fr, Fi)
void addHarmonics (Rotation R, int lmax, int nq, Coords* hq, double* sqn, double* p, double* Fr, double* Fi);
// Harmonic (l, m), -l <= m <= l, at point q of a table filled by addHarmonics
void harmonicValue (int l, int m, int lmax, int nq, int q, double* Fr, double* Fi, double* re, double* im);
// Matrix W of order l with F_lm = sum_n W_mn G_ln, by quadrature with weights wq
void projectHarmonics (int l, int lmax, int nq, double* wq, double* Fr, double* Fi, double* Gr, double* Gi, double* Wr, double* Wi);

#endif
