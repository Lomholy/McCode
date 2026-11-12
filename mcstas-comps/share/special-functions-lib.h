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
* Associated Legendre functions. Used by Texture_process.
*
* scaled_associated_legendre is the associated Legendre function scaled to
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
double associated_legendre (int l, int m, double x);
// N_lm (x) for -l <= m <= l
double scaled_associated_legendre (int l, int m, double x);
// scale = sqrt ((l - m)! / (l + m)!) for 0 <= m <= l <= lmax
void legendre_scale_factors (int lmax, double* scale);
// N_lm (x) for 0 <= m <= l <= lmax into values, from the scale of legendre_scale_factors
void all_scaled_associated_legendre (int lmax, double x, double* scale, double* values);

#endif
