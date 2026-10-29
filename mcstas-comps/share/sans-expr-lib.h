/*******************************************************************************
*
*  McStas, neutron ray-tracing package
*
* Library: share/sans-expr-lib.h
*
* %Identification
* Written by: Daniel Lomholt Christensen
* Date: 28/10/2026
* Origin: Niels Bohr Institute @ UCPH
*
* Small-angle scattering from a differential cross-section given as tinyexpr
* expressions,
*
*   dSigma/dOmega(q) = scale * F(q)^2 * S(q) + background,
*
* with F and S in q = |q| and the components qx, qy, qz of q = k_i - k_f [AA^-1],
* and the generic constants p1..p6. Used by Sans_expr and the Union SANS processes.
*
* The total scattering cross-section
*
*   Sigma_s(k) = 2 pi / k^2 int_0^min(q_max, 2k) q dq <dSigma/dOmega(q(q,phi))>_phi
*
* is integrated over the part of the Ewald sphere of k with |q| <= SANS_EXPR_Q_MAX,
* and tabulated for wavelengths 1..30 AA. For expressions in q only it does not
* depend on the direction of k; otherwise it is tabulated also over the direction of
* k (n_dir_theta polar angles and n_dir_phi azimuths), each node computed the first
* time it is needed, and interpolated linearly in wavelength and direction.
*
* Expressions are meant for small-angle scattering: sans_expr_check_large_q stops
* the simulation if scale*F^2*S at |q| > SANS_EXPR_Q_MAX exceeds 1e-3 of its largest
* value below.
*
* Usage: %include "sans-expr-lib" in the SHARE section of the component.
*
*******************************************************************************/

%include "tinyexpr"

#ifndef SANS_EXPR_LIB_H
#define SANS_EXPR_LIB_H

// Number of variables of the expressions: q, qx, qy, qz, p1..p6
#define SANS_EXPR_NVARS 10

// Largest |q| [AA^-1] of the small-angle scattering
#ifndef SANS_EXPR_Q_MAX
#define SANS_EXPR_Q_MAX 0.5
#endif

// Tables of int q dq <dSigma/dOmega>_phi over the Ewald sphere, for wavelengths 1..30 AA:
// iso for expressions in q only, otherwise I over directions of k in the frame of the
// component (n_polar polar angles from 0 to 180 deg, n_azim azimuths), computed when first needed
struct sans_expr_table {
  te_expr* Fexpr;
  te_expr* Sexpr;
  double* qv; // q, qx, qy, qz bound in the expressions
  double dsdo_scale;
  double dsdo_background;
  int Nphi;   // azimuths around k in the integral
  int q_only;
  double iso[31];
  int n_polar;
  int n_azim;
  double* node_I;
  char* filled;
};

te_expr* sans_expr_compile (const char* expr_str, double* qv, double* pars, const char* component, const char* comp_name, const char* which);
int sans_expr_is_isotropic (const char* expr_str, double* qv, double* pars);
double sans_expr_dSigdOmega (struct sans_expr_table* t, double q, double qx, double qy, double qz);
double sans_expr_largest (te_expr* Fexpr, te_expr* Sexpr, double* qv, double dsdo_scale, int q_only, double q_low, double q_high, int n_q,
                          double* at_q);
void sans_expr_check_large_q (te_expr* Fexpr, te_expr* Sexpr, double* qv, double dsdo_scale, int q_only, const char* component,
                              const char* comp_name);
void sans_expr_table_init (struct sans_expr_table* t, te_expr* Fexpr, te_expr* Sexpr, double* qv, double dsdo_scale, double dsdo_background,
                           int q_only, int n_phi, int n_dir_theta, int n_dir_phi, const char* component, const char* comp_name);
void sans_expr_table_free (struct sans_expr_table* t);
double sans_expr_Idsdw (struct sans_expr_table* t, double k, Coords khat, Coords e1, Coords e2);
double sans_expr_node (struct sans_expr_table* t, int ilam, int it, int ip);
double sans_expr_Idsdw_at (struct sans_expr_table* t, int ilam, double kx, double ky, double kz);
double sans_expr_sigma_s (struct sans_expr_table* t, double kx, double ky, double kz);

#endif
