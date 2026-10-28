/*******************************************************************************
*
*  McStas, neutron ray-tracing package
*
* Library: share/sans-expr-lib.c
*
* Small-angle scattering from tinyexpr expressions; see sans-expr-lib.h.
*
*******************************************************************************/

#ifndef SANS_EXPR_LIB_C
#define SANS_EXPR_LIB_C

// Variables of the expressions: q, qx, qy, qz bound to qv[0..3] and p1..p6 bound to pars[0..5]
static const char *const sans_expr_names[SANS_EXPR_NVARS] = { "q", "qx", "qy", "qz", "p1", "p2", "p3", "p4", "p5", "p6" };

static void
sans_expr_addresses (double *qv, double *pars, double **addresses) {
  for (int i = 0; i < 4; i++)
    addresses[i] = &qv[i];
  for (int i = 0; i < 6; i++)
    addresses[4 + i] = &pars[i];
}

te_expr *
sans_expr_compile (const char *expr_str, double *qv, double *pars,
                   const char *component, const char *comp_name, const char *which) {
  double *addresses[SANS_EXPR_NVARS];
  sans_expr_addresses (qv, pars, addresses);
  int err;
  te_expr *out = te_compile_vars (expr_str, sans_expr_names, addresses, SANS_EXPR_NVARS, -1, &err);
  if (!out)
    exit (fprintf (stderr,
                   "%s: %s: failed to parse %s expression \"%s\" "
                   "(tinyexpr error near character %d)\n",
                   component, comp_name, which, expr_str, err));
  return out;
}

/* Returns 1 if expr_str needs none of qx, qy, qz, i.e. depends on |q| only */
int
sans_expr_is_isotropic (const char *expr_str, double *qv, double *pars) {
  double *addresses[SANS_EXPR_NVARS];
  sans_expr_addresses (qv, pars, addresses);
  for (int var = 1; var <= 3; var++)
    if (te_depends_on (expr_str, sans_expr_names, addresses, SANS_EXPR_NVARS, var))
      return 0;
  return 1;
}

/* scale F^2 S with the largest magnitude at the values of |q| from q_low to q_high (log
 * spaced) and, for expressions in qx, qy, qz, in 200 directions spread over the sphere plus
 * the six axis directions; at_q is set to its (q, qx, qy, qz) */
double
sans_expr_largest (te_expr* Fexpr, te_expr* Sexpr, double* qv, double dsdo_scale, int q_only, double q_low, double q_high,
                   int n_q, double* at_q) {
  int n_fib = q_only ? 0 : 200;
  int n_dir = q_only ? 1 : n_fib + 6;
  double axes[6][3] = { { 1, 0, 0 }, { -1, 0, 0 }, { 0, 1, 0 }, { 0, -1, 0 }, { 0, 0, 1 }, { 0, 0, -1 } };
  double golden_angle = PI * (3.0 - sqrt (5.0));
  double largest = 0;
  for (int i = 0; i < n_q; i++) {
    double q = q_low * pow (q_high / q_low, (double)i / (n_q - 1));
    for (int j = 0; j < n_dir; j++) {
      double u[3] = { 1, 0, 0 };
      if (j < n_fib) {
        u[2] = 1.0 - 2.0 * (j + 0.5) / n_fib;
        double ur = sqrt (fmax (0.0, 1.0 - u[2] * u[2]));
        u[0] = ur * cos (j * golden_angle);
        u[1] = ur * sin (j * golden_angle);
      } else if (!q_only) {
        for (int c = 0; c < 3; c++)
          u[c] = axes[j - n_fib][c];
      }
      qv[0] = q;
      for (int c = 0; c < 3; c++)
        qv[c + 1] = q * u[c];
      double F = te_eval (Fexpr);
      double ds = dsdo_scale * F * F * te_eval (Sexpr);
      if (fabs (ds) > fabs (largest)) {
        largest = ds;
        for (int c = 0; c < 4; c++)
          at_q[c] = qv[c];
      }
    }
  }
  return largest;
}

/* Stops the simulation if scale F^2 S at |q| > SANS_EXPR_Q_MAX (up to 4 pi AA^-1, backscattering
 * at 1 AA) exceeds 1e-3 of its largest magnitude at |q| <= SANS_EXPR_Q_MAX: the total
 * cross-section only includes |q| <= SANS_EXPR_Q_MAX */
void
sans_expr_check_large_q (te_expr* Fexpr, te_expr* Sexpr, double* qv, double dsdo_scale, int q_only, const char* component, const char* comp_name) {
  double q_limit = SANS_EXPR_Q_MAX;
  double fraction = 1e-3;
  double at_low[4], at_high[4];
  double low = sans_expr_largest (Fexpr, Sexpr, qv, dsdo_scale, q_only, 1e-4, q_limit, 200, at_low);
  double high = sans_expr_largest (Fexpr, Sexpr, qv, dsdo_scale, q_only, q_limit * 1.0001, 4.0 * PI, 400, at_high);
  if (fabs (high) > fraction * fabs (low))
    exit (fprintf (stderr,
                   "%s: %s: ERROR: scale*F^2*S at |q| > %g AA^-1 reaches %g m^-1 sr^-1 at |q| = %g AA^-1 "
                   "(qx, qy, qz = %g, %g, %g), more than %g of its largest value %g m^-1 sr^-1 at |q| <= %g AA^-1. "
                   "The component is for small-angle scattering, and its total cross-section only includes "
                   "|q| <= %g AA^-1; give expressions that decay at large q, also along qx, qy and qz.\n",
                   component, comp_name, q_limit, high, at_high[0], at_high[1], at_high[2], at_high[3], fraction, low, q_limit, q_limit));
}

/* int q dq <dSigma/dOmega>_phi for k_i = k khat [AA^-1], with
 * q(q,phi) = q^2/2k khat - qperp (cos(phi) e1 + sin(phi) e2), qperp = q sqrt(1 - q^2/4k^2),
 * where e1, e2 are the theta and phi unit vectors of khat (x and y for khat along z). */
double
sans_expr_Idsdw (struct sans_expr_table *t, double k, Coords khat, Coords e1, Coords e2) {
  // Only |q| <= SANS_EXPR_Q_MAX contributes (see sans_expr_check_large_q). Log grids from
  // both ends of the Ewald sphere resolve features near forward and near backward
  // scattering: u = q from 0 to sqrt(2) k, and u = Q = sqrt(4k^2 - q^2) from 0 to sqrt(2) k,
  // with q dq = Q dQ; the backward half is used only where |q| <= SANS_EXPR_Q_MAX reaches it.
  double q_end = fmin (SANS_EXPR_Q_MAX, 2.0 * k);
  double q_split = sqrt (2.0) * k;
  double Qminl = log10 (0.0000001);
  double Qmidl = log10 (q_split);
  int istp = floor ((Qmidl - Qminl) * 300.0 + 0.5);
  double step = (Qmidl - Qminl) / istp;
  double Isq = 0.0;
  for (int half = 0; half < 2; half++) {
    if (half == 1 && q_end <= q_split)
      break;
    double u_lo = half == 0 ? 0.0 : sqrt (fmax (0.0, 4.0 * k * k - q_end * q_end));
    double u_hi = half == 0 ? fmin (q_end, q_split) : q_split;
    double umin = 0.0;
    for (int kkk = 0; kkk <= istp; kkk++) {
      double umax = pow (10.0, Qminl + kkk * step);
      double a = fmax (umin, u_lo);
      double b = fmin (umax, u_hi);
      umin = umax;
      if (b <= a) {
        if (a >= u_hi)
          break;
        continue;
      }
      double u = 0.5 * (a + b);
      double q, qpar, qperp;
      if (half == 0) {
        q = u;
        qpar = q * q / (2.0 * k);
        qperp = q * sqrt (fmax (0.0, 1.0 - q * q / (4.0 * k * k)));
      } else {
        q = sqrt (fmax (0.0, 4.0 * k * k - u * u));
        qpar = 2.0 * k - u * u / (2.0 * k);
        qperp = q * u / (2.0 * k);
      }
      double phi_sum = 0.0;
      for (int jjj = 0; jjj < t->Nphi; jjj++) {
        double phi = 2.0 * PI * (jjj + 0.5) / t->Nphi;
        double c = -qperp * cos (phi);
        double s = -qperp * sin (phi);
        phi_sum += sans_expr_dSigdOmega (t, q, qpar * khat.x + c * e1.x + s * e2.x, qpar * khat.y + c * e1.y + s * e2.y,
                                         qpar * khat.z + c * e1.z + s * e2.z);
      }
      Isq += phi_sum / t->Nphi * u * (b - a);
      if (b >= u_hi)
        break;
    }
  }
  return Isq;
}

/* Direction table node of wavelength ilam [AA], polar index it and azimuth index ip */
double
sans_expr_node (struct sans_expr_table *t, int ilam, int it, int ip) {
  long index = ((long)ilam * t->n_polar + it) * t->n_azim + ip;
  if (!t->filled[index]) {
    double theta = PI * it / (t->n_polar - 1);
    double phi = 2.0 * PI * ip / t->n_azim;
    Coords khat = coords_set (sin (theta) * cos (phi), sin (theta) * sin (phi), cos (theta));
    Coords e1 = coords_set (cos (theta) * cos (phi), cos (theta) * sin (phi), -sin (theta));
    Coords e2 = coords_set (-sin (phi), cos (phi), 0);
    t->I[index] = sans_expr_Idsdw (t, 2.0 * PI / ilam, khat, e1, e2);
    t->filled[index] = 1;
  }
  return t->I[index];
}

/* int q dq <dSigma/dOmega>_phi at wavelength ilam [AA] for k_i along (kx, ky, kz),
 * interpolated linearly in the polar angle and azimuth of k_i */
double
sans_expr_Idsdw_at (struct sans_expr_table *t, int ilam, double kx, double ky, double kz) {
  if (t->q_only)
    return t->iso[ilam];
  double k = sqrt (kx * kx + ky * ky + kz * kz);
  double ft = acos (fmax (-1.0, fmin (1.0, kz / k))) / PI * (t->n_polar - 1);
  int it = (int)floor (ft);
  if (it > t->n_polar - 2)
    it = t->n_polar - 2;
  ft -= it;
  double phi = atan2 (ky, kx);
  if (phi < 0)
    phi += 2.0 * PI;
  double fp = phi / (2.0 * PI) * t->n_azim;
  int ip = (int)floor (fp);
  fp -= ip;
  ip %= t->n_azim;
  int ip2 = (ip + 1) % t->n_azim;
  // Nodes with zero weight are not computed
  double I = 0;
  if ((1 - ft) * (1 - fp) > 0)
    I += (1 - ft) * (1 - fp) * sans_expr_node (t, ilam, it, ip);
  if ((1 - ft) * fp > 0)
    I += (1 - ft) * fp * sans_expr_node (t, ilam, it, ip2);
  if (ft * (1 - fp) > 0)
    I += ft * (1 - fp) * sans_expr_node (t, ilam, it + 1, ip);
  if (ft * fp > 0)
    I += ft * fp * sans_expr_node (t, ilam, it + 1, ip2);
  return I;
}

/* Total scattering cross-section [m^-1] for k_i = (kx, ky, kz) [AA^-1] in the local frame,
 * interpolated linearly in wavelength */
double
sans_expr_sigma_s (struct sans_expr_table *t, double kx, double ky, double kz) {
  double k0 = sqrt (kx * kx + ky * ky + kz * kz);
  double lambda = 2.0 * PI / k0;
  double I;
  if (lambda <= 1.0)
    I = sans_expr_Idsdw_at (t, 1, kx, ky, kz);
  else if (lambda >= 30.0)
    I = sans_expr_Idsdw_at (t, 30, kx, ky, kz);
  else {
    int Ilam = (int)floor (lambda);
    double w = lambda - Ilam;
    I = (1 - w) * sans_expr_Idsdw_at (t, Ilam, kx, ky, kz);
    if (w > 0)
      I += w * sans_expr_Idsdw_at (t, Ilam + 1, kx, ky, kz);
  }
  return 2 * PI * I / (k0 * k0);
}

double
sans_expr_dSigdOmega (struct sans_expr_table *t, double q, double qx, double qy, double qz) {
  t->qv[0] = q;
  t->qv[1] = qx;
  t->qv[2] = qy;
  t->qv[3] = qz;
  double F = te_eval (t->Fexpr);
  return t->dsdo_scale * F * F * te_eval (t->Sexpr) + t->dsdo_background;
}

/* Sets up the table for the compiled expressions, and computes it in full for expressions in
 * q only; for expressions in qx, qy, qz the direction nodes are computed when first needed */
void
sans_expr_table_init (struct sans_expr_table *t, te_expr *Fexpr, te_expr *Sexpr, double *qv, double dsdo_scale, double dsdo_background,
                      int q_only, int n_phi, int n_dir_theta, int n_dir_phi, const char *component, const char *comp_name) {
  if (n_phi < 1)
    exit (fprintf (stderr, "%s: %s: n_phi must be at least 1.\n", component, comp_name));
  t->Fexpr = Fexpr;
  t->Sexpr = Sexpr;
  t->qv = qv;
  t->dsdo_scale = dsdo_scale;
  t->dsdo_background = dsdo_background;
  t->q_only = q_only;
  t->Nphi = q_only ? 1 : n_phi;
  t->n_polar = n_dir_theta;
  t->n_azim = n_dir_phi;
  t->I = NULL;
  t->filled = NULL;
  if (q_only) {
    // Expressions in q only: the integrals do not depend on the direction of k
    t->iso[0] = 0;
    for (int ilam = 1; ilam <= 30; ilam++)
      t->iso[ilam] = sans_expr_Idsdw (t, 2.0 * PI / ilam, coords_set (0, 0, 1), coords_set (1, 0, 0), coords_set (0, 1, 0));
  } else {
    if (n_dir_theta < 2 || n_dir_phi < 1)
      exit (fprintf (stderr, "%s: %s: n_dir_theta must be at least 2 and n_dir_phi at least 1.\n", component, comp_name));
    long n_nodes = 31L * n_dir_theta * n_dir_phi;
    t->I = malloc (n_nodes * sizeof (double));
    t->filled = calloc (n_nodes, 1);
    if (!t->I || !t->filled)
      exit (fprintf (stderr, "%s: %s: could not allocate the Sigma_s table.\n", component, comp_name));
  }
}

void
sans_expr_table_free (struct sans_expr_table *t) {
  free (t->I);
  free (t->filled);
  t->I = NULL;
  t->filled = NULL;
}

#endif
