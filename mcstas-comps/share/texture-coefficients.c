/*******************************************************************************
*
*  McStas, neutron ray-tracing package
*
* Program: share/texture-coefficients.c
*
* %Identification
* Written by: Daniel Lomholt Christensen
* Date: 12/11/2026
* Origin: Niels Bohr Institute @ UCPH
*
* Writes the Fourier coefficients of an orientation distribution function (ODF) to the
* file that Texture_process reads as fcoef_fn, with the five columns l m n Real(Clmn) Imag(Clmn).
* The ODF is given by a list of orientations, either an orientation list such as an EBSD
* map, or orientations drawn from three tinyexpr expressions.
*
* The orientations rotate the Cartesian crystal axes into the axes of the component.
* The crystal axes are those of Texture_process: by default those of NCrystal, x || a and
* b in the xy plane; mtex_axes = 1 takes those of MTEX, z || c and x || a*; ax ... cz place
* the lattice vectors along given directions in the component. Give the same lattice,
* mtex_axes and ax ... cz as to Texture_process.
*
* The ODF is the orientations smoothed by a Gaussian kernel on the rotation group, with the
* Fourier coefficients kappa_l = exp(-l(l+1) sigma^2 / 2), averaged over the proper rotations
* of the lattice. Each orientation has the same weight, and the computation time grows as
* the number of orientations times lmax^4.
*
* ebsd: each row starts with the Bunge Euler angles phi1 Phi phi2 in radians of one
* orientation, crystal to sample; further columns are ignored, and rows that do not start
* with three numbers are skipped.
*
* rot_x_expr, rot_y_expr, rot_z_expr: the angles in deg of rotations of the crystal axes
* about the x, y and z axes of the component, composed as ROTATED (rx, ry, rz) composes
* them: g = Rx (rx) Ry (ry) Rz (rz). The expressions use the random numbers u1, u2, u3 in
* [0, 1), drawn anew for each of the n_orientations orientations from the generator seeded
* by seed, the constants p1 ... p6, and pi. The noise of the drawn texture falls as
* 1/sqrt(n_orientations) and grows as orientation_sigma decreases, about as
* orientation_sigma^-3 in the number of orientations needed.
*
* Build:  cc -O2 -o texture-coefficients texture-coefficients.c -lm
* (in this folder, which has tinyexpr and special-functions-lib)
*
* Usage:  texture-coefficients out=FILE a=A [name=value ...]
*
*   out                  file to write
*   a, b, c              lattice parameters [AA]; b and c default to the one before
*   alpha, beta, gamma   lattice angles [deg]; alpha defaults to 90, beta and gamma to the one before
*   mtex_axes            0 (default) or 1, as in Texture_process
*   ax, ay, az, bx, by, bz, cx, cy, cz   directions of a, b, c in the component, as in Texture_process
*   orientation_sigma    standard deviation of the kernel [deg], default 10
*   lmax                 cut-off of the series; negative (default) gives the l where the kernel has fallen to 1e-3
*   ebsd                 orientation list
*   rot_x_expr, rot_y_expr, rot_z_expr   rotation expressions, instead of ebsd
*   n_orientations       number of orientations drawn from the expressions, default 10000
*   p1 ... p6            constants of the expressions, default 0
*   seed                 seed of the random numbers of the expressions, default 1
*
* Example, a fibre texture with a tilted from x by a Gaussian of 5 deg and turned about x:
*   texture-coefficients out=fibre.txt a=4.04958 mtex_axes=1 p1=5 \
*     rot_x_expr="360*u3" rot_y_expr="p1*sqrt(-2*log(1-u1))*cos(2*pi*u2)" rot_z_expr="0"
*
*******************************************************************************/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "tinyexpr.h"
#include "tinyexpr.c"
#include "special-functions-lib.c"

typedef struct {
  double x, y, z;
} Vector;

typedef double RotationMatrix[3][3];

Vector
vector_set (double x, double y, double z) {
  Vector v;
  v.x = x;
  v.y = y;
  v.z = z;
  return v;
}

double
dot_product (Vector a, Vector b) {
  return a.x * b.x + a.y * b.y + a.z * b.z;
}

Vector
cross_product (Vector b, Vector c) {
  return vector_set (b.y * c.z - c.y * b.z, b.z * c.x - c.z * b.x, b.x * c.y - c.x * b.y);
}

Vector
vector_scale (Vector v, double scale) {
  return vector_set (v.x * scale, v.y * scale, v.z * scale);
}

double
vector_length (Vector v) {
  return sqrt (v.x * v.x + v.y * v.y + v.z * v.z);
}

Vector
rotate_vector (RotationMatrix R, Vector v) {
  return vector_set (R[0][0] * v.x + R[0][1] * v.y + R[0][2] * v.z, R[1][0] * v.x + R[1][1] * v.y + R[1][2] * v.z,
                     R[2][0] * v.x + R[2][1] * v.y + R[2][2] * v.z);
}

// Right-handed orthonormal basis: basis[0] || v1, basis[1] in the plane of v1 and v2
int
orthonormal_basis (Vector v1, Vector v2, Vector* basis) {
  Vector normal = cross_product (v1, v2);
  if (vector_length (v1) == 0 || vector_length (normal) < 1e-12 * vector_length (v1) * vector_length (v2)) {
    return 0;
  }
  basis[0] = vector_scale (v1, 1.0 / vector_length (v1));
  basis[2] = vector_scale (normal, 1.0 / vector_length (normal));
  basis[1] = cross_product (basis[2], basis[0]);
  return 1;
}

// Rotation R that takes crystal direction c1 onto component direction l1, and c2 into the plane of l1 and l2
void
rotation_between_directions (Vector c1, Vector c2, Vector l1, Vector l2, RotationMatrix R) {
  Vector crystal_basis[3], component_basis[3];
  double e[3][3], f[3][3];
  int i, j, k;

  if (!orthonormal_basis (c1, c2, crystal_basis) || !orthonormal_basis (l1, l2, component_basis)) {
    fprintf (stderr, "texture-coefficients: the directions of the two lattice vectors must be nonzero and not parallel.\n");
    exit (1);
  }
  for (k = 0; k < 3; k++) {
    e[k][0] = crystal_basis[k].x;
    e[k][1] = crystal_basis[k].y;
    e[k][2] = crystal_basis[k].z;
    f[k][0] = component_basis[k].x;
    f[k][1] = component_basis[k].y;
    f[k][2] = component_basis[k].z;
  }
  // R = sum_k f_k e_k^T
  for (i = 0; i < 3; i++) {
    for (j = 0; j < 3; j++) {
      R[i][j] = f[0][i] * e[0][j] + f[1][i] * e[1][j] + f[2][i] * e[2][j];
    }
  }
}

// Lattice vectors along the three directions, with the lengths of lattice. The angles between the directions
// must be those between the lattice vectors, and the directions must be right-handed.
void
lattice_along_directions (Vector* lattice, Vector* directions) {
  double tolerance = 1e-4; // rad
  double angle_lattice, angle_directions;
  Vector unit[3];
  int i, i_next, i_last;
  char name[3] = { 'a', 'b', 'c' };

  for (i = 0; i < 3; i++) {
    unit[i] = vector_scale (directions[i], 1.0 / vector_length (directions[i]));
  }
  for (i = 0; i < 3; i++) {
    i_next = (i + 1) % 3;
    i_last = (i + 2) % 3;
    angle_lattice = acos (dot_product (lattice[i_next], lattice[i_last]) / (vector_length (lattice[i_next]) * vector_length (lattice[i_last])));
    angle_directions = acos (fmax (-1.0, fmin (1.0, dot_product (unit[i_next], unit[i_last]))));
    if (fabs (angle_directions - angle_lattice) > tolerance) {
      fprintf (stderr, "texture-coefficients: the directions of %c and %c are %g deg apart, but %g deg in the lattice.\n", name[i_next], name[i_last],
               angle_directions * 180 / M_PI, angle_lattice * 180 / M_PI);
      exit (1);
    }
  }
  if (dot_product (unit[0], cross_product (unit[1], unit[2])) <= 0) {
    fprintf (stderr, "texture-coefficients: the directions of a, b and c must be right-handed.\n");
    exit (1);
  }
  for (i = 0; i < 3; i++) {
    lattice[i] = vector_scale (unit[i], vector_length (lattice[i]));
  }
}

// Lattice vectors in the Cartesian crystal frame of Texture_process
void
crystal_lattice (double a, double b, double c, double alpha, double beta, double gamma, int mtex_axes, Vector* directions, Vector* lattice) {
  double cx, cy;
  int i, i_first = -1, i_second = -1, n_given = 0;
  RotationMatrix R;

  // NCrystal: x || a and b in the xy plane, so z || c*
  cx = c * cos (beta * (M_PI / 180));
  cy = c * (cos (alpha * (M_PI / 180)) - cos (beta * (M_PI / 180)) * cos (gamma * (M_PI / 180))) / sin (gamma * (M_PI / 180));
  lattice[0] = vector_set (a, 0, 0);
  lattice[1] = vector_set (b * cos (gamma * (M_PI / 180)), b * sin (gamma * (M_PI / 180)), 0);
  lattice[2] = vector_set (cx, cy, sqrt (c * c - cx * cx - cy * cy));

  for (i = 0; i < 3; i++) {
    if (vector_length (directions[i]) > 0) {
      n_given++;
      if (i_first < 0) {
        i_first = i;
      } else if (i_second < 0) {
        i_second = i;
      }
    }
  }
  if (n_given == 1) {
    fprintf (stderr, "texture-coefficients: give the directions of two or three of the lattice vectors a, b and c, or none.\n");
    exit (1);
  }
  if (n_given && mtex_axes) {
    fprintf (stderr, "texture-coefficients: the directions of the lattice vectors cannot be combined with mtex_axes = 1.\n");
    exit (1);
  }
  if (n_given == 3) {
    lattice_along_directions (lattice, directions);
  } else if (mtex_axes || n_given == 2) {
    if (mtex_axes) {
      rotation_between_directions (lattice[2], cross_product (lattice[1], lattice[2]), vector_set (0, 0, 1), vector_set (1, 0, 0), R);
    } else {
      rotation_between_directions (lattice[i_first], lattice[i_second], directions[i_first], directions[i_second], R);
    }
    for (i = 0; i < 3; i++) {
      lattice[i] = rotate_vector (R, lattice[i]);
    }
  }
}

// Proper rotations of the Cartesian crystal frame that map the lattice onto itself
int
lattice_rotations (Vector* lattice, RotationMatrix* rotations) {
  double A[3][3], A_inverse[3][3], AM[3][3];
  double volume, deviation, determinant;
  int M[3][3], code, i, j, k, n_rotations = 0;
  Vector row;

  volume = dot_product (lattice[0], cross_product (lattice[1], lattice[2]));
  for (j = 0; j < 3; j++) {
    A[0][j] = lattice[j].x;
    A[1][j] = lattice[j].y;
    A[2][j] = lattice[j].z;
    row = vector_scale (cross_product (lattice[(j + 1) % 3], lattice[(j + 2) % 3]), 1.0 / volume);
    A_inverse[j][0] = row.x;
    A_inverse[j][1] = row.y;
    A_inverse[j][2] = row.z;
  }

  // integer matrices M with entries -1, 0, 1 in the lattice basis, R = A M A^-1
  for (code = 0; code < 19683; code++) {
    k = code;
    for (i = 0; i < 3; i++) {
      for (j = 0; j < 3; j++) {
        M[i][j] = k % 3 - 1;
        k /= 3;
      }
    }
    for (i = 0; i < 3; i++) {
      for (j = 0; j < 3; j++) {
        AM[i][j] = A[i][0] * M[0][j] + A[i][1] * M[1][j] + A[i][2] * M[2][j];
      }
    }
    for (i = 0; i < 3; i++) {
      for (j = 0; j < 3; j++) {
        rotations[n_rotations][i][j] = AM[i][0] * A_inverse[0][j] + AM[i][1] * A_inverse[1][j] + AM[i][2] * A_inverse[2][j];
      }
    }
    deviation = 0;
    for (i = 0; i < 3; i++) {
      for (j = 0; j < 3; j++) {
        deviation += fabs (rotations[n_rotations][0][i] * rotations[n_rotations][0][j] + rotations[n_rotations][1][i] * rotations[n_rotations][1][j]
                           + rotations[n_rotations][2][i] * rotations[n_rotations][2][j] - (i == j));
      }
    }
    determinant = rotations[n_rotations][0][0] * (rotations[n_rotations][1][1] * rotations[n_rotations][2][2] - rotations[n_rotations][1][2] * rotations[n_rotations][2][1])
                  - rotations[n_rotations][0][1] * (rotations[n_rotations][1][0] * rotations[n_rotations][2][2] - rotations[n_rotations][1][2] * rotations[n_rotations][2][0])
                  + rotations[n_rotations][0][2] * (rotations[n_rotations][1][0] * rotations[n_rotations][2][1] - rotations[n_rotations][1][1] * rotations[n_rotations][2][0]);
    if (deviation < 1e-6 && determinant > 0) {
      n_rotations++;
    }
  }
  return n_rotations;
}

// Gauss-Legendre nodes and weights of order n on [-1, 1]
void
gauss_legendre_quadrature (int n, double* nodes, double* weights) {
  int i, degree, iteration;
  double z, z_previous, p_n, p_n1, p_n2, derivative = 1;

  for (i = 0; i < (n + 1) / 2; i++) {
    // Newton iteration for the i-th root of P_n, from an asymptotic first guess
    z = cos (M_PI * (i + 0.75) / (n + 0.5));
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

// Adds N_lm (cos theta) exp (i m phi), 0 <= m <= l <= lmax, at the directions rotation * points to (harmonic_re, harmonic_im),
// with (l, m) at point p at (l * (lmax + 1) + m) * n_points + p. legendre_values is a work array of (lmax + 1)^2 values.
void
add_rotated_harmonics (RotationMatrix rotation, int lmax, int n_points, Vector* points, double* scale, double* legendre_values, double* harmonic_re,
                       double* harmonic_im) {
  int point, l, m, lm;
  double phi, cos_m_phi, sin_m_phi;
  Vector direction;

  for (point = 0; point < n_points; point++) {
    direction = rotate_vector (rotation, points[point]);
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

// Harmonic (l, m), -l <= m <= l, at one point of a table filled by add_rotated_harmonics, with N_l,-m = (-1)^m conj (N_lm)
void
stored_harmonic (int l, int m, int lmax, int n_points, int point, double* harmonic_re, double* harmonic_im, double* re, double* im) {
  int index = (l * (lmax + 1) + abs (m)) * n_points + point;
  double sign = (m < 0 && abs (m) % 2) ? -1 : 1;

  *re = sign * harmonic_re[index];
  *im = (m < 0 ? -sign : sign) * harmonic_im[index];
}

// Matrix W of degree l, at W[(l + m) * (2l + 1) + l + n], with rotated_lm = sum_n W_mn reference_ln:
// W_mn = (2l + 1) / (4 pi) sum_points weight * rotated_lm * conj (reference_ln)
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
      matrix_re[(l + m) * size + l + n] = size / (4 * M_PI) * sum_re;
      matrix_im[(l + m) * size + l + n] = size / (4 * M_PI) * sum_im;
    }
  }
}

// Writes the Fourier coefficients of the ODF of the orientations, smoothed by the kernel exp (-l (l + 1) sigma^2 / 2)
// and averaged over the lattice rotations: C_lmn = (2l + 1) kappa_l [avg_i W (g_i) avg_s W (s)]_mn, where
// N_lm (g h) = sum_n W (g)_mn N_ln (h). Each row is l, n, m, Re, Im with Re + i Im = s conj (C_lmn), and s = -1 when the
// negative ones of m and n sum to an odd number: the inverse of readFourierCoef of Texture_process.
void
write_orientation_coefficients (RotationMatrix* orientations, long n_orientations, double sigma, int lmax, Vector* lattice, char* source,
                                char* out_file) {
  int n_nodes, n_phi, n_points, n_harmonics, n_symmetry, i, j, l, m, n, k, size, sign;
  double *nodes, *node_weights, *weights, *scale, *legendre_values;
  double *orientation_re, *orientation_im, *reference_re, *reference_im, *symmetry_re, *symmetry_im;
  double *orientation_matrix_re, *orientation_matrix_im, *symmetry_matrix_re, *symmetry_matrix_im;
  double kappa, sin_theta, product_re, product_im;
  Vector* points;
  RotationMatrix identity = { { 1, 0, 0 }, { 0, 1, 0 }, { 0, 0, 1 } }, symmetry[48];
  FILE* file;

  if (lmax < 0) {
    lmax = (int)ceil (sqrt (2 * log (1000.0)) / sigma);
  }
  n_nodes = lmax + 2;
  n_phi = 2 * lmax + 4;
  n_points = n_nodes * n_phi;
  n_harmonics = (lmax + 1) * (lmax + 1) * n_points;

  nodes = malloc (n_nodes * sizeof (double));
  node_weights = malloc (n_nodes * sizeof (double));
  weights = malloc (n_points * sizeof (double));
  points = malloc (n_points * sizeof (Vector));
  scale = malloc ((lmax + 1) * (lmax + 1) * sizeof (double));
  legendre_values = malloc ((lmax + 1) * (lmax + 1) * sizeof (double));
  orientation_re = calloc (n_harmonics, sizeof (double));
  orientation_im = calloc (n_harmonics, sizeof (double));
  reference_re = calloc (n_harmonics, sizeof (double));
  reference_im = calloc (n_harmonics, sizeof (double));
  symmetry_re = calloc (n_harmonics, sizeof (double));
  symmetry_im = calloc (n_harmonics, sizeof (double));
  orientation_matrix_re = malloc ((2 * lmax + 1) * (2 * lmax + 1) * sizeof (double));
  orientation_matrix_im = malloc ((2 * lmax + 1) * (2 * lmax + 1) * sizeof (double));
  symmetry_matrix_re = malloc ((2 * lmax + 1) * (2 * lmax + 1) * sizeof (double));
  symmetry_matrix_im = malloc ((2 * lmax + 1) * (2 * lmax + 1) * sizeof (double));
  if (!nodes || !node_weights || !weights || !points || !scale || !legendre_values || !orientation_re || !orientation_im || !reference_re
      || !reference_im || !symmetry_re || !symmetry_im || !orientation_matrix_re || !orientation_matrix_im || !symmetry_matrix_re
      || !symmetry_matrix_im) {
    fprintf (stderr, "texture-coefficients: out of memory for lmax = %d.\n", lmax);
    exit (1);
  }
  file = fopen (out_file, "w");
  if (!file) {
    fprintf (stderr, "texture-coefficients: cannot write %s.\n", out_file);
    exit (1);
  }

  // quadrature points on the sphere, exact for the products of two harmonics of degree <= lmax
  gauss_legendre_quadrature (n_nodes, nodes, node_weights);
  for (i = 0; i < n_nodes; i++) {
    sin_theta = sqrt (1 - nodes[i] * nodes[i]);
    for (j = 0; j < n_phi; j++) {
      points[i * n_phi + j] = vector_set (sin_theta * cos (2 * M_PI * j / n_phi), sin_theta * sin (2 * M_PI * j / n_phi), nodes[i]);
      weights[i * n_phi + j] = node_weights[i] * 2 * M_PI / n_phi;
    }
  }
  legendre_scale_factors (lmax, scale);

  for (i = 0; i < n_orientations; i++) {
    add_rotated_harmonics (orientations[i], lmax, n_points, points, scale, legendre_values, orientation_re, orientation_im);
  }
  n_symmetry = lattice_rotations (lattice, symmetry);
  for (i = 0; i < n_symmetry; i++) {
    add_rotated_harmonics (symmetry[i], lmax, n_points, points, scale, legendre_values, symmetry_re, symmetry_im);
  }
  add_rotated_harmonics (identity, lmax, n_points, points, scale, legendre_values, reference_re, reference_im);
  for (k = 0; k < n_harmonics; k++) {
    orientation_re[k] /= n_orientations;
    orientation_im[k] /= n_orientations;
    symmetry_re[k] /= n_symmetry;
    symmetry_im[k] /= n_symmetry;
  }

  fprintf (file, "# Fourier coefficients of the ODF of %s: %ld orientations, %d lattice rotations, sigma = %g deg, lmax = %d\n", source,
           n_orientations, n_symmetry, sigma * 180 / M_PI, lmax);
  fprintf (file, "# written by texture-coefficients for fcoef_fn of Texture_process\n");
  fprintf (file, "# l m n Real(Clmn) Imag(Clmn)\n");
  for (l = 0; l <= lmax; l++) {
    size = 2 * l + 1;
    kappa = exp (-l * (l + 1) * sigma * sigma / 2);
    harmonic_projection_matrix (l, lmax, n_points, weights, orientation_re, orientation_im, reference_re, reference_im, orientation_matrix_re,
                                orientation_matrix_im);
    harmonic_projection_matrix (l, lmax, n_points, weights, symmetry_re, symmetry_im, reference_re, reference_im, symmetry_matrix_re, symmetry_matrix_im);
    for (n = -l; n <= l; n++) {
      for (m = -l; m <= l; m++) {
        product_re = 0;
        product_im = 0;
        for (k = 0; k < size; k++) {
          product_re += orientation_matrix_re[(l + m) * size + k] * symmetry_matrix_re[k * size + l + n]
                        - orientation_matrix_im[(l + m) * size + k] * symmetry_matrix_im[k * size + l + n];
          product_im += orientation_matrix_re[(l + m) * size + k] * symmetry_matrix_im[k * size + l + n]
                        + orientation_matrix_im[(l + m) * size + k] * symmetry_matrix_re[k * size + l + n];
        }
        sign = ((m < 0 ? -m : 0) + (n < 0 ? -n : 0)) % 2 ? -1 : 1;
        fprintf (file, "%d %d %d %.17g %.17g\n", l, n, m, sign * (size * kappa * product_re), -sign * (size * kappa * product_im));
      }
    }
  }
  fclose (file);
  printf ("texture-coefficients: %s: %ld orientations, %d lattice rotations, sigma = %g deg, lmax = %d, written to %s\n", source, n_orientations,
          n_symmetry, sigma * 180 / M_PI, lmax, out_file);

  free (nodes);
  free (node_weights);
  free (weights);
  free (points);
  free (scale);
  free (legendre_values);
  free (orientation_re);
  free (orientation_im);
  free (reference_re);
  free (reference_im);
  free (symmetry_re);
  free (symmetry_im);
  free (orientation_matrix_re);
  free (orientation_matrix_im);
  free (symmetry_matrix_re);
  free (symmetry_matrix_im);
}

// Orientations of the rows of filename, Bunge Euler angles g = Rz (phi1) Rx (Phi) Rz (phi2)
long
read_ebsd (char* filename, RotationMatrix** orientations) {
  char line[4096];
  char* end;
  double angle[3], c1, s1, c2, s2, C, S;
  long n_orientations = 0, capacity = 1024;
  int i;
  RotationMatrix* g;
  FILE* file = fopen (filename, "r");

  *orientations = malloc (capacity * sizeof (RotationMatrix));
  if (!file || !*orientations) {
    fprintf (stderr, "texture-coefficients: cannot read %s.\n", filename);
    exit (1);
  }
  while (fgets (line, sizeof (line), file)) {
    end = line;
    for (i = 0; i < 3; i++) {
      char* start = end + strspn (end, " \t,;");
      angle[i] = strtod (start, &end);
      if (end == start) {
        break;
      }
    }
    if (i < 3) {
      continue;
    }
    if (n_orientations == capacity) {
      capacity *= 2;
      *orientations = realloc (*orientations, capacity * sizeof (RotationMatrix));
      if (!*orientations) {
        fprintf (stderr, "texture-coefficients: out of memory reading %s.\n", filename);
        exit (1);
      }
    }
    c1 = cos (angle[0]);
    s1 = sin (angle[0]);
    C = cos (angle[1]);
    S = sin (angle[1]);
    c2 = cos (angle[2]);
    s2 = sin (angle[2]);
    g = *orientations + n_orientations;
    (*g)[0][0] = c1 * c2 - s1 * s2 * C;
    (*g)[0][1] = -c1 * s2 - s1 * c2 * C;
    (*g)[0][2] = s1 * S;
    (*g)[1][0] = s1 * c2 + c1 * s2 * C;
    (*g)[1][1] = -s1 * s2 + c1 * c2 * C;
    (*g)[1][2] = -c1 * S;
    (*g)[2][0] = s2 * S;
    (*g)[2][1] = c2 * S;
    (*g)[2][2] = C;
    n_orientations++;
  }
  fclose (file);
  if (n_orientations < 1) {
    fprintf (stderr, "texture-coefficients: %s must have rows starting with the Euler angles phi1 Phi phi2.\n", filename);
    exit (1);
  }
  return n_orientations;
}

// Uniform random number in [0, 1), splitmix64
double
random_uniform (unsigned long long* state) {
  unsigned long long z = (*state += 0x9E3779B97F4A7C15ULL);
  z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
  z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
  z ^= z >> 31;
  return (z >> 11) * (1.0 / 9007199254740992.0);
}

// n_orientations orientations from the expressions of the angles (deg) about x, y and z, composed as ROTATED:
// g = Rx (rx) Ry (ry) Rz (rz), active, crystal axes to component axes
RotationMatrix*
draw_expression_orientations (char** expressions, long n_orientations, double* constants, unsigned long long seed) {
  const char* names[9] = { "u1", "u2", "u3", "p1", "p2", "p3", "p4", "p5", "p6" };
  double values[9];
  double* addresses[9];
  double angle[3], cx, sx, cy, sy, cz, sz;
  te_expr* compiled[3];
  int i_value, axis, error;
  long i_orientation;
  RotationMatrix* orientations;

  for (i_value = 0; i_value < 9; i_value++) {
    addresses[i_value] = values + i_value;
  }
  for (i_value = 0; i_value < 6; i_value++) {
    values[3 + i_value] = constants[i_value];
  }
  for (axis = 0; axis < 3; axis++) {
    compiled[axis] = te_compile_vars (expressions[axis], names, addresses, 9, -1, &error);
    if (!compiled[axis]) {
      fprintf (stderr, "texture-coefficients: cannot parse the rotation expression \"%s\" at character %d.\n", expressions[axis], error);
      exit (1);
    }
  }
  if (n_orientations < 1) {
    fprintf (stderr, "texture-coefficients: n_orientations must be positive.\n");
    exit (1);
  }
  orientations = malloc (n_orientations * sizeof (RotationMatrix));
  if (!orientations) {
    fprintf (stderr, "texture-coefficients: out of memory for %ld orientations.\n", n_orientations);
    exit (1);
  }

  for (i_orientation = 0; i_orientation < n_orientations; i_orientation++) {
    values[0] = random_uniform (&seed);
    values[1] = random_uniform (&seed);
    values[2] = random_uniform (&seed);
    for (axis = 0; axis < 3; axis++) {
      angle[axis] = te_eval (compiled[axis]) * (M_PI / 180);
      if (!isfinite (angle[axis])) {
        fprintf (stderr, "texture-coefficients: the rotation expression \"%s\" is not finite for u1 = %g, u2 = %g, u3 = %g.\n", expressions[axis],
                 values[0], values[1], values[2]);
        exit (1);
      }
    }
    // the active rotation, the transpose of the coordinate transformation of ROTATED (rot_set_rotation)
    cx = cos (angle[0]);
    sx = sin (angle[0]);
    cy = cos (angle[1]);
    sy = sin (angle[1]);
    cz = cos (angle[2]);
    sz = sin (angle[2]);
    orientations[i_orientation][0][0] = cy * cz;
    orientations[i_orientation][1][0] = sx * sy * cz + cx * sz;
    orientations[i_orientation][2][0] = sx * sz - cx * sy * cz;
    orientations[i_orientation][0][1] = -cy * sz;
    orientations[i_orientation][1][1] = cx * cz - sx * sy * sz;
    orientations[i_orientation][2][1] = sx * cz + cx * sy * sz;
    orientations[i_orientation][0][2] = sy;
    orientations[i_orientation][1][2] = -sx * cy;
    orientations[i_orientation][2][2] = cx * cy;
  }

  for (axis = 0; axis < 3; axis++) {
    te_free (compiled[axis]);
  }
  return orientations;
}

int
main (int argc, char** argv) {
  char *out_file = NULL, *ebsd_file = NULL, *expressions[3] = { NULL, NULL, NULL };
  double a = 0, b = 0, c = 0, alpha = 0, beta = 0, gamma = 0, sigma_deg = 10, constants[6] = { 0 };
  double direction_values[9] = { 0 };
  const char* direction_names[9] = { "ax", "ay", "az", "bx", "by", "bz", "cx", "cy", "cz" };
  const char* constant_names[6] = { "p1", "p2", "p3", "p4", "p5", "p6" };
  int mtex_axes = 0, lmax = -1, i_arg, i, n_expressions;
  long n_orientations = 10000;
  unsigned long long seed = 1;
  char *name, *value;
  Vector directions[3], lattice[3];
  RotationMatrix* orientations;

  for (i_arg = 1; i_arg < argc; i_arg++) {
    name = argv[i_arg];
    value = strchr (name, '=');
    if (!value) {
      fprintf (stderr, "texture-coefficients: arguments are name=value, not \"%s\". See the top of texture-coefficients.c.\n", name);
      return 1;
    }
    *value++ = '\0';
    if (!strcmp (name, "out")) {
      out_file = value;
    } else if (!strcmp (name, "ebsd")) {
      ebsd_file = value;
    } else if (!strcmp (name, "rot_x_expr")) {
      expressions[0] = value;
    } else if (!strcmp (name, "rot_y_expr")) {
      expressions[1] = value;
    } else if (!strcmp (name, "rot_z_expr")) {
      expressions[2] = value;
    } else if (!strcmp (name, "a")) {
      a = atof (value);
    } else if (!strcmp (name, "b")) {
      b = atof (value);
    } else if (!strcmp (name, "c")) {
      c = atof (value);
    } else if (!strcmp (name, "alpha")) {
      alpha = atof (value);
    } else if (!strcmp (name, "beta")) {
      beta = atof (value);
    } else if (!strcmp (name, "gamma")) {
      gamma = atof (value);
    } else if (!strcmp (name, "mtex_axes")) {
      mtex_axes = atoi (value);
    } else if (!strcmp (name, "orientation_sigma")) {
      sigma_deg = atof (value);
    } else if (!strcmp (name, "lmax")) {
      lmax = atoi (value);
    } else if (!strcmp (name, "n_orientations")) {
      n_orientations = atol (value);
    } else if (!strcmp (name, "seed")) {
      seed = strtoull (value, NULL, 10);
    } else {
      for (i = 0; i < 9 && strcmp (name, direction_names[i]); i++)
        ;
      if (i < 9) {
        direction_values[i] = atof (value);
        continue;
      }
      for (i = 0; i < 6 && strcmp (name, constant_names[i]); i++)
        ;
      if (i == 6) {
        fprintf (stderr, "texture-coefficients: unknown parameter %s. See the top of texture-coefficients.c.\n", name);
        return 1;
      }
      constants[i] = atof (value);
    }
  }

  n_expressions = (expressions[0] != NULL) + (expressions[1] != NULL) + (expressions[2] != NULL);
  if (!out_file || a <= 0) {
    fprintf (stderr, "texture-coefficients: give out=FILE and the lattice parameter a. See the top of texture-coefficients.c.\n");
    return 1;
  }
  if (n_expressions != 0 && n_expressions != 3) {
    fprintf (stderr, "texture-coefficients: give all three of rot_x_expr, rot_y_expr and rot_z_expr, or none.\n");
    return 1;
  }
  if ((ebsd_file != NULL) + (n_expressions == 3) != 1) {
    fprintf (stderr, "texture-coefficients: give the orientations in exactly one of ebsd and the rotation expressions.\n");
    return 1;
  }
  if (sigma_deg <= 0) {
    fprintf (stderr, "texture-coefficients: orientation_sigma must be positive.\n");
    return 1;
  }

  // lattice defaults of Texture_process
  if (!b) {
    b = a;
  }
  if (!c) {
    c = b;
  }
  if (!alpha) {
    alpha = 90;
  }
  if (!beta) {
    beta = alpha;
  }
  if (!gamma) {
    gamma = beta;
  }
  for (i = 0; i < 3; i++) {
    directions[i] = vector_set (direction_values[3 * i], direction_values[3 * i + 1], direction_values[3 * i + 2]);
  }
  crystal_lattice (a, b, c, alpha, beta, gamma, mtex_axes, directions, lattice);

  if (ebsd_file) {
    n_orientations = read_ebsd (ebsd_file, &orientations);
    write_orientation_coefficients (orientations, n_orientations, sigma_deg * (M_PI / 180), lmax, lattice, ebsd_file, out_file);
  } else {
    orientations = draw_expression_orientations (expressions, n_orientations, constants, seed);
    write_orientation_coefficients (orientations, n_orientations, sigma_deg * (M_PI / 180), lmax, lattice, "rotation expressions", out_file);
  }
  free (orientations);
  return 0;
}
