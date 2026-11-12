#!/usr/bin/env python3
"""
McStas, neutron ray-tracing package

Program: share/texture-coefficients.py
Written by: Daniel Lomholt Christensen
Date: 12/11/2026
Origin: Niels Bohr Institute @ UCPH

Writes the Fourier coefficients of an orientation distribution function (ODF) to the
file that Texture_process reads as fcoef_fn, with the five columns l m n Real(Clmn) Imag(Clmn).
When --out ends in .bin, the file is binary instead, about 3.5 times smaller: only
Real(Clmn) and Imag(Clmn) of each row, as two little-endian 8-byte doubles, without header,
in the row order of the text file, from which Texture_process recovers l, m, n and lmax.
The ODF is given by a list of orientations, either an orientation list such as an EBSD
map, or orientations drawn from three expressions. Needs Python 3 and numpy.

The orientations rotate the Cartesian crystal axes into the axes of the component.
The crystal axes are those of Texture_process: by default those of NCrystal, x || a and
b in the xy plane; --mtex_axes 1 takes those of MTEX, z || c and x || a*; --ax ... --cz place
the lattice vectors along given directions in the component. Give the same lattice,
mtex_axes and ax ... cz as to Texture_process.

The ODF is the orientations smoothed by a Gaussian kernel on the rotation group, with the
Fourier coefficients kappa_l = exp(-l(l+1) sigma^2 / 2), averaged over the proper rotations
of the lattice. Each orientation has the same weight, and the computation time grows as
the number of orientations times lmax^4.

--ebsd: each row starts with the Bunge Euler angles phi1 Phi phi2 in radians of one
orientation, crystal to sample; further columns are ignored, and rows that do not start
with three numbers are skipped.

--rot_x_expr, --rot_y_expr, --rot_z_expr: the angles in deg of rotations of the crystal axes
about the x, y and z axes of the component, composed as ROTATED (rx, ry, rz) composes
them: g = Rx (rx) Ry (ry) Rz (rz). The expressions are evaluated by Python with numpy, in
the random numbers u1, u2, u3 in [0, 1), drawn anew for each of the n_orientations
orientations from the generator seeded by --seed, the constants p1 ... p6, pi, and the
functions sqrt, exp, log, sin, cos, tan, asin, acos, atan, atan2, abs; ^ is a power.
The noise of the drawn texture falls as 1/sqrt(n_orientations) and grows as
orientation_sigma decreases, about as orientation_sigma^-3 in the number of orientations needed.

Example, a fibre texture with a tilted from x by a Gaussian of 5 deg and turned about x:
  python3 texture-coefficients.py --out fibre.txt --a 4.04958 --mtex_axes 1 --p1 5 \\
    --rot_x_expr "360*u3" --rot_y_expr "p1*sqrt(-2*log(1-u1))*cos(2*pi*u2)" --rot_z_expr "0"
"""

import argparse
import itertools
import sys

import numpy as np


def fail(message):
    sys.exit("texture-coefficients: " + message)


def orthonormal_basis(v1, v2):
    """Rows of a right-handed orthonormal basis: e0 || v1, e1 in the plane of v1 and v2."""
    normal = np.cross(v1, v2)
    if np.linalg.norm(v1) == 0 or np.linalg.norm(normal) < 1e-12 * np.linalg.norm(v1) * np.linalg.norm(v2):
        fail("the directions of the two lattice vectors must be nonzero and not parallel.")
    e0 = v1 / np.linalg.norm(v1)
    e2 = normal / np.linalg.norm(normal)
    return np.array([e0, np.cross(e2, e0), e2])


def rotation_between_directions(c1, c2, l1, l2):
    """Rotation that takes crystal direction c1 onto component direction l1, and c2 into the plane of l1 and l2."""
    return orthonormal_basis(np.asarray(l1, float), np.asarray(l2, float)).T @ orthonormal_basis(c1, c2)


def lattice_along_directions(lattice, directions):
    """Lattice vectors along the three directions, with the lengths of lattice."""
    unit = directions / np.linalg.norm(directions, axis=1)[:, None]
    for i, name in enumerate("abc"):
        i_next, i_last = (i + 1) % 3, (i + 2) % 3
        angle_lattice = np.arccos(lattice[i_next] @ lattice[i_last] / (np.linalg.norm(lattice[i_next]) * np.linalg.norm(lattice[i_last])))
        angle_directions = np.arccos(np.clip(unit[i_next] @ unit[i_last], -1, 1))
        if abs(angle_directions - angle_lattice) > 1e-4:
            fail(f"the directions of {'abc'[i_next]} and {'abc'[i_last]} are {np.degrees(angle_directions):g} deg apart, "
                 f"but {np.degrees(angle_lattice):g} deg in the lattice.")
    if unit[0] @ np.cross(unit[1], unit[2]) <= 0:
        fail("the directions of a, b and c must be right-handed.")
    return unit * np.linalg.norm(lattice, axis=1)[:, None]


def crystal_lattice(a, b, c, alpha, beta, gamma, mtex_axes, directions):
    """Rows a, b, c of the lattice vectors in the Cartesian crystal frame of Texture_process."""
    alpha, beta, gamma = np.radians([alpha, beta, gamma])
    # NCrystal: x || a and b in the xy plane, so z || c*
    cx = c * np.cos(beta)
    cy = c * (np.cos(alpha) - np.cos(beta) * np.cos(gamma)) / np.sin(gamma)
    lattice = np.array([[a, 0, 0], [b * np.cos(gamma), b * np.sin(gamma), 0], [cx, cy, np.sqrt(c * c - cx * cx - cy * cy)]])

    given = [i for i in range(3) if np.linalg.norm(directions[i]) > 0]
    if len(given) == 1:
        fail("give the directions of two or three of the lattice vectors a, b and c, or none.")
    if given and mtex_axes:
        fail("the directions of the lattice vectors cannot be combined with mtex_axes = 1.")
    if len(given) == 3:
        return lattice_along_directions(lattice, directions)
    if mtex_axes:
        rotation = rotation_between_directions(lattice[2], np.cross(lattice[1], lattice[2]), [0, 0, 1], [1, 0, 0])
    elif len(given) == 2:
        rotation = rotation_between_directions(lattice[given[0]], lattice[given[1]], directions[given[0]], directions[given[1]])
    else:
        return lattice
    return lattice @ rotation.T


def lattice_rotations(lattice):
    """Proper rotations of the Cartesian crystal frame that map the lattice onto itself."""
    A = lattice.T
    # integer matrices M with entries -1, 0, 1 in the lattice basis, R = A M A^-1
    M = np.array(list(itertools.product((-1, 0, 1), repeat=9)), float).reshape(-1, 3, 3)
    R = A @ M @ np.linalg.inv(A)
    deviation = np.abs(np.transpose(R, (0, 2, 1)) @ R - np.eye(3)).sum(axis=(1, 2))
    return R[(deviation < 1e-6) & (np.linalg.det(R) > 0)]


def harmonic_sums(rotations, lmax, points):
    """Sum over the rotations g of N_lm (g h) = scaled P_l^m (cos theta) exp (i m phi) at the points h,
    for 0 <= m <= l <= lmax, at row l * (lmax + 1) + m, with N_lm = sqrt ((l - m)! / (l + m)!) P_l^m."""
    total = np.zeros(((lmax + 1) ** 2, len(points)), complex)
    chunk = max(1, 200000 // len(points))
    for start in range(0, len(rotations), chunk):
        directions = np.einsum("gij,pj->gpi", rotations[start:start + chunk], points)
        x = np.clip(directions[..., 2], -1, 1)
        phi = np.arctan2(directions[..., 1], directions[..., 0])
        sin_theta = np.sqrt((1 - x) * (1 + x))
        p_mm = np.ones_like(x)
        for m in range(lmax + 1):
            if m > 0:
                p_mm = -(2 * m - 1) * sin_theta * p_mm
            phase = np.exp(1j * m * phi)
            # P_l^m by recurrence in the degree, from P_m^m and P_(m+1)^m
            p_previous, p_l = np.zeros_like(x), p_mm
            for l in range(m, lmax + 1):
                if l > m:
                    p_previous, p_l = p_l, (x * (2 * l - 1) * p_l - (l + m - 1) * p_previous) / (l - m)
                scale = np.sqrt(1 / np.prod(np.arange(l - m + 1, l + m + 1, dtype=float)))
                total[l * (lmax + 1) + m] += scale * (p_l * phase).sum(axis=0)
    return total


def harmonic_rows(table, l, lmax):
    """Rows m = -l ... l of degree l of a table of harmonic_sums, with N_l,-m = (-1)^m conj (N_lm)."""
    positive = table[l * (lmax + 1):l * (lmax + 1) + l + 1]
    negative = [(-1) ** m * np.conj(positive[m]) for m in range(l, 0, -1)]
    return np.vstack(negative + [positive])


def orientation_coefficients(orientations, sigma, lmax, lattice):
    """C_lmn = (2l + 1) kappa_l [avg_i W (g_i) avg_s W (s)]_mn at [l][l + m, l + n], where N_lm (g h) = sum_n W (g)_mn N_ln (h),
    with W from the quadrature of the harmonics on the sphere, exact for the products of two harmonics of degree <= lmax."""
    n_phi = 2 * lmax + 4
    nodes, node_weights = np.polynomial.legendre.leggauss(lmax + 2)
    azimuth = 2 * np.pi * np.arange(n_phi) / n_phi
    sin_theta = np.sqrt(1 - nodes ** 2)
    points = np.stack([np.outer(sin_theta, np.cos(azimuth)), np.outer(sin_theta, np.sin(azimuth)), np.repeat(nodes[:, None], n_phi, axis=1)], -1).reshape(-1, 3)
    weights = np.repeat(node_weights * 2 * np.pi / n_phi, n_phi)

    symmetry = lattice_rotations(lattice)
    orientation_table = harmonic_sums(orientations, lmax, points) / len(orientations)
    symmetry_table = harmonic_sums(symmetry, lmax, points) / len(symmetry)
    reference_table = harmonic_sums(np.eye(3)[None], lmax, points)

    coefficients = []
    for l in range(lmax + 1):
        size = 2 * l + 1
        reference = np.conj(harmonic_rows(reference_table, l, lmax)).T
        orientation_matrix = size / (4 * np.pi) * (harmonic_rows(orientation_table, l, lmax) * weights) @ reference
        symmetry_matrix = size / (4 * np.pi) * (harmonic_rows(symmetry_table, l, lmax) * weights) @ reference
        coefficients.append(size * np.exp(-l * (l + 1) * sigma ** 2 / 2) * orientation_matrix @ symmetry_matrix)
    return coefficients, len(symmetry)


def write_coefficients(filename, coefficients, comment):
    """Rows l, n, m, Re, Im with Re + i Im = s conj (C_lmn), and s = -1 when the negative ones of m and n sum to an odd number:
    the inverse of readFourierCoef of Texture_process."""
    rows = []
    for l, C in enumerate(coefficients):
        for n, m in itertools.product(range(-l, l + 1), repeat=2):
            sign = -1 if (max(-m, 0) + max(-n, 0)) % 2 else 1
            value = C[l + m, l + n]
            rows.append((l, n, m, sign * value.real, -sign * value.imag))
    if filename.endswith(".bin"):
        np.asarray(rows)[:, 3:].astype("<f8").tofile(filename)
        return
    header = comment + "\nwritten by texture-coefficients.py for fcoef_fn of Texture_process\nl m n Real(Clmn) Imag(Clmn)"
    np.savetxt(filename, rows, fmt=["%d", "%d", "%d", "%.17g", "%.17g"], header=header)


def read_ebsd(filename):
    """Orientations of the rows of filename, Bunge Euler angles g = Rz (phi1) Rx (Phi) Rz (phi2)."""
    angles = []
    with open(filename) as file:
        for line in file:
            fields = line.replace(",", " ").replace(";", " ").split()
            try:
                angles.append([float(field) for field in fields[:3]])
            except ValueError:
                continue
            if len(angles[-1]) < 3:
                angles.pop()
    if not angles:
        fail(f"{filename} must have rows starting with the Euler angles phi1 Phi phi2.")
    phi1, Phi, phi2 = np.array(angles).T
    c1, s1, C, S, c2, s2 = np.cos(phi1), np.sin(phi1), np.cos(Phi), np.sin(Phi), np.cos(phi2), np.sin(phi2)
    return np.stack([np.stack([c1 * c2 - s1 * s2 * C, -c1 * s2 - s1 * c2 * C, s1 * S], -1),
                     np.stack([s1 * c2 + c1 * s2 * C, -s1 * s2 + c1 * c2 * C, -c1 * S], -1),
                     np.stack([s2 * S, c2 * S, C], -1)], -2)


def draw_expression_orientations(expressions, n_orientations, constants, seed):
    """Orientations from the expressions of the angles (deg) about x, y and z, composed as ROTATED:
    g = Rx (rx) Ry (ry) Rz (rz), active, crystal axes to component axes."""
    if n_orientations < 1:
        fail("n_orientations must be positive.")
    u = np.random.default_rng(seed).random((3, n_orientations))
    names = {"u1": u[0], "u2": u[1], "u3": u[2], "pi": np.pi, "sqrt": np.sqrt, "exp": np.exp, "log": np.log, "sin": np.sin,
             "cos": np.cos, "tan": np.tan, "asin": np.arcsin, "acos": np.arccos, "atan": np.arctan, "atan2": np.arctan2, "abs": np.abs}
    names.update({f"p{i + 1}": value for i, value in enumerate(constants)})
    angles = []
    for expression in expressions:
        try:
            angle = eval(expression.replace("^", "**"), {"__builtins__": {}}, names)
        except Exception as error:
            fail(f'cannot evaluate the rotation expression "{expression}": {error}')
        angle = np.radians(np.broadcast_to(angle, (n_orientations,)).astype(float))
        if not np.isfinite(angle).all():
            fail(f'the rotation expression "{expression}" is not finite for some u1, u2, u3.')
        angles.append(angle)
    cx, sx, cy, sy, cz, sz = np.cos(angles[0]), np.sin(angles[0]), np.cos(angles[1]), np.sin(angles[1]), np.cos(angles[2]), np.sin(angles[2])
    # the transpose of the coordinate transformation of ROTATED (rot_set_rotation)
    return np.stack([np.stack([cy * cz, -cy * sz, sy], -1),
                     np.stack([sx * sy * cz + cx * sz, cx * cz - sx * sy * sz, -sx * cy], -1),
                     np.stack([sx * sz - cx * sy * cz, sx * cz + cx * sy * sz, cx * cy], -1)], -2)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--out", required=True, help="file to write; binary when it ends in .bin")
    parser.add_argument("--a", type=float, required=True, help="lattice parameter a [AA]")
    parser.add_argument("--b", type=float, default=0, help="lattice parameter b [AA], default a")
    parser.add_argument("--c", type=float, default=0, help="lattice parameter c [AA], default b")
    parser.add_argument("--alpha", type=float, default=0, help="lattice angle [deg], default 90")
    parser.add_argument("--beta", type=float, default=0, help="lattice angle [deg], default alpha")
    parser.add_argument("--gamma", type=float, default=0, help="lattice angle [deg], default beta")
    parser.add_argument("--mtex_axes", type=int, default=0, help="0: crystal axes of NCrystal, 1: of MTEX, as in Texture_process")
    for axis in ["ax", "ay", "az", "bx", "by", "bz", "cx", "cy", "cz"]:
        parser.add_argument("--" + axis, type=float, default=0, help="direction of a lattice vector in the component, as in Texture_process")
    parser.add_argument("--orientation_sigma", type=float, default=10, help="standard deviation of the kernel [deg]")
    parser.add_argument("--lmax", type=int, default=-1, help="cut-off of the series; negative: the l where the kernel has fallen to 1e-3")
    parser.add_argument("--ebsd", help="orientation list")
    parser.add_argument("--rot_x_expr", help="angle about x [deg], applied last")
    parser.add_argument("--rot_y_expr", help="angle about y [deg]")
    parser.add_argument("--rot_z_expr", help="angle about z [deg], applied first")
    parser.add_argument("--n_orientations", type=int, default=10000, help="number of orientations drawn from the expressions")
    for i in range(1, 7):
        parser.add_argument(f"--p{i}", type=float, default=0, help="constant of the expressions")
    parser.add_argument("--seed", type=int, default=1, help="seed of the random numbers of the expressions")
    args = parser.parse_args()

    expressions = [args.rot_x_expr, args.rot_y_expr, args.rot_z_expr]
    n_expressions = sum(expression is not None for expression in expressions)
    if n_expressions not in (0, 3):
        fail("give all three of rot_x_expr, rot_y_expr and rot_z_expr, or none.")
    if (args.ebsd is not None) + (n_expressions == 3) != 1:
        fail("give the orientations in exactly one of ebsd and the rotation expressions.")
    if args.a <= 0 or args.orientation_sigma <= 0:
        fail("a and orientation_sigma must be positive.")

    # lattice defaults of Texture_process
    b = args.b or args.a
    c = args.c or b
    alpha = args.alpha or 90
    beta = args.beta or alpha
    gamma = args.gamma or beta
    directions = np.array([[args.ax, args.ay, args.az], [args.bx, args.by, args.bz], [args.cx, args.cy, args.cz]])
    lattice = crystal_lattice(args.a, b, c, alpha, beta, gamma, args.mtex_axes, directions)

    if args.ebsd:
        orientations, source = read_ebsd(args.ebsd), args.ebsd
    else:
        constants = [args.p1, args.p2, args.p3, args.p4, args.p5, args.p6]
        orientations, source = draw_expression_orientations(expressions, args.n_orientations, constants, args.seed), "rotation expressions"

    sigma = np.radians(args.orientation_sigma)
    lmax = args.lmax if args.lmax >= 0 else int(np.ceil(np.sqrt(2 * np.log(1000.0)) / sigma))
    coefficients, n_symmetry = orientation_coefficients(orientations, sigma, lmax, lattice)
    summary = (f"{source}: {len(orientations)} orientations, {n_symmetry} lattice rotations, "
               f"sigma = {args.orientation_sigma:g} deg, lmax = {lmax}")
    write_coefficients(args.out, coefficients, "Fourier coefficients of the ODF of " + summary)
    print(f"texture-coefficients: {summary}, written to {args.out}")


if __name__ == "__main__":
    main()
