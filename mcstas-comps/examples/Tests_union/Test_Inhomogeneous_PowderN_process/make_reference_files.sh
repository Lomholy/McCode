#!/bin/sh
# Writes the ncmat templates and the Powder_process reference files (.laz) used by
# Test_Inhomogeneous_PowderN_process.instr. Needs NCrystal's command line tools.
set -e
cd "$(dirname "$0")"
rm -f Al_300K.laz Al_a101_300K.laz Al_600K.laz Al_gamma92_300K.laz

A=4.04958
A101=$(python3 -c "print(f'{4.04958*1.01:.12g}')")

# Al_sg225.ncmat from the NCrystal data library, with the lattice parameter as a template
python3 - <<'PY'
import NCrystal as NC
text = NC.createTextData("Al_sg225.ncmat").rawData
assert "  cubic 4.04958\n" in text
for name, cell in (("Al_lattice_x.ncmat", "${hvs(x,0,p2,p1)}"), ("Al_lattice_t.ncmat", "${hvs(t,p3,p2,p1)}"), ("Al_a101.ncmat", f"{4.04958*1.01:.12g}")):
    open(name, "w").write(text.replace("  cubic 4.04958\n", f"  cubic {cell}\n"))

# Cubic for x < 0 and gamma = 92 deg (lower symmetry) for x > 0. The template keeps
# @SPACEGROUP 225, which Inhomogeneous_PowderN_process removes. The reference for the
# distorted cell has no @SPACEGROUP, as NCrystal rejects it with space group 225.
lengths = "  lengths 4.04958 4.04958 4.04958\n"
open("Al_angle_x.ncmat", "w").write(text.replace("  cubic 4.04958\n", lengths + "  angles 90 90 ${hvs(x,0,92,90)}\n"))
distorted = text.replace("  cubic 4.04958\n", lengths + "  angles 90 90 92\n")
start = distorted.index("@SPACEGROUP")
open("Al_gamma92.ncmat", "w").write(distorted[:start] + distorted[distorted.index("@", start + 1):])
PY

ncrystal_ncmat2hkl "Al_sg225.ncmat;temp=300K" --format=laz -o Al_300K.laz --quiet
ncrystal_ncmat2hkl "Al_a101.ncmat;temp=300K" --format=laz -o Al_a101_300K.laz --quiet
ncrystal_ncmat2hkl "Al_sg225.ncmat;temp=600K" --format=laz -o Al_600K.laz --quiet
# ncrystal_ncmat2hkl needs a space group number, so this .laz is written directly from NCrystal
python3 - <<'PY'
import NCrystal as NC
import numpy as np
cfg = "Al_gamma92.ncmat;temp=300K"
info = NC.createInfo(cfg)
with open("Al_gamma92_300K.laz", "w") as f:
    f.write(f"# Reflections of \"{cfg}\" from NCrystal {NC.__version__}, in PowderN .laz format\n")
    f.write(f"# Vc {info.structure_info['volume']!r} [Aa^3]\n")
    f.write("# column_h  1\n# column_k  2\n# column_l  3\n")
    f.write("# column_j  4  multiplicity\n# column_d  5  [d-spacing in Aa]\n# column_F2 6  [norm of scattering factor |F|^2 in barn]\n")
    for hkl in info.hklObjects():
        h, k, l = (int(np.atleast_1d(v)[0]) for v in (hkl.h, hkl.k, hkl.l))
        f.write(f"{h:4d} {k:4d} {l:4d} {hkl.mult:3d} {float(hkl.d)!r} {float(hkl.f2)!r}\n")
PY
rm Al_a101.ncmat Al_gamma92.ncmat
