#!/bin/sh
# Writes the .ncmat templates and the Single_crystal_process reference files (.lau) used
# by Test_Inhomogeneous_single_crystal_process.instr, from Al_sg225.ncmat in the NCrystal
# data library. Needs NCrystal's command line tools and Python module.
set -e
cd "$(dirname "$0")"
rm -f Al_300K.lau Al_600K.lau Al_a101_300K.lau Al_gamma92_300K.lau

python3 - <<'PY'
import NCrystal as NC
text = NC.createTextData("Al_sg225.ncmat").rawData
assert "  cubic 4.04958\n" in text
lengths = "  lengths 4.04958 4.04958 4.04958\n"
open("Al_lattice_x.ncmat", "w").write(text.replace("  cubic 4.04958\n", "  cubic ${hvs(x,0,p2,p1)}\n"))
open("Al_angle_x.ncmat", "w").write(text.replace("  cubic 4.04958\n", lengths + "  angles 90 90 ${hvs(x,0,92,90)}\n"))
open("Al_a101.ncmat", "w").write(text.replace("  cubic 4.04958\n", f"  cubic {4.04958*1.01:.12g}\n"))
distorted = text.replace("  cubic 4.04958\n", lengths + "  angles 90 90 92\n")
start = distorted.index("@SPACEGROUP")
open("Al_gamma92.ncmat", "w").write(distorted[:start] + distorted[distorted.index("@", start + 1):])
PY

ncrystal_ncmat2hkl "Al_sg225.ncmat;temp=300K" --format=lau -o Al_300K.lau --quiet
ncrystal_ncmat2hkl "Al_sg225.ncmat;temp=600K" --format=lau -o Al_600K.lau --quiet
ncrystal_ncmat2hkl "Al_a101.ncmat;temp=300K" --format=lau -o Al_a101_300K.lau --quiet
# Single_crystal_process builds its own lattice (a along y, b along z) from lattice_* lines in
# the header, ignoring the lattice vectors given to it; without them it uses the given vectors,
# which have a along x and b in the xy plane, as Inhomogeneous_single_crystal_process
for f in Al_300K.lau Al_600K.lau Al_a101_300K.lau; do
  grep -v "^# lattice_" "$f" > "$f.tmp" && mv "$f.tmp" "$f"
done

# ncrystal_ncmat2hkl needs a space group number, so this .lau is written directly from NCrystal:
# every (h,k,l) of each family and its negative, with the family's |F|^2
python3 - <<'PY'
import NCrystal as NC
import numpy as np
cfg = "Al_gamma92.ncmat;temp=300K"
info = NC.createInfo(cfg)
with open("Al_gamma92_300K.lau", "w") as f:
    f.write(f"# Reflections of \"{cfg}\" from NCrystal {NC.__version__}, in Single_crystal .lau format\n")
    f.write("# column_h  1\n# column_k  2\n# column_l  3\n# column_j  4\n# column_d  5  [d-spacing in Aa]\n# column_F2 6  [norm of scattering factor |F|^2 in barn]\n")
    for fam in info.hklObjects():
        for h, k, l in zip(np.atleast_1d(fam.h), np.atleast_1d(fam.k), np.atleast_1d(fam.l)):
            for sign in (1, -1):
                f.write(f"{sign*int(h):4d} {sign*int(k):4d} {sign*int(l):4d} 1 {float(fam.d)!r} {float(fam.f2)!r}\n")
PY
rm Al_a101.ncmat Al_gamma92.ncmat
