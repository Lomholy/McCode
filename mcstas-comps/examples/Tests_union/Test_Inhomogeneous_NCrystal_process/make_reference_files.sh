#!/bin/sh
# Writes the .ncmat template and reference material used by
# Test_Inhomogeneous_NCrystal_process.instr, from Al_sg225.ncmat in the NCrystal data
# library. Needs NCrystal's Python module.
set -e
cd "$(dirname "$0")"
python3 - <<'PY'
import NCrystal as NC
text = NC.createTextData("Al_sg225.ncmat").rawData
assert "  cubic 4.04958\n" in text
lengths = "  lengths 4.04958 4.04958 4.04958\n"
# Cubic for x < 0 and gamma = 92 deg for x > 0; @SPACEGROUP 225 is kept, and removed by the process
open("Al_angle_x.ncmat", "w").write(text.replace("  cubic 4.04958\n", lengths + "  angles 90 90 ${hvs(x,0,92,90)}\n"))
# The gamma = 92 deg reference for NCrystal_process, without @SPACEGROUP, which it would violate
distorted = text.replace("  cubic 4.04958\n", lengths + "  angles 90 90 92\n")
start = distorted.index("@SPACEGROUP")
open("Al_gamma92.ncmat", "w").write(distorted[:start] + distorted[distorted.index("@", start + 1):])
PY
