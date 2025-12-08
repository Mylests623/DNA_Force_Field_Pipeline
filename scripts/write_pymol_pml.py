#!/usr/bin/env python3
import sys
from pathlib import Path
import numpy as np

if len(sys.argv) < 3:
    print("Usage: write_pymol_pml.py <in_pdb> <out_pml> <out_png> [width height dpi]")
    sys.exit(1)
in_pdb = Path(sys.argv[1])
out_pml = Path(sys.argv[2])
out_png = Path(sys.argv[3])
width = int(sys.argv[4]) if len(sys.argv) > 4 else 1600
height = int(sys.argv[5]) if len(sys.argv) > 5 else 1200
dpi = int(sys.argv[6]) if len(sys.argv) > 6 else 150

# read b-factors to set a sensible min/max
# crude parse: read column 61-66 (tempFactor) from ATOM lines
bfactors = []
with open(in_pdb, "r") as f:
    for line in f:
        if line.startswith(("ATOM","HETATM")) and len(line) >= 66:
            try:
                bf = float(line[60:66].strip())
                bfactors.append(bf)
            except:
                pass
if len(bfactors) == 0:
    vmin = 0.0
    vmax = 1.0
else:
    vmin = float(np.percentile(bfactors, 2))
    vmax = float(np.percentile(bfactors, 98))
    if vmax <= vmin:
        vmax = max(vmin+1.0, vmax)

pml = f"""
# PyMOL heatmap script - generated
reinitialize
load {in_pdb}
hide everything, all
# show cartoon for nucleic acid backbone and sticks for bases (adjust as desired)
show cartoon, polymer and not solvent
show sticks, (polymer and not solvent)
# color by B-factor using blue-white-red (low -> high)
spectrum b, blue_white_red, minimum={vmin:.4f}, maximum={vmax:.4f}
# optionally set cartoon thickness and stick scale
set cartoon_sampling, 10
set cartoon_transparency, 0.35
bg_color white
set ray_trace_mode, 1
# orient nicely
orient
# ray and save image
png {out_png}, width={width}, height={height}, dpi={dpi}, ray=1
quit
"""
out_pml.write_text(pml)

# Call PyMOL headless to render (if available)
import subprocess
try:
    subprocess.run(["pymol", "-cq", str(out_pml)], check=True)
    print("Rendered PNG:", out_png)
except FileNotFoundError:
    print("PyMOL not found — created PML script only:", out_pml)
except subprocess.CalledProcessError:
    print("PyMOL call failed; created PML script only:", out_pml)
