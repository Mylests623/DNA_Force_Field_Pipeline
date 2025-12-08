from openmm.app import PDBFile, Modeller
from openmm.unit import Quantity, nanometer, Vec3
import sys

inp, outp = sys.argv[1:]

def normalize_positions(pos):
    out = []
    for p in pos:
        if hasattr(p, "unit"):
            out.append(Vec3(p[0].value_in_unit(nanometer),
                            p[1].value_in_unit(nanometer),
                            p[2].value_in_unit(nanometer)))
        else:
            out.append(Vec3(float(p[0]), float(p[1]), float(p[2])))
    return Quantity(out, nanometer)

def fix_resnames(res):
    n = res.name.upper()
    if n.startswith("DA"): return "DA"
    if n.startswith("DT"): return "DT"
    if n.startswith("DG"): return "DG"
    if n.startswith("DC"): return "DC"
    return res.name

pdb = PDBFile(inp)
mod = Modeller(pdb.topology, pdb.positions)

# Renaming
for chain in mod.topology.chains():
    for res in chain.residues():
        res.name = fix_resnames(res)

# Delete problematic termini
to_delete = []
for chain in mod.topology.chains():
    res = list(chain.residues())
    first, last = res[0], res[-1]

    for atom in first.atoms():
        if atom.name in ("P", "OP1", "OP2", "O5'"):
            to_delete.append(atom)
    for atom in last.atoms():
        if atom.name == "O3'":
            to_delete.append(atom)

mod.delete(to_delete)

# Fix positions
mod.positions = normalize_positions(mod.positions)

# Add hydrogens
from openmm.app import ForceField
ff = ForceField("amber14-all.xml", "amber14/tip3p.xml")
mod.addHydrogens(ff)

with open(outp, "w") as f:
    PDBFile.writeFile(mod.topology, mod.positions, f)
