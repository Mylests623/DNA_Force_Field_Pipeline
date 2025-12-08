#!/usr/bin/env python3
import sys
import numpy as np
from openmm.app import PDBFile
from openmm.unit import nanometer
from collections import OrderedDict
from pathlib import Path

# ------------------------------------------------------------
# LOAD INPUTS (via Snakemake or CLI)
# ------------------------------------------------------------
try:
    snakemake  # type: ignore
    forces_path = snakemake.input.forces
    topo_path = snakemake.input.topo
    out_tsv = snakemake.output.tsv
    out_pdb = snakemake.output.bpfb_pdb
except NameError:
    if len(sys.argv) < 3:
        print("Usage: aggregate_forces.py <forces.npy> <state_topology.pdb> [out.tsv] [out_bfactor.pdb]")
        sys.exit(1)
    forces_path = sys.argv[1]
    topo_path = sys.argv[2]
    out_tsv = sys.argv[3] if len(sys.argv) > 3 else "dna_forces_per_residue.tsv"
    out_pdb = sys.argv[4] if len(sys.argv) > 4 else "dna_forces_bfactor.pdb"

# ------------------------------------------------------------
# LOAD DATA
# ------------------------------------------------------------
forces = np.load(forces_path)  # (N,3) numpy array, unitless
pdb = PDBFile(topo_path)

atom_list = list(pdb.topology.atoms())
n_atoms = len(atom_list)

if forces.shape[0] != n_atoms:
    raise SystemExit(f"ERROR: forces rows {forces.shape[0]} != atoms in PDB {n_atoms}")

# ------------------------------------------------------------
# FORCE MAGNITUDES → pN
# ------------------------------------------------------------
# 1 kJ/mol/nm = 1.66053906660 pN
KJPM_PER_NM_TO_pN = 1.66053906660

magnitudes_kjpm_per_nm = np.linalg.norm(forces, axis=1)
magnitudes_pN = magnitudes_kjpm_per_nm * KJPM_PER_NM_TO_pN

# ------------------------------------------------------------
# GROUP BY RESIDUE
# ------------------------------------------------------------
res_map = OrderedDict()

for idx, atom in enumerate(atom_list):
    res = atom.residue
    res_id = (res.chain.id, res.id, res.name)
    res_map.setdefault(res_id, []).append(idx)

# Average per residue
res_avg_pN = OrderedDict()
for res_id, idx_list in res_map.items():
    res_avg_pN[res_id] = float(np.mean(magnitudes_pN[idx_list]))

# ------------------------------------------------------------
# WRITE TSV
# ------------------------------------------------------------
with open(out_tsv, "w") as f:
    f.write("chain\tresid\tresname\tavg_force_pN\n")
    for (chainid, resid, resname), avg in res_avg_pN.items():
        f.write(f"{chainid}\t{resid}\t{resname}\t{avg:.6f}\n")

print("Wrote per-residue TSV to:", out_tsv)

# ------------------------------------------------------------
# WRITE PDB WITH B-FACTORS = avg force per residue
# ------------------------------------------------------------
out_lines = []

pos_nm = pdb.positions.value_in_unit(nanometer)
atom_index = 0

for atom, pos in zip(atom_list, pos_nm):
    res = atom.residue
    key = (res.chain.id, res.id, res.name)

    bval = res_avg_pN.get(key, 0.0)

    serial = atom_index + 1
    name = atom.name
    resname = res.name
    chainid = res.chain.id

    # Convert residue ID safely
    try:
        resseq = int(res.id)
    except:
        resseq = 0

    x, y, z = float(pos[0]), float(pos[1]), float(pos[2])
    occupancy = 1.00
    tempFactor = bval
    element = (atom.element.symbol if atom.element else "X")

    line = (
        "ATOM  {serial:5d} {name:^4s}{altloc:1s}{resname:>3s} {chain:1s}"
        "{resseq:4d}{icode:1s}   "
        "{x:8.3f}{y:8.3f}{z:8.3f}"
        "{occ:6.2f}{temp:6.2f}          "
        "{elem:>2s}\n"
    ).format(
        serial=serial, name=name, altloc='', resname=resname, chain=chainid,
        resseq=resseq, icode='', x=x, y=y, z=z,
        occ=occupancy, temp=tempFactor, elem=element
    )

    out_lines.append(line)
    atom_index += 1

# ------------------------------------------------------------
# INSERT TER AND END
# ------------------------------------------------------------
final_lines = []
for i, atom in enumerate(atom_list):
    final_lines.append(out_lines[i])

    if i == len(atom_list) - 1:
        final_lines.append("TER\n")
        break

    curr = atom_list[i]
    nxt = atom_list[i+1]

    if nxt.residue.chain.id != curr.residue.chain.id or nxt.residue.id != curr.residue.id:
        final_lines.append("TER\n")

final_lines.append("END\n")

with open(out_pdb, "w") as f:
    f.writelines(final_lines)

print("Wrote PDB with per-residue force in B-factor:", out_pdb)
