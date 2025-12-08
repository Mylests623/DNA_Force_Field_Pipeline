#!/usr/bin/env python3
import sys
import numpy as np
from Bio.PDB import PDBParser, PDBIO

# --- snakemake inputs ---
try:
    snakemake  # type: ignore
    full_pdb = snakemake.input.full
    dna_sim_pdb = snakemake.input.dna
    force_tsv = snakemake.input.tsv
    out_pdb = snakemake.output.pdb
except NameError:
    if len(sys.argv) < 4:
        print("Usage: map_forces_to_full_complex.py <full_complex.pdb> <dna_sim_topology.pdb> <force.tsv> [out.pdb]")
        sys.exit(1)
    full_pdb = sys.argv[1]
    dna_sim_pdb = sys.argv[2]
    force_tsv = sys.argv[3]
    out_pdb = sys.argv[4] if len(sys.argv) > 4 else "full_complex_with_dna_forces.pdb"

# ------------------------------------------------------------
# 1) Load force table (robust TSV parser)
# ------------------------------------------------------------
force_data = {}
with open(force_tsv) as f:
    header = next(f)
    for line in f:
        parts = line.strip().split("\t")
        if len(parts) < 4:
            continue
        chain = parts[0]
        resid = parts[1]          # keep EXACT STRING, do NOT convert to int
        resname = parts[2]
        force = float(parts[3])
        force_key = (chain, resid, resname)
        force_data[force_key] = force

# ------------------------------------------------------------
# 2) Load full complex + DNA sim structures using permissive parser
# ------------------------------------------------------------
parser = PDBParser(PERMISSIVE=True, QUIET=True)
full_struct = parser.get_structure("full", full_pdb)
dna_struct = parser.get_structure("dna", dna_sim_pdb)

# ------------------------------------------------------------
# 3) Build mapping from DNA-only PDB atom -> residue force
# ------------------------------------------------------------
dna_force_per_atom = {}

for model in dna_struct:
    for chain in model:
        for residue in chain:
            chain_id = chain.id
            resid_str = residue.id[1]   # this is a string because PERMISSIVE=True
            resname = residue.resname

            key = (chain_id, str(resid_str), resname)
            if key not in force_data:
                continue

            force_val = force_data[key]

            for atom in residue:
                dna_force_per_atom[(chain_id, str(resid_str), atom.name)] = force_val

# ------------------------------------------------------------
# 4) Write B-factors into the full complex for matching atoms
# ------------------------------------------------------------
for model in full_struct:
    for chain in model:
        for residue in chain:
            chain_id = chain.id
            resid_str = residue.id[1]          # again, string not integer
            resname = residue.resname

            for atom in residue:
                key = (chain_id, str(resid_str), atom.name)
                if key in dna_force_per_atom:
                    atom.bfactor = dna_force_per_atom[key]
                else:
                    atom.bfactor = 0.0

# ------------------------------------------------------------
# 5) Output final annotated PDB
# ------------------------------------------------------------
io = PDBIO()
io.set_structure(full_struct)
io.save(out_pdb)

print("Wrote:", out_pdb)
