#!/usr/bin/env python3
"""
force_extract.py
Usage: python3 force_extract.py \
         --system output/system.xml \
         --pdb    output/state_topology.pdb \
         --out-tsv dna_forces_per_residue.tsv \
         --out-pdb dna_forces_bfactor.pdb

Produces:
 - tsv of average force per residue (pN)
 - pdb with B-factor set to residue avg force for visualization
"""

from openmm import XmlSerializer, LangevinIntegrator, Platform, unit
from openmm.app import PDBFile, Simulation
import numpy as np
from scipy.constants import N_A
import argparse
import sys

def kj_per_nm_per_mol_to_pN_factor():
    # Build the quantity: 1 kilojoule / nanometer  (per molecule: divide by N_A)
    # then convert to piconewton numerically
    # Use unit.kilojoule (not per mole), then divide by Avogadro constant
    q = (1.0 * unit.kilojoule / unit.nanometer) / N_A
    return q.value_in_unit(unit.piconewton)  # numeric factor (≈1.66054)

def vec_magnitude_kjpm_per_nm_to_pN(vec):
    # vec is an openmm quantity Vec3 or sequence with units kJ/(mol·nm) components OR floats
    # We will try to read components in the force units (kJ/mol/nm), compute magnitude, convert
    # to pN using factor above.
    try:
        # guarded attempt: each component is a Quantity (has value_in_unit)
        fx = vec[0].value_in_unit(unit.kilojoule_per_mole/unit.nanometer)
        fy = vec[1].value_in_unit(unit.kilojoule_per_mole/unit.nanometer)
        fz = vec[2].value_in_unit(unit.kilojoule_per_mole/unit.nanometer)
        mag_kjpm_per_nm = np.sqrt(fx*fx + fy*fy + fz*fz)
    except Exception:
        # fallback: numeric floats (already stripped of units)
        fx, fy, fz = float(vec[0]), float(vec[1]), float(vec[2])
        mag_kjpm_per_nm = np.sqrt(fx*fx + fy*fy + fz*fz)

    # convert the kJ/(mol·nm) magnitude to pN
    factor = kj_per_nm_per_mol_to_pN_factor()
    return mag_kjpm_per_nm * factor

def write_pdb_with_bfactor(inp_pdb, out_pdb_path, per_atom_bfactors):
    # per_atom_bfactors: list/array length == number atoms
    # writes a PDB where B-factor column is replaced with provided values
    with open(inp_pdb, 'r') as inf, open(out_pdb_path, 'w') as outf:
        atom_idx = 0
        for line in inf:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                if atom_idx >= len(per_atom_bfactors):
                    # safety: if mismatch, write original line
                    outf.write(line)
                else:
                    bf = per_atom_bfactors[atom_idx]
                    # PDB b-factor field columns 61-66 (1-based). Build new formatted line.
                    # We'll preserve most fields and replace the B-factor field.
                    new_line = line[:60] + f"{bf:6.2f}" + line[66:]
                    outf.write(new_line)
                atom_idx += 1
            else:
                outf.write(line)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--system", "-s", required=True, help="OpenMM system XML (output/system.xml)")
    parser.add_argument("--pdb", "-p", required=True, help="PDB file with matching topology/positions")
    parser.add_argument("--out-tsv", "-t", default="dna_forces_per_residue.tsv")
    parser.add_argument("--out-pdb", "-o", default="dna_forces_bfactor.pdb")
    parser.add_argument("--min-iters", type=int, default=50)
    parser.add_argument("--platform", default="CPU")
    args = parser.parse_args()

    print("Loading system:", args.system)
    with open(args.system, 'r') as f:
        system = XmlSerializer.deserialize(f.read())

    print("Loading PDB:", args.pdb)
    pdb = PDBFile(args.pdb)

    # sanity check
    n_sys = system.getNumParticles()
    n_pdb = sum(1 for _ in pdb.positions)
    topology_atom_count = sum(1 for _ in pdb.topology.atoms())
    print("System particles:", n_sys, "PDB positions:", n_pdb, "topology atoms:", topology_atom_count)
    if n_sys != n_pdb or n_sys != topology_atom_count:
        print("ERROR: system and PDB atom counts do not match. Aborting.", file=sys.stderr)
        sys.exit(1)

    # integrator & simulation
    integrator = LangevinIntegrator(300*unit.kelvin, 1.0/unit.picosecond, 0.002*unit.picoseconds)
    platform = Platform.getPlatformByName(args.platform)
    sim = Simulation(pdb.topology, system, integrator, platform)
    sim.context.setPositions(pdb.positions)

    print("Minimizing energy ({} iterations)...".format(args.min_iters))
    sim.minimizeEnergy(maxIterations=args.min_iters)

    print("Requesting forces from context...")
    # ask for forces
    state = sim.context.getState(getForces=True)
    forces = state.getForces(asNumpy=True)  # returns an array-like of Vec3 (with units kJ/mol/nm)

    # Build per-atom magnitudes in pN
    n_atoms = len(forces)
    print("Number of forces returned:", n_atoms)
    force_mags_pN = np.empty(n_atoms, dtype=float)
    for i, f in enumerate(forces):
        mag_pN = vec_magnitude_kjpm_per_nm_to_pN(f)
        force_mags_pN[i] = float(mag_pN)

    # Map atoms -> residues
    residues = list(pdb.topology.residues())
    atom_to_res_index = []
    for r_idx, res in enumerate(residues):
        for atom in res.atoms():
            atom_to_res_index.append(r_idx)
    atom_to_res_index = np.array(atom_to_res_index, dtype=int)
    if len(atom_to_res_index) != n_atoms:
        print("WARNING: atom->res mapping length doesn't match atom count. Attempting recovery...")
        # fallback: assume atom ordering matches topology atoms()
        atom_to_res_index = np.empty(n_atoms, dtype=int)
        atom_idx = 0
        for r_idx, res in enumerate(residues):
            for _ in res.atoms():
                if atom_idx < n_atoms:
                    atom_to_res_index[atom_idx] = r_idx
                atom_idx += 1

    # aggregate per-residue
    res_force_list = [[] for _ in residues]
    for atom_idx, res_idx in enumerate(atom_to_res_index):
        res_force_list[res_idx].append(force_mags_pN[atom_idx])

    res_avg_force = []
    for r_idx, vals in enumerate(res_force_list):
        if len(vals) == 0:
            avg = 0.0
        else:
            avg = float(np.mean(vals))
        res_avg_force.append(avg)

    # Write TSV
    with open(args.out_tsv, 'w') as out:
        out.write("res_index\tchain\tresname\tresid\tavg_force_pN\n")
        for i, res in enumerate(residues):
            chain = res.chain.id if res.chain is not None else ""
            out.write(f"{i}\t{chain}\t{res.name}\t{res.id}\t{res_avg_force[i]:.6f}\n")

    print("Wrote per-residue forces to:", args.out_tsv)

    # Also write a PDB with B-factor set to residue average (for visualization).
    # Build per-atom b-factors: each atom gets its residue's avg
    per_atom_bfactors = np.empty(n_atoms, dtype=float)
    # iterate atoms in pdb.topology to assign
    atom_idx = 0
    for res_idx, res in enumerate(residues):
        avg = res_avg_force[res_idx]
        for atom in res.atoms():
            per_atom_bfactors[atom_idx] = avg
            atom_idx += 1

    print("Writing PDB with B-factor = avg force (pN):", args.out_pdb)
    write_pdb_with_bfactor(args.pdb, args.out_pdb, per_atom_bfactors)
    print("Done.")

if __name__ == "__main__":
    main()
