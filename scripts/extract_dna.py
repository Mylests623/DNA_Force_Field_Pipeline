from openmm.app import PDBFile
from openmm.app import Modeller

input_file = "output/validated.pdb"
output_file = "output/dna_only.pdb"

pdb = PDBFile(input_file)

# DNA residue names
dna_names = {"DA", "DT", "DG", "DC", "A", "T", "G", "C"}

# Build list of residues to keep
dna_residues = []

for chain in pdb.topology.chains():
    for res in chain.residues():
        if res.name.upper() in dna_names:
            dna_residues.append(res)

if len(dna_residues) == 0:
    raise ValueError("ERROR: No DNA residues found in validated PDB.")

# Create modeller that contains **ONLY DNA**
modeller = Modeller(pdb.topology, pdb.positions)
modeller.delete(
    atom for atom in pdb.topology.atoms()
    if atom.residue not in dna_residues
)

with open(output_file, "w") as f:
    PDBFile.writeFile(modeller.topology, modeller.positions, f)

print("DNA extraction complete. Output:", output_file)
