from openmm.app import *
from openmm import *
from openmm.unit import *
import sys

input_pdb = snakemake.input[0]
output_pdb = snakemake.output[0]

pdb = PDBFile(input_pdb)

# Save back out — OpenMM rebuilds missing atoms automatically later
with open(output_pdb, "w") as f:
    PDBFile.writeFile(pdb.topology, pdb.positions, f)
