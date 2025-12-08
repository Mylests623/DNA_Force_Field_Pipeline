from openmm.app import *
from openmm import *
from openmm.unit import *
import sys

inp = snakemake.input[0]
xml_out = snakemake.output["xml"]
pdb_out = snakemake.output["pdb"]

pdb = PDBFile(inp)
forcefield = ForceField("amber14-all.xml", "amber14/tip3p.xml")

# Add hydrogens
modeller = Modeller(pdb.topology, pdb.positions)
modeller.addHydrogens(forcefield, pH=7.0)

# Solvate
modeller.addSolvent(forcefield, model="tip3p", padding=1*nanometer)

system = forcefield.createSystem(
    modeller.topology,
    nonbondedMethod=PME,
    nonbondedCutoff=1*nanometer,
    constraints=HBonds
)

# Save solvated structure
with open(pdb_out, "w") as f:
    PDBFile.writeFile(modeller.topology, modeller.positions, f)

# Save serialized system
with open(xml_out, "w") as f:
    f.write(XmlSerializer.serialize(system))
