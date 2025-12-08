from openmm.app import *
from openmm import *
from openmm.unit import *
import numpy as np

xml = snakemake.input["xml"]
pdbfile = snakemake.input["pdb"]

forces_out = snakemake.output["forces"]
positions_out = snakemake.output["positions"]
topo_out = snakemake.output["topo"]

# Load system
system = XmlSerializer.deserialize(open(xml).read())
pdb = PDBFile(pdbfile)

integrator = LangevinIntegrator(300*kelvin, 1/picosecond, 0.002*picoseconds)
platform = Platform.getPlatformByName("CPU")

simulation = Simulation(pdb.topology, system, integrator, platform)
simulation.context.setPositions(pdb.positions)

# Compute instantaneous forces
state = simulation.context.getState(getPositions=True, getForces=True)

forces = state.getForces(asNumpy=True).value_in_unit(kilojoules_per_mole/nanometer)
positions = state.getPositions(asNumpy=True).value_in_unit(nanometer)

np.save(forces_out, forces)
np.save(positions_out, positions)

# Save coordinates
with open(topo_out, "w") as f:
    PDBFile.writeFile(pdb.topology, state.getPositions(), f)
