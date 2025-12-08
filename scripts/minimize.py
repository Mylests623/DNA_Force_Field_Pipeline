from openmm import *
from openmm.app import *
from openmm.unit import *

# Load system
system = XmlSerializer.deserialize(open("output/system.xml").read())

# Load coordinates
pdb = PDBFile("output/solvated.pdb")

integrator = LangevinIntegrator(300*kelvin, 1/picosecond, 0.002*picoseconds)
platform = Platform.getPlatformByName("CPU")

simulation = Simulation(pdb.topology, system, integrator, platform)
simulation.context.setPositions(pdb.positions)

# Before minimization
state = simulation.context.getState(getEnergy=True)
print("Initial potential:", state.getPotentialEnergy())

# Minimize
simulation.minimizeEnergy(maxIterations=200)

state = simulation.context.getState(getEnergy=True)
print("After minimization:", state.getPotentialEnergy())
