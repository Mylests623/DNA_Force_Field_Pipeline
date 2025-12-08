import matplotlib.pyplot as plt

# Load per-residue forces
residues = []
forces = []

with open("dna_forces_per_residue.txt") as f:
    next(f)  # skip header
    for line in f:
        res, force = line.strip().split("\t")
        residues.append(res)
        forces.append(float(force))

plt.figure(figsize=(12,4))
plt.bar(residues, forces, color='firebrick')
plt.xticks(rotation=90)
plt.ylabel("Average Force (pN)")
plt.title("DNA Stress / Force Heat Map")
plt.tight_layout()
plt.show()
