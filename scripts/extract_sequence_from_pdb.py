from Bio.PDB import PDBParser
import sys

pdb_path = snakemake.input[0]
out = snakemake.output[0]

parser = PDBParser(QUIET=True)
structure = parser.get_structure("dna", pdb_path)

seq = []
three_to_one = {
    "DA":"A","DT":"T","DG":"G","DC":"C",
    "A":"A","T":"T","G":"G","C":"C"
}

for model in structure:
    for chain in model:
        s = ""
        for res in chain:
            name = res.get_resname().strip()
            if name in three_to_one:
                s += three_to_one[name]
        if s:
            seq.append((chain.id, s))

with open(out, "w") as f:
    for cid, s in seq:
        f.write(f">{cid}\n{s}\n")
