import os

pdb_file = "galP1_delP2_RNAP_CRP_PRc_2nd"
INPUT = f"input/{pdb_file}.pdb"


rule all:
    input:
        # f"output/{pdb_file}_system.xml",
        # f"output/{pdb_file}_state_forces.npy",
        # f"output/{pdb_file}_state_positions.npy",
        # f"output/{pdb_file}_state_topology.pdb",
        # f"output/{pdb_file}_dna_forces_per_residue.tsv",
        # f"output/{pdb_file}_dna_forces_bfactor.pdb",
        # f"output/{pdb_file}_full_complex_with_dna_forces.pdb"
        f"output/{pdb_file}_reference_system.xml",
        f"output/{pdb_file}_reference_forces.npy",
        f"output/{pdb_file}_reference_positions.npy",
        f"output/{pdb_file}_reference_topology.pdb",
        f"output/{pdb_file}_reference_forces_per_residue.tsv",
        f"output/{pdb_file}_reference_forces_bfactor.pdb"


rule validate_pdb:
    input:
        INPUT
    output:
        "output/validated.pdb"
    script:
        "scripts/validate_pdb.py"

rule extract_dna:
    input:
        "output/validated.pdb"
    output:
        "output/dna_only.pdb"
    script:
        "scripts/extract_dna.py"

rule extract_dna_sequence:
    input:
        "output/dna_only.pdb"
    output:
        "output/dna_sequence.fasta"
    script:
        "scripts/extract_sequence_from_pdb.py"

rule build_linear_dna:
    input:
        "output/dna_sequence.fasta"
    output:
        "output/reference_linear.pdb"
    conda:
        "envs/3dna.yaml"
    shell:
        """
        # Using 3DNA ‘rebuild’ tool; assumes installed in PATH
        x3dna-dssr --rebuild --seqfile {input} --o {output} --model=bdna
        """

rule build_reference_system:
    input:
        "output/reference_linear.pdb"
    output:
        xml=f"output/{pdb_file}_reference_system.xml",
        pdb=f"output/{pdb_file}_reference_solvated.pdb"
    script:
        "scripts/build_dna_system.py"

rule compute_reference_forces:
    input:
        xml=f"output/{pdb_file}_reference_system.xml",
        pdb=f"output/{pdb_file}_reference_solvated.pdb"
    output:
        forces=f"output/{pdb_file}_reference_forces.npy",
        positions=f"output/{pdb_file}_reference_positions.npy",
        topo=f"output/{pdb_file}_reference_topology.pdb"
    script:
        "scripts/compute_dna_forces.py"

rule aggregate_reference_forces:
    input:
        forces=f"output/{pdb_file}_reference_forces.npy",
        topo=f"output/{pdb_file}_reference_topology.pdb"
    output:
        tsv=f"output/{pdb_file}_reference_forces_per_residue.tsv",
        pdb=f"output/{pdb_file}_reference_forces_bfactor.pdb"
    script:
        "scripts/aggregate_forces.py"

rule build_system:
    input:
        "output/dna_only.pdb"
    output:
        xml=f"output/{pdb_file}_system.xml",
        pdb="output/solvated.pdb"
    script:
        "scripts/build_dna_system.py"

rule compute_forces:
    input:
        xml=f"output/{pdb_file}_system.xml",
        pdb="output/solvated.pdb"
    output:
        forces=f"output/{pdb_file}_state_forces.npy",
        positions=f"output/{pdb_file}_state_positions.npy",
        topo=f"output/{pdb_file}_state_topology.pdb"
    script:
        "scripts/compute_dna_forces.py"

rule aggregate_forces:
    input:
        forces=f"output/{pdb_file}_state_forces.npy",
        topo=f"output/{pdb_file}_state_topology.pdb"
    output:
        tsv=f"output/{pdb_file}_dna_forces_per_residue.tsv",
        bpfb_pdb=f"output/{pdb_file}_dna_forces_bfactor.pdb"
    script:
        "scripts/aggregate_forces.py"

# rule pymol_heatmap:
#     input:
#         pdb="output/dna_forces_bfactor.pdb"
#     output:
#         pml="output/pymol_heatmap.pml",
#         png="output/dna_force_heatmap.png"
#     params:
#         # optional: user can tune PyMOL render size and camera options
#         width=1600,
#         height=1200,
#         dpi=150
#     run:
#         # create a small pml script (with computed min/max embedded) then call pymol to render
#         shell("python3 scripts/write_pymol_pml.py {input.pdb} {output.pml} {output.png} {params.width} {params.height} {params.dpi}")
rule map_forces_to_full_complex:
    input:
        full= INPUT,
        dna=f"output/{pdb_file}_state_topology.pdb",
        tsv=f"output/{pdb_file}_dna_forces_per_residue.tsv"
    output:
        pdb=f"output/{pdb_file}_full_complex_with_dna_forces.pdb"
    script:
        "scripts/map_forces_to_full_complex.py"
