# def get_todays_date():
#     from datetime import datetime
#     date = datetime.today().strftime('%Y-%m-%d')
#     return date

rule all:
    input:
        auspice = "auspice/ebola_nord-kivu_manuscript.json"

rule files:
    params:
        manually_aligned_sequences = "data/manual_alignment.fasta",
        metadata = "data/metadata.tsv",
        dropped_strains = "config/dropped_strains.txt",
        reference = "config/reference.gb",
        colors = "config/colors.tsv",
        lat_longs = "config/lat_longs.tsv",
        auspice_config = "config/auspice_config.json",
        root_name = "MK007329.1",
        clades = "config/clades.tsv"

files = rules.files.params

rule filter:
    message:
        """
        Filtering manual alignment to exclude explicitly dropped strains from {input.exclude}, including the reference!
        """
    input:
        sequences = files.manually_aligned_sequences,
        metadata = files.metadata,
        exclude = files.dropped_strains
    output:
        sequences = "results/filtered.fasta"
    shell:
        """
        augur filter \
            --sequences {input.sequences} \
            --metadata {input.metadata} \
            --exclude {input.exclude} \
            --output {output.sequences} \
        """

rule remove_potential_sequencing_errors:
    message:
        """
        Removing SNPs which flank regions of Ns as they are probably be sequencing errors
        """
    input:
        sequences = rules.filter.output.sequences
    output:
        sequences = "results/filtered_errors_removed.fasta"
    shell:
        """
        python scripts/remove_SNPs_flanking_Ns.py --input {input.sequences} --output {output.sequences}
        """

rule remove_potential_adar_mutations:
    message:
        """
        Removing SNPs which look like ADAR mutations
        """
    input:
        sequences = rules.remove_potential_sequencing_errors.output.sequences
    output:
        sequences = "results/filtered_errors_removed_no_adar.fasta"
    shell:
        """
        python scripts/remove_ADAR_edits.py --input {input.sequences} --output {output.sequences}
        """

rule tree:
    message: "Building tree"
    input:
        alignment = rules.remove_potential_adar_mutations.output.sequences
    output:
        tree = "results/tree_raw.nwk"
    shell:
        """
        augur tree \
            --alignment {input.alignment} \
            --output {output.tree} \
            --nthreads auto
        """

rule reroot:
    message:
        """
        Rerooting tree based on reference sequence
        """
    input:
        tree = rules.tree.output.tree,
    output:
        tree = "results/tree_rerooted.nwk",
    params:
        root = files.root_name
    shell:
        """
        augur refine \
            --tree {input.tree} \
            --output-tree {output.tree} \
            --root {params.root}
        """

rule prune_outgroup:
    message: "Pruning the outgroup (reference sequence MK007329.1) from the tree"
    input:
        tree = rules.reroot.output.tree
    output:
        tree = "results/tree_rerooted_pruned.nwk"
    params:
        root = files.root_name
    run:
        from Bio import Phylo
        T = Phylo.read(input[0], "newick")
        outgroup = [c for c in T.find_clades() if str(c.name) == params[0]][0]
        T.prune(outgroup)
        T.ladderize()
        Phylo.write(T, output[0], "newick")

rule refine:
    message:
        """
        Refining tree
          - estimate timetree
          - use {params.coalescent} coalescent timescale
          - estimate {params.date_inference} node dates
        """
    input:
        tree = rules.prune_outgroup.output.tree,
        alignment = rules.remove_potential_adar_mutations.output.sequences,
        metadata = files.metadata,
    output:
        tree = "results/tree.nwk",
        node_data = "results/branch_lengths.json"
    params:
        coalescent = "skyline",
        date_inference = "marginal"
    shell:
        """
        augur refine \
            --tree {input.tree} \
            --alignment {input.alignment} \
            --metadata {input.metadata} \
            --output-tree {output.tree} \
            --output-node-data {output.node_data} \
            --timetree \
            --keep-root \
            --coalescent {params.coalescent} \
            --date-confidence \
            --date-inference {params.date_inference} \
            --keep-polytomies
        """

## This rule shouldn't be necessary, but I think there's a bug when using Treetime's `keep-root` functionality
rule ladderize:
    message: "Ladderizing the refined tree"
    input:
        tree = rules.refine.output.tree
    output:
        tree = "results/tree_ladderized.nwk"
    run:
        from Bio import Phylo
        T = Phylo.read(input[0], "newick")
        T.ladderize()
        Phylo.write(T, output[0], "newick")


rule ancestral:
    message: "Reconstructing ancestral sequences and mutations"
    input:
        tree = rules.refine.output.tree,
        alignment = rules.remove_potential_adar_mutations.output.sequences,
    output:
        node_data = "results/nt_muts.json"
    params:
        inference = "joint"
    shell:
        """
        augur ancestral \
            --tree {input.tree} \
            --alignment {input.alignment} \
            --output-node-data {output.node_data} \
            --inference {params.inference}
        """

rule translate:
    message: "Translating amino acid sequences"
    input:
        tree = rules.refine.output.tree,
        node_data = rules.ancestral.output.node_data,
        reference = files.reference
    output:
        node_data = "results/aa_muts.json"
    shell:
        """
        augur translate \
            --tree {input.tree} \
            --ancestral-sequences {input.node_data} \
            --reference-sequence {input.reference} \
            --output {output.node_data} \
        """
rule clades:
    message: " Labeling clades as specified in config/clades.tsv"
    input:
        tree = rules.refine.output.tree,
        aa_muts = rules.translate.output.node_data,
        nuc_muts = rules.ancestral.output.node_data,
        clades = files.clades
    output:
        clade_data = "results/clades.json"
    shell:
        """
        augur clades --tree {input.tree} \
            --mutations {input.nuc_muts} {input.aa_muts} \
            --clades {input.clades} \
            --output {output.clade_data}
        """

rule traits:
    message: "Inferring ancestral traits for {params.columns!s}"
    input:
        tree = rules.refine.output.tree,
        metadata = files.metadata
    output:
        node_data = "results/traits.json",
    params:
        columns = "health_zone"
    shell:
        """
        augur traits \
            --tree {input.tree} \
            --metadata {input.metadata} \
            --output {output.node_data} \
            --columns {params.columns} \
            --confidence
        """

rule export:
    message: "Exporting data files for auspice"
    input:
        tree = rules.ladderize.output.tree,
        metadata = files.metadata,
        branch_lengths = rules.refine.output.node_data,
        traits = rules.traits.output.node_data,
        nt_muts = rules.ancestral.output.node_data,
        aa_muts = rules.translate.output.node_data,
        colors = files.colors,
        lat_longs = files.lat_longs,
        clades = rules.clades.output.clade_data,
        auspice_config = files.auspice_config,
    output:
        auspice = rules.all.input.auspice
    shell:
        """
        augur export v2 \
            --tree {input.tree} \
            --metadata {input.metadata} \
            --node-data {input.branch_lengths} {input.traits} {input.nt_muts} {input.clades} {input.aa_muts} \
            --colors {input.colors} \
            --lat-longs {input.lat_longs} \
            --auspice-config {input.auspice_config} \
            --output {output.auspice}
        """

rule clean:
    message: "Removing directories: {params}"
    params:
        "results ",
        "auspice"
    shell:
        "rm -rfv {params}"
