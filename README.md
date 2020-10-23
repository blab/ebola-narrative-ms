# Operationalizing genomic epidemiology during the Nord-Kivu Ebola outbreak, Democratic Republic of the Congo.

This repository represents the analysis and code behind Kinganda-Lusamaki _et al._, currently available at https://www.medrxiv.org/content/10.1101/2020.06.08.20125567v1


# Jupyter notebooks

The analysis in the paper was performed using Jupyter notebooks.
The analysis was performed in a [Nextstrain conda environment](https://nextstrain.org/docs/getting-started/local-installation#install-augur--auspice-with-conda-recommended) which will give you the necessary python libraries.

Figure 1 is produced by `notebooks/lags.ipynb`.


# Extended methods

### Generation of manual alignment `data/manual_alignment.fasta`

The sequence alignment was generated using a banded-alignment approach with manual editing in [AliView](http://ormbunkar.se/aliview/). The dataset was partitioned using a script from [nextstrain/ncov](https://github.com/nextstrain/ncov/blob/master/scripts/partition-sequences.py) (`python partition-sequences.py --sequences data/sequences.fasta --sequences-per-group 102 --output-dir ignore/bands/`). Alignment was performed using a modified version of `augur` which added the `--maxiterate 1000` argument to [mafft](https://mafft.cbrc.jp/alignment/software/) (`augur align --nthreads 2 --sequences results/0.fasta --reference-sequence config/reference.gb --output results/0.aligned.fasta --fill-gaps`).
Each of the resulting alignments were checked by hand in AliView, with no clear misalignments identified.
Alignments were combined using a script from [nextstrain/ncov](https://github.com/nextstrain/ncov/blob/master/scripts/partition-sequences.py) (`python combine-and-dedup-fastas.py --input ignore/bands/*.aligned.fasta --output data/aln.fasta`).
ADAR edits, which were noticed upon manual inspection were removed using a custom script based upon [Dudas et al](https://www.nature.com/articles/nature22040) via `python scripts/remove_ADAR_edits.py --input data/aln.fasta --output results/aln.noAdar.fasta`.
A number of (singlteon) SNPs flanking regions of `N`s were also detected upon manual checking of the alignment which are most probably spurious sequencing errors.
These were removed via `scripts/remove_SNPs_flanking_Ns.py --input aln.noAdar.fasta --output results/manual_alignment.fasta`.
The final alignment consisted of 808 sequences.

