# Operationalizing genomic epidemiology during the Nord-Kivu Ebola outbreak, Democratic Republic of the Congo.

This repository represents the analysis and code behind Kinganda-Lusamaki _et al._, currently available at https://www.medrxiv.org/content/10.1101/2020.06.08.20125567v1


# Figures

Analysis in the paper was performed using Jupyter notebooks and Nextstrain.
I suggest installing the [Nextstrain conda environment](https://nextstrain.org/docs/getting-started/local-installation#install-augur--auspice-with-conda-recommended) as well as jupyter (e.g. `pip install jupyterlab`).
Full resolution figures can be found in `./figures`.

**Figure 1**, "Progress of genomic surveillance over the course of the outbreak", is produced by `./notebooks/lags.ipynb`.

**Figure 2**, "Broad scale spatiotemporal dynamics of the EVD outbreak in Nord Kivu province DRC", comprises screenshots from Auspice of the `./auspice/ebola-narrative-ms_full-build.json` dataset.

**Figure 3**, "Transmission patterns within and between health zones", as well as corresponding supplementary figures, are produced by `./notebooks/chains.ipynb`

**Figure 4**, "Initial genomic evidence for an infection relapse event" was created similarly to Figure 2. 

**Supplementary Figures 2 & 3**, "Patterns of transmission between health zones (HZ)" and "Viral movements are robust to subsampling" are produced by the same notebook as Figure 3 (`./notebooks/chains.ipynb`).

**Supplementary Figure 4**, "The outbreak was sustained by movement between health zones", is produced by `./notebooks/exploded_tree.ipynb`.

**Supplementary Figures 5 & 6**, "Genomic characterization of transmission after unsafe burial of a pastor" and "Secondary transmission associated with infection of a motorcycle taxi driver", are produced from Auspice screenshots of the `./auspice/ebola-narrative-ms_full-build.json` dataset with additional manual annotations.

# Extended methods

### Generation of manual alignment `data/manual_alignment.fasta`

The sequence alignment was generated using a banded-alignment approach with manual editing in [AliView](http://ormbunkar.se/aliview/). The dataset was partitioned using a script from [nextstrain/ncov](https://github.com/nextstrain/ncov/blob/master/scripts/partition-sequences.py) (`python partition-sequences.py --sequences data/sequences.fasta --sequences-per-group 102 --output-dir ignore/bands/`). Alignment was performed using a modified version of `augur` which added the `--maxiterate 1000` argument to [mafft](https://mafft.cbrc.jp/alignment/software/) (`augur align --nthreads 2 --sequences results/0.fasta --reference-sequence config/reference.gb --output results/0.aligned.fasta --fill-gaps`).
Each of the resulting alignments were checked by hand in AliView, with no clear misalignments identified.
Alignments were combined using a script from [nextstrain/ncov](https://github.com/nextstrain/ncov/blob/master/scripts/partition-sequences.py) (`python combine-and-dedup-fastas.py --input ignore/bands/*.aligned.fasta --output data/aln.fasta`).
ADAR edits, which were noticed upon manual inspection were removed using a custom script based upon [Dudas et al](https://www.nature.com/articles/nature22040) via `python scripts/remove_ADAR_edits.py --input data/aln.fasta --output results/aln.noAdar.fasta`.
A number of (singlteon) SNPs flanking regions of `N`s were also detected upon manual checking of the alignment which are most probably spurious sequencing errors.
These were removed via `scripts/remove_SNPs_flanking_Ns.py --input aln.noAdar.fasta --output results/manual_alignment.fasta`.
The final alignment consisted of 808 sequences.

