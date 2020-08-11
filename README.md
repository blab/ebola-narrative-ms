# ebola-narratives-ms



### Generation of manual alignment `data/manual_alignment.fasta`

I employed a banded-alignment approach with manual editing in AliView.
* Created partitioned sequences using `python ~/github/nextstrain/ncov/scripts/partition-sequences.py --sequences data/sequences.fasta --sequences-per-group 102 --output-dir ignore/bands/` 
* modified `augur align` to employ mafft using `--maxiterate 1000` for better restuls (sidenote: we should allow this as an option in `augur align`)
* `augur align --nthreads 2 --sequences results/0.fasta --reference-sequence config/reference.gb --output results/0.aligned.fasta --fill-gaps`
* checked each (n=7) alignment for obvious misalignments in AliView. There were none, however noticed three things, both of which should be fixed
    * Presence of ADAR edits
    * Presence of (single) SNPs flanking regions of Ns.
    * The outgroup is much more divergent than the reference. Should use the latter for outgroup rooting!
* Combined the 7 alignments via `python ~/github/nextstrain/ncov/scripts/combine-and-dedup-fastas.py --input ignore/bands/*.aligned.fasta --output data/manual_alignment.fasta`
* Manually renamed the reference to `MK007329.1` & removed the outgroup
* n=809 sequences (including reference + outgroup)


