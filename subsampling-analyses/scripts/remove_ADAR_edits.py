import argparse
from Bio import AlignIO, SeqIO

# From Dudas et al, 2018:
#
# As noticed by Tong et al.38, Park et al.13 and other studies, some EBOV isolates contain clusters of T-to-C mutations within
# relatively short stretches of the genome. Interferon-inducible adenosine deaminases acting on RNA (ADAR) are known to induce
# adenosine to inosine hypermutations in double-stranded RNA42. ADARs have been suggested to act on RNAs from numerous groups
# of viruses43. When negative sense single stranded RNA virus genomes are edited by ADARs, A-to-G hypermutations seem to
# preferentially occur on the negative strand, which results in U/T-to-C mutations on the positive strand44–46. Multiple
# T-to-C mutations are introduced simultaneously via ADAR-mediated RNA editing which would interfere with molecular clock
# estimates and, by extension, the tree topology. We thus designate four or more T-to-C mutations within 300 nucleotides of
# each other as a putative hypermutation tract, whenever there is evidence that all T-to-C mutations within such stretches
# were introduced at the same time, i.e. every T-to-C mutation in a stretch occurred on a single branch. We detect a total
# of 15 hypermutation patterns with up to 13 T-to-C mutations within 35 to 145 nucleotides. Of these patterns, 11 are unique
# to a single genome and 4 are shared across multiple isolates, suggesting that occasionally viruses survive hypermutation
# are transmitted47. Putative tracts of T-to-C hypermutation almost exclusively occur within non-coding intergenic regions,
# where their effects on viral fitness are presumably minimal. In each case we mask out these sites as ambiguous nucleotides
# but leave the first T-to-C mutation unmasked to provide phylogenetic information on the relatedness of these sequences.


def collapse(window, n, potential, seq):
    # look backwards in a greedy fashion and see if there's a stretch within the allowed window
    adar = []
    open_idx = len(potential)-1
    while True:
        # we're done once we get so close to the start of the potentials that there's not room for enough
        if open_idx < n-1:
            break
        for close_idx in range(open_idx-1,-1,-1):
            if potential[open_idx][0]-potential[close_idx][0] > window:
                stretch = [potential[idx] for idx in range(close_idx+1, open_idx+1)]
                if len(stretch) > n:
                    adar.append(stretch)
                open_idx = close_idx
                break

            if close_idx == 0:
                stretch = [potential[idx] for idx in range(0, open_idx+1)]
                if len(stretch) > n:
                    adar.append(stretch)
                open_idx = 0
                break

    return adar


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("--input", help="Alignment to detect & remove ADAR edits")
    parser.add_argument("--output", help="File to write ADAR removed alignment")
    parser.add_argument("--onlyTC", action="store_true", help="Don't consider AG mutations")
    parser.add_argument("--window", default=300, help="window length")
    parser.add_argument("--n", default=4, help="n ADAR mutations in window required to call")
    args = parser.parse_args()

    aln = AlignIO.read(args.input, "fasta")
    ref = [a for a in aln if a.id=='MK007329.1'][0]

    for seq in aln:
        potential_AG = []
        potential_TC = []

        assert(len(seq) == len(ref))
        for idx in range(0, len(ref)):
            wt = ref[idx]
            alt = seq[idx]
            if wt=="T" and alt=="C":
                potential_TC.append([idx, "TC"])
            elif wt=="A" and alt=="G" and not args.onlyTC:
                potential_AG.append([idx, "AG"])

        adar_TC = collapse(args.window, args.n, potential_TC, seq)
        adar_AG = collapse(args.window, args.n, potential_AG, seq)
        
        if len(adar_TC) or len(adar_AG):
            print(f"\n{seq.id} detected {len(adar_TC)} region(s) of T->C and {len(adar_AG)} region(s) of A->G ADAR hypermutation")
            mutable_seq = seq.seq.tomutable()
            for stretch in adar_TC:
                print(f"\tT->C @ positions {','.join([str(x[0]) for x in stretch])} (turning all but first into Ns)")
                for x in stretch[1:]:
                    mutable_seq[x[0]]="N"
            for stretch in adar_AG:
                print(f"\tA->G @ positions {','.join([str(x[0]) for x in stretch])} (turning all but first into Ns)")
                for x in stretch[1:]:
                    mutable_seq[x[0]]="N"
            seq.seq = mutable_seq.toseq()

    SeqIO.write(aln, args.output, "fasta")