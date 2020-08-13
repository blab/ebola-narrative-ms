import argparse
from Bio import AlignIO, SeqIO

# On inspection of manual alignment I noticed a number of SNPs flanking regions of Ns. These SNPs were singletons
# (not seen in other samples in the alignment) and are most probably spurious. This script will remove them. The expected
# change in the tree will only be terminal branch lengths since they don't seem to be transmitted, but this will improve
# the overall temporal inference


def check_flanks(seq, ref, nstart, nstop):
    # start and stop idxs are Ns
    # print(seq.id, nstart, nstop, nstop+1-nstart)
    x = []
    if nstart-1 >=0 and seq[nstart-1]!=ref[nstart-1]:
        x.append(nstart-1)
        # print(seq.id, "SNP BEFORE",ref[nstart-1], seq[nstart-1])
    if (nstop+1)<len(seq) and seq[nstop+1]!=ref[nstop+1]:
        # print(seq.id, "SNP AFTER",ref[nstop+1], seq[nstop+1])
        x.append(nstop+1)
    return x

def change_to_N(seq, posns):
    mutable_seq = seq.seq.tomutable()
    for p in posns:
        mutable_seq[p] = 'N'
    seq.seq = mutable_seq.toseq()
    


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("--input", help="Alignment to detect & remove ADAR edits")
    parser.add_argument("--output", help="File to write ADAR removed alignment")
    parser.add_argument("--window", default=2, help="min number of Ns to check either side of")
    # parser.add_argument("--look", default=1, help="number of bases either side of Ns to consider")
    args = parser.parse_args()

    aln = AlignIO.read(args.input, "fasta")
    ref = [a for a in aln if a.id=='MK007329.1'][0]
    n_changes = 0


    # for seq in [a for a in aln if a.id=="GOM2467" or a.id=="MAN9412"]:
    for seq in aln:
        inside = False
        start = 0
        sequencing_error_positions = []
        for pos in range(0, len(ref)):
            # open stretch of Ns
            if not inside and seq[pos]=="N":
                inside = True
                start = pos
            if inside and seq[pos]!="N":
                inside = False
                if pos-start >= args.window:
                    sequencing_error_positions += check_flanks(seq, ref, start, pos-1)
        if inside and (len(ref)-1)-start >= args.window:
            sequencing_error_positions += check_flanks(seq, ref, start, len(ref)-1)
        if len(sequencing_error_positions):
            print(f"{seq.id} had {len(sequencing_error_positions)} SNP(s) flanking regions of Ns. Changing them to Ns.")
            change_to_N(seq, sequencing_error_positions)
            n_changes += len(sequencing_error_positions)
    
    print(f"--- Overall there were {n_changes} bases changed to Ns ---")
    SeqIO.write(aln, args.output, "fasta")

