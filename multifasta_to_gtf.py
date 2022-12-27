from Bio import SeqIO
from collections import defaultdict
import argparse

'''
Create generic GTF records and new fasta from a fasta/multifasta, with the intention
of them being appended to a reference. Inspiried by cellranger's "Add a
marker gene to the FASTA and GTF" docs:

https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/tutorial_mr#marker
'''

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input",
                    help="Input fasta or multifasta")
parser.add_argument("-g", "--gtf",
                    help="Output GTF")
parser.add_argument("-f", "--fasta",
                    help="Output fasta/multifasta")
args = parser.parse_args()


with open(args.input) as f:

    # Will use this to count duplicate gene names. Default value is 0.
    seen = defaultdict(int)


    for fasta_record in SeqIO.parse(f, "fasta"):
        # Parse fasta/multifasta (my ids are just the gene names)
        seq = fasta_record.seq
        length = len(seq)
        id = fasta_record.id

        # Add gene to dictionary (if not in there already) and increase counter
        seen[id] +=1

        # If there are multiples of this gene, append ".2" etc to name (e.g. TRAV5.2)
        if seen[id] > 1:
            gene = id + '.' + str(seen[id])
        else:
            gene = id

        # Generate gtf record
        gtf_record = '{gene}\tunknown\texon\t1\t{length}\t.\t+\t.\t' \
                     'gene_id "{gene}"; transcript_id "{gene}"; ' \
                     'gene_name "{gene}"; ' \
                     'gene_biotype "protein_coding";\n'.format(gene = gene,
                                                               length = length)
        
        # Append to output gtf
        with open(args.gtf, 'a') as g:
            g.write(gtf_record)

        # Output new multifasta with updated (e.g. TRAV5.2) gene names
        with open(args.fasta, 'a') as f:
            f.write('>' + gene + '\n')
            f.write(str(seq) + '\n')
