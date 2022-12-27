from Bio import SeqIO

'''
Create generic GTF records from a multifasta, with the intention
of them being appended to a reference GTF. Inspiried by cellranger's "Add a 
marker gene to the FASTA and GTF" docs:

https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/tutorial_mr#marker
'''

in_multifasta = "path/to/in.fasta"  # Input fasta or multifasta
out_gtf = "path/to/out.gtf" # Desired output file

with open(in_multifasta) as f:

    for fasta_record in SeqIO.parse(f, "fasta"):
        gene = fasta_record.id  # My fasta ids are just the gene names
        length = len(fasta_record.seq)

        # GTF record
        outline = '{gene}\tunknown\texon\t1\t{length}\t.\t+\t.\tgene_id "{gene}"; ' \
               'transcript_id "{gene}"; gene_name "{gene}"; ' \
               'gene_biotype "protein_coding";\n'.format(gene = gene, 
                                                         length = length)
        
        # Append to output GTF
        with open(out_gtf, 'a') as outfile:
            outfile.write(outline)
