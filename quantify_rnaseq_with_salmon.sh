#!/bin/bash

# Define genome fasta file and genome annotation file
genome_fa="genome.fa"
annotation_gtf="annotation.gtf"

# Generate transcriptome fasta file
gffread -w transcriptome.fa -g $genome_fa $annotation_gtf

# Compress transcriptome and genome fasta files
gzip -k transcriptome.fa
gzip -k $genome_fa

# Define chromosome names as decoy sequences
grep "^>" <(gunzip -c $genome_fa.gz) | cut -d " " -f 1 > decoys.txt

# Remove fasta formatting
sed -i.bak -e 's/>//g' decoys.txt

# Concatenate transcriptome and genome, transcriptome MUST come first
cat transcriptome.fa.gz $genome_fa.gz > gentrome.fa.gz

# Generate the Salmon index
salmon index -t gentrome.fa.gz -d decoys.txt -i salmon_index

# Loop through sequencing reads and perform Salmon quantification
# may be single end
for i in /path/to/fastq_files/*.fastq.gz; do
	name=$(basename ${i} .fastq.gz);
    salmon quant -l A -i salmon_index --seqBias --gcBias --posBias -r /path/to/fastq_files/${name}.fastq.gz -o salmon_quants/${name}
done


# or paired end:
#for i in /path/to/paired_fastq_files/*_1.fastq.gz; do
#	name=$(basename ${i} _1.fastq.gz);
#   salmon quant -l A -i salmon_index --seqBias --gcBias --posBias \
#    -1 /path/to/paired_fastq_files/${name}.fastq.gz -2 /path/to/paired_fastq_files/${name}_2.fastq.gz -o salmon_quants/${name}
#done
