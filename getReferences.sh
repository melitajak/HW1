#!/bin/bash

# download the reference genome in FASTA format 
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/GRCm38.primary_assembly.genome.fa.gz -O /home/bioinformatikai/HW1/refs/mm10.fa.gz

# download the reference transcriptome in FASTA format 
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/gencode.vM25.transcripts.fa.gz -O /home/bioinformatikai/HW1/refs/mm10_rna.fa.gz

# download the GTF file for the reference genome
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/gencode.vM25.primary_assembly.annotation.gtf.gz -O /home/bioinformatikai/HW1/refs/mm10.gtf.gz


prefetch -O /home/bioinformatikai/HW1/inputs SRR8985047 SRR8985048 SRR8985051 SRR8985052

fastq-dump --outdir /home/bioinformatikai/HW1/inputs/ --gzip /home/bioinformatikai/HW1/inputs/SRR8985047 /home/bioinformatikai/HW1/inputs/SRR8985048 /home/bioinformatikai/HW1/inputs/SRR8985051 /home/bioinformatikai/HW1/inputs/SRR8985052
