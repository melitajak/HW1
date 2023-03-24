#!/bin/bash

threads=6
refs_dir="/home/bioinformatikai/HW1/refs"
genome_file="${refs_dir}/mm10.fa.gz"
genome_index="${refs_dir}/mm10"
gtf_file="${refs_dir}/mm10.gtf.gz"
rna_file="${refs_dir}/mm10_rna.fa.gz"
inputs_dir="/home/bioinformatikai/HW1/inputs"
outputs_dir="/home/bioinformatikai/HW1/outputs"
results_dir="/home/bioinformatikai/HW1/results"

# Reference genome preparation
if [ ! -f "${genome_index}.1.ht2" ]; then
    gunzip -c ${genome_file} | hisat2-build ${refs_dir}/mm10.fa ${genome_index}
fi

# Run FASTQC analysis 
mkdir -p ${outputs_dir}/raw_data
fastqc -t ${threads} ${inputs_dir}/*fastq.gz -o ${outputs_dir}/raw_data

# Generate MULTIQC 
multiqc -o ${outputs_dir}/raw_data ${outputs_dir}/raw_data

# Run standard FASTQ trimming and re-run FASTQC on cleaned FASTQ files
for i in ${inputs_dir}/*_1.fastq.gz;
do
  R1=${i};
  R2="${inputs_dir}/"$(basename ${i} _1.fastq.gz)"_2.fastq.gz";
  trim_galore -j ${threads} -o ${outputs_dir} --paired --length 20 ${R1} ${R2}
done

fastqc -t ${threads} ${outputs_dir}/*fq.gz -o ${outputs_dir}

# Create MultiQC plots for raw and processed data
multiqc -o ${outputs_dir} ${outputs_dir}

# Mapping, QC and quantification 

# Map each sample to the reference genome and remove unmapped reads and incorrectly mapped reads
for fq1 in ${outputs_dir}/*_1_val_1.fq.gz; do
  fq2=${outputs_dir}/$(basename ${fq1} _1_val_1.fq.gz)_2_val_2.fq.gz
  name=$(basename ${fq1} _1_val_1.fq.gz)
  hisat2 -p ${threads} -x ${genome_index} -1 ${fq1} -2 ${fq2} -S ${outputs_dir}/${name}.sam
  samtools view -@ ${threads} -b -F 4 -bS -o ${outputs_dir}/${name}.bam ${outputs_dir}/${name}.sam
  samtools sort -@ ${threads} -o ${outputs_dir}/${name}.sorted.bam ${outputs_dir}/${name}.bam
  samtools index -@ ${threads} ${outputs_dir}/${name}.sorted.bam
done

# Deduplicate BAM files and index them
for bam in ${outputs_dir}/*.sorted.bam; do
  name=$(basename ${bam} .sorted.bam)
  samtools sort -@ ${threads} -n -o ${outputs_dir}/${name}.namesort.bam ${bam}
  samtools fixmate -@ ${threads} -m ${outputs_dir}/${name}.namesort.bam ${outputs_dir}/${name}.fixmate.bam
  samtools sort -@ ${threads} -o ${outputs_dir}/${name}.positionsort.bam ${outputs_dir}/${name}.fixmate.bam
  samtools markdup -@ ${threads} -r ${outputs_dir}/${name}.positionsort.bam ${outputs_dir}/${name}.markdup.bam
  samtools index -@ ${threads} ${outputs_dir}/${name}.markdup.bam
done

# Run stringtie quantification on deduplicated BAM files
for bam in "${outputs_dir}"/*.markdup.bam; do
  name=$(basename "${bam}" .markdup.bam)
  stringtie -p ${threads} -G <(gunzip -c ${gtf_file}) -eB -o ${outputs_dir}/${name}.gtf -A ${outputs_dir}/${name}.tab ${bam}
done


# Create correlation and PCA plots using Deeptools
multiBamSummary bins --outFileName ${results_dir}/mapped.npz -p 6 --outRawCounts ${results_dir}/raw_counts.tsv -b ${outputs_dir}/*.markdup.bam
plotCorrelation -in ${results_dir}/mapped.npz -c pearson -p heatmap -o ${results_dir}/correlation_heatmap.pdf --removeOutliers
plotPCA -in ${results_dir}/mapped.npz -o ${results_dir}/pca_plot.pdf
