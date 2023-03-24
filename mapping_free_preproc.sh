#!/bin/bash

ref_dir="/home/bioinformatikai/HW1/refs"
input_dir="/home/bioinformatikai/HW1/inputs"
output_dir="/home/bioinformatikai/HW1/outputs"
genome="mm10"
gtf_file="${ref_dir}/${genome}.gtf.gz"
transcripts_file="${ref_dir}/${genome}_rna.fa.gz"
index_dir="${ref_dir}/${genome}_salmon_index"

# Create index if it does not exist
if [ ! -d "${index_dir}" ]; then
    salmon index -t "${transcripts_file}" -i "${index_dir}" -p 6
fi

## Data QC
mkdir -p "${output_dir}/fastqc"
fastqc -o "${output_dir}/fastqc" "${input_dir}"/*.fastq.gz
multiqc "${output_dir}/fastqc" -o "${output_dir}"

mkdir -p "${output_dir}/trimming"
for file in "${input_dir}"/*_1.fastq.gz; do
    sample=$(basename "${file}" _1.fastq.gz)
    trim_galore -j 6 --paired --length 20 -o "${output_dir}/trimming" "${input_dir}/${sample}_1.fastq.gz" "${input_dir}/${sample}_2.fastq.gz"
done

mkdir -p "${output_dir}/fastqc_trimmed"
fastqc -o "${output_dir}/fastqc_trimmed" "${output_dir}/trimming"/*_val_*.fq.gz

## Quantification
mkdir -p "${output_dir}/quant"
for file in "${output_dir}/trimming"/*_1_val_1.fq.gz; do
    sample=$(basename "${file}" _1_val_1.fq.gz)
    salmon quant -i "${index_dir}" -l A -1 "${file}" -2 "${file/_1_val_1/_2_val_2}" -o "${output_dir}/quant/${sample}" --validateMappings
done