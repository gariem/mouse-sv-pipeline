#!/usr/bin/env nextflow

params.file= "DBA_2J-mapped_200_2_20_all.vcf"
params.input_dir = "./data/merged"
params.input_files = "${params.input_dir}/${params.file}"
params.out_dir = './data/reports/raw'
params.previous_dir = './data/previous'
params.validated_dir = './data/validated'

params.merge_mappings = "./data/input/merge_mappings.txt"
params.metrics_mappings = "./data/input/metrics_mappings.txt"
params.genome_sizes = "./data/input/genome_size.txt"


Channel.fromPath(params.input_files).set{input_files}

process generate_data_file {

    input:
        file vcf_file from input_files

    output:
        file 'data.csv' into data_file

    """
    bcftools query -f'%CHROM,%SVTYPE,%SVLEN\\n' ${vcf_file} | awk -F',' 'BEGIN {OFS = FS} {abs=\$3;sub("^-", "",abs); print \$1,\$2,abs}' > "data.csv"
    """
}

