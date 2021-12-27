#!/usr/bin/env nextflow

params.strain = "DBA_2J"
params.input_dir = "./data/input"
params.input_files = "${params.input_dir}/${params.strain}-*.vcf"

params.mappings = "./data/input/mappings.txt"

Channel.fromPath(params.input_files).set{input_files}

process bed_from_vcf {

    input:
        file vcf_file from input_files
        file mappings_file from file(params.mappings)

    output:
        file "${params.strain}-*.bed" into vcf_files
    
    
    script:

    caller = vcf_file.getName().toString().tokenize('-').get(1).tokenize('.').get(0)
    
    """
    while read -r line
    do
        TYPE="\$(echo \$line | cut -d ':' -f 1)"
        MAPPING="\$(echo \$line | cut -d ':' -f 2)"
        
        /home/egarcia/appdir/bcftools/bin/bcftools query -i"\$MAPPING" -f'%CHROM\\t%POS0\\t%END0\\t%SVLEN\\n' ${vcf_file} | awk -F'\\t' 'BEGIN {OFS = FS} \$1 ~/^[0-9]*\$|^X\$/{print \$1,\$2,\$3,\$4}' >> "${params.strain}-${caller}-\$TYPE.bed"
    done < $mappings_file
    """
}


vcf_files.view()