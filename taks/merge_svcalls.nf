#!/usr/bin/env nextflow

params.strain = "DBA_2J"
params.input_dir = "./data/input"
params.input_files = "${params.input_dir}/${params.strain}-*.vcf"

params.mappings = "./data/input/mappings.txt"

Channel.fromPath(params.input_files).set{input_files}

process mapped_bed_from_vcf {

    input:
        file vcf_file from input_files
        file mappings_file from file(params.mappings)

    output:
        file "${params.strain}-*.bed" into bed_files
    
    
    //Generate a set of BED files by type according to the values in data/input/mappings.txt 
    script:

    caller = vcf_file.getName().toString().tokenize('-').get(1).tokenize('.').get(0)
    
    """
    while read -r line
    do
        TYPE="\$(echo \$line | cut -d ':' -f 1)"
        MAPPING="\$(echo \$line | cut -d ':' -f 2)"
        
        /home/egarcia/appdir/bcftools/bin/bcftools query -i"\$MAPPING" -f'%CHROM\\t%POS0\\t%END0\\t%SVLEN\\n' ${vcf_file} | \
        awk -F'\\t' 'BEGIN {OFS = FS} \$1 ~/^[0-9]*\$|^X\$/{print \$1,\$2,\$3,\$4}' >> "${params.strain}-${caller}-\$TYPE.bed"
    done < $mappings_file
    """
}


process mapped_vcf_from_bed {

    input:
        file bed_file from bed_files.flatten()
    
    output:
        file '*.vcf' into vcf_files

    script:

    name_parts = bed_file.getName().toString().tokenize('.').get(0).tokenize('-')
    strain = name_parts.get(0)
    caller = name_parts.get(1)
    type = name_parts.get(2)

    """
        /home/egarcia/workspace/github/SURVIVOR/Debug/SURVIVOR bedtovcf ${bed_file} ${type} '${strain}-${caller}-${type}.vcf'
    """

}