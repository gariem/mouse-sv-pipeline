#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process size_data_from_vcf {

    input:
        file vcf_file
        
    output:
        file "*.csv"

    script:

    simple_name = vcf_file.name.replace(".vcf","")

    """
    bcftools query -f'%CHROM,%SVTYPE,%SVLEN\\n' ${vcf_file} | awk -F',' 'BEGIN {OFS = FS} \$1 ~/^[0-9]*\$|^X\$/{abs=\$3;sub("^-", "",abs); print \$1,\$2,abs}' > "${simple_name}.pacbio.csv"
    """
}

process size_data_from_previous {

    input:
        file vcf_file

    output:
        file "*.csv"

    script:
    
    simple_name = vcf_file.name.replace(".vcf","")
    strain = simple_name.tokenize('-').get(0)

    previous_dir = file(params.previous_dir)
    """
    for file in ${previous_dir}/${strain}*; do
        TYPE="\$(echo \$file | cut -d '.' -f 2)"
        cat \$file | awk -v type=\$TYPE '{print \$1","type","\$4}' >> "${simple_name}.ilumina.csv"
    done
    """
}

process calculate_size_distribution {

    publishDir file(params.out_dir), mode: "copy",  saveAs: {filename -> filename.tokenize('-').get(0) + '/' + filename.replace(".sizes.png","")+'/figure_sizes.png' }
    
    input:
        tuple val(simple_name), file(pacbio_data), file(ilumina_data)
        file size_distribution_script
    
    output:
        file "*.png"

    script:

    strain = simple_name.tokenize('-').get(0)

    """
    grep -v BND ${pacbio_data} > pacbio.csv

    python ${size_distribution_script} -s ${strain} -i1 pacbio.csv -i2 ${ilumina_data} -o "${simple_name}.sizes.png"

    rm pacbio.csv
    """
}