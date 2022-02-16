#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process generate_len_csv_from_vcf {

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

process generate_len_csv_from_previous_tsv {

    input:
        file vcf_file
        file previous_dir

    output:
        file "*.tsv"

    script:
    
    simple_name = vcf_file.name.replace(".vcf","")
    strain = simple_name.tokenize('-').get(0)
    """
    for file in ${previous_dir}/${strain}*; do
        TYPE="\$(echo \$file | cut -d '.' -f 2)"
        cat \$file | awk -v type=\$TYPE '{print \$1","type","\$4}' >> "${simple_name}.ilumina.tsv"
    done
    """
}

process calculate_size_distribution {

    input:
        file pacbio_data
        file ilumina_data
        file size_distribution_script
    
    output:
        file "*.png"

    script:

    simple_name = pacbio_data.name.replace(".pacbio.csv","")
    strain = simple_name.tokenize('-').get(0)

    """
    python ${size_distribution_script} -s ${strain} -i1 ${pacbio_data} -i2 ${ilumina_data} -o "${simple_name}.sizes.png"
    """
}