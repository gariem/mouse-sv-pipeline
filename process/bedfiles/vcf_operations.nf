#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process bed_from_vcf {

    input: 
        file vcf_file
        each type
        val suffix
    
    output:
        file "*.bed"

    script:

    simple_name = vcf_file.name.replace(".vcf","")
    strain = simple_name.tokenize('-').get(0)

    """
    bcftools query -i"SVTYPE='${type}'" -f'%CHROM\\t%POS0\\t%END0\\t%SVLEN\\n' ${vcf_file} | \
            awk -F'\\t' 'BEGIN {OFS = FS} {print \$1,\$2,\$3,\$4}' > "${simple_name}__${type}${suffix}.bed"
    """
}

process filter_prob_dysgu {

    input: 
        file vcf_file
        each prob

    output:
        file '*.vcf'

    script:

    simple_name = vcf_file.name.replace(".vcf","")
    prob_percent = Double.parseDouble(prob)/100

    """
    bcftools view -i'PROB>=${prob_percent}' ${vcf_file} > "${simple_name}.p${prob}.vcf"
    """
}

process clean_hets {

    input: 
        file vcf_file

    output:
        file '*.vcf'

    script:
    
    simple_name = vcf_file.name.replace(".all.vcf","")

    """
    bcftools view -i'GT="HOM"' ${vcf_file} > "${simple_name}.hom.vcf"
    """
}



process clean_regions {

    input: 
        file vcf_file

    output:
        file '*.vcf'

    script:
    simple_name = vcf_file.name.replace(".vcf","")

    """
    output=\$(bcftools view -i'CHROM="1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,X"' ${vcf_file} > ${simple_name}.all.vcf || true)
    
    """
}
