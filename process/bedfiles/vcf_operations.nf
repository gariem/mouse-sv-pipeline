#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process bed_from_vcf {

    input: 
        file vcf_file
        each type
    
    output:
        tuple val(strain), val(type), val(simple_name), file('*.bed')

    script:

    simple_name = vcf_file.name.replace(".vcf","")
    strain = simple_name.tokenize('-').get(0)

    """
    bcftools query -i"SVTYPE='${type}'" -f'%CHROM\\t%POS0\\t%END0\\t%SVLEN\\n' ${vcf_file} | \
            awk -F'\\t' 'BEGIN {OFS = FS} {print \$1,\$2,\$3,\$4}' > "${simple_name}__${type}.bed"
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

process filter_read_depth {
    input: 
        file vcf_file
        each read_depth
    
    output:
        file '*.vcf'

    script:

    simple_name = vcf_file.name.replace(".vcf","")

    """
    bcftools view -i'DP>=${read_depth}' ${vcf_file} > "${simple_name}.dp${read_depth}.vcf"
    """
}


process filter_diff_allele_depth {
        
    input: 
        file vcf_file
        each limits
    
    output:
        file '*.vcf'

    script:

    simple_name = vcf_file.name.replace(".vcf","")
    outname = "${simple_name}.ad_${limits[0]}to${limits[1]}.vcf".replace("-","_")

    """
    bcftools view -i'(AD[0:1] - AD[0:0])>=${limits[0]} && (AD[0:1] - AD[0:0])<=${limits[1]}' ${vcf_file} > "${outname}.vcf"
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

process big_inversions {

    publishDir file(params.out_dir), mode: "copy",  saveAs: {filename -> filename.tokenize('-').get(0) + '/' + filename.replace(".bed","")+'/big_invs.bed' }

    input:
        file vcf_file
    
    output:
        file "*.bed" optional true

    script:

    simple_name = vcf_file.name.replace(".vcf","")

    """
    bcftools query -i"SVTYPE='INV' && SVLEN>1000000" -f'%CHROM\\t%POS0\\t%END0\\t%SVLEN\\n' ${vcf_file} | \
            awk -F'\\t' 'BEGIN {OFS = FS} {print \$1,\$2,\$3,\$4}' > big_invs

    if [ -s big_invs ]; then
        mv big_invs ${simple_name}.bed
    fi
    """
}
