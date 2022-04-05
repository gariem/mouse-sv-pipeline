#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// Generate VCF Files from output BED in previous proccess (Mapped VCF)
process survivor_vcf_from_bed {

    input:
        tuple val(strain), val(type), file(input_files), file(validation_file)
    
    output:
        tuple val(strain), val(type), file("*.vcf")

    script:

    if(input_files[0].name.contains("minigraph")){
        minigraph_file = input_files[0]
        pbsv_file = input_files[1]
    }else{
        minigraph_file = input_files[1]
        pbsv_file = input_files[0]
    }

    """
    SURVIVOR bedtovcf ${minigraph_file} ${type} '${strain}-minigraph-${type}.vcf'
    SURVIVOR bedtovcf ${pbsv_file} ${type} '${strain}-pbsv-${type}.vcf'
    """
}

// Join by-type VCF Files into single per-caller VCF File
process join_survivor_vcfs {

    input:
        tuple val(strain), val(caller), val(types), file(vcf_files)

    output:
        tuple val(strain), val(caller), file("*.vcf")

    script:

    """
    for FILE in ${vcf_files}
    do
        if [[ -f ${strain}-${caller}-joint.vcf ]]
        then
            bcftools view --no-header \$FILE >> "${strain}-${caller}-joint.vcf"
        else
            bcftools view \$FILE > "${strain}-${caller}-joint.vcf"
        fi
    done
    """
}

// Mege VCF Files using SURVIVOR
process merge_survivor_vcfs {

    //publishDir file("${params.output_dir}"), mode: "move", pattern: "*.vcf"

    input:
        tuple val(strain), val(callers), file(vcf_files)

    output:
        tuple val(strain), file("*.vcf")

    """
    echo "${vcf_files}" | tr ' ' '\n' > '${strain}-merged-inputlist.txt'

    SURVIVOR merge '${strain}-merged-inputlist.txt' 200 1 1 1 0 20 "${strain}-survivor_200_1_20.vcf"
    """ 
}