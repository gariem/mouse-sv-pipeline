#!/usr/bin/env nextflow

params.strain = "DBA_2J"

params.input_dir = "./data/input"
params.input_files = "${params.input_dir}/${params.strain}-*.vcf"
params.output_dir = "./data/merged"

params.mappings = "./data/input/mappings.txt"

params.bcftools = '/home/egarcia/appdir/bcftools/bin/bcftools'
params.survivor = '/home/egarcia/workspace/github/SURVIVOR/Debug/SURVIVOR'

params.max_dist = '100'
params.min_callers = '1'
params.same_type = '1'
params.min_size = '30'

Channel.fromPath(params.input_files).set{input_files}

// Generate BED files by type according to the values in data/input/mappings.txt 
process mapped_bed_from_vcf {

    input:
        file vcf_file from input_files
        file mappings_file from file(params.mappings)

    output:
        file "${params.strain}-*.bed" into bed_files
    
    script:

    caller = vcf_file.getName().toString().tokenize('-').get(1).tokenize('.').get(0)
    
    """
    while read -r line
    do
        TYPE="\$(echo \$line | cut -d ':' -f 1)"
        MAPPING="\$(echo \$line | cut -d ':' -f 2)"
        
        ${params.bcftools} query -i"\$MAPPING" -f'%CHROM\\t%POS0\\t%END0\\t%SVLEN\\n' ${vcf_file} | \
        awk -F'\\t' 'BEGIN {OFS = FS} \$1 ~/^[0-9]*\$|^X\$/{print \$1,\$2,\$3,\$4}' >> "${params.strain}-${caller}-\$TYPE.bed"
    done < $mappings_file
    """
}


// Generate VCF Files from output BED in previous proccess (Mapped VCF)
process mapped_vcf_from_bed {

    input:
        file bed_file from bed_files.flatten()
    
    output:
        file "${params.strain}-*.vcf" into vcf_files

    script:

    name_parts = bed_file.getName().toString().tokenize('.').get(0).tokenize('-')
    strain = name_parts.get(0)
    caller = name_parts.get(1)
    type = name_parts.get(2)

    """
    ${params.survivor} bedtovcf ${bed_file} ${type} '${strain}-${caller}-${type}.vcf'
    """

}


// Transform channel, group by caller (from file name)
vcf_files.map{ file ->
    def caller = file.name.toString().tokenize('-').get(1)
    return tuple(caller, file)
}
.groupTuple()
.set{ grouped_vcfs }


// Join by-type VCF Files into single per-caller VCF File
process join_mapped_vcfs {

    input:
        set caller, file(mapped_vcf) from grouped_vcfs

    output:
        file '*.vcf' into final_vcfs

    """
    for FILE in ${mapped_vcf}
    do
        if [[ -f ${params.strain}-${caller}-joint.vcf ]]
        then
            ${params.bcftools} view --no-header \$FILE >> '${params.strain}-${caller}-joint.vcf'
        else
            ${params.bcftools} view \$FILE > '${params.strain}-${caller}-joint.vcf'
        fi
        
    done
    """
}


process merge_mapped_vcfs {

    publishDir file("${params.output_dir}"), mode: "copy", pattern: "*.vcf"

    input:
        file (vcf_to_merge) from final_vcfs.collect()

    output:
        file "${params.strain}-merged_*.vcf" into merged_vcf


    """
        echo "${vcf_to_merge}" | tr " " "\n" > '${params.strain}-merged-inputlist.txt'

        ${params.survivor} merge '${params.strain}-merged-inputlist.txt' ${params.max_dist} ${params.min_callers} ${params.same_type} 1 0 ${params.min_size} '${params.strain}-merged_${params.max_dist}_${params.min_callers}_${params.min_size}.vcf'
    """ 

}