#!/usr/bin/env nextflow

params.strain = "*"

params.input_dir = "./data/input"
params.input_files = "${params.input_dir}/${params.strain}*.vcf"
params.output_dir = "./data/merged"

params.mappings = "./data/input/merge_mappings.txt"
params.filter_hets = '0'

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
        val filter_hets from params.filter_hets

    output:
        file "*.bed" into bed_files
    
    script:

    fileName = vcf_file.getName().toString().tokenize('.').get(0)
    caller = fileName.tokenize('-').get(1)
    strain = fileName.tokenize('-').get(0)
    """
    while read -r line
    do
        TYPE="\$(echo \$line | cut -d ':' -f 1)"
        
        if [ "${filter_hets}" == "1" ]; then
            MAPPING="\$(echo "\$line && GT='HOM'" | cut -d ':' -f 2)"
        else
            MAPPING="\$(echo \$line | cut -d ':' -f 2)"
        fi

        bcftools query -i"\$MAPPING" -f'%CHROM\\t%POS0\\t%END0\\t%SVLEN\\n' ${vcf_file} | \
        awk -F'\\t' 'BEGIN {OFS = FS} \$1 ~/^[0-9]*\$|^X\$/{print \$1,\$2,\$3,\$4}' >> "${strain}-${caller}-\$TYPE.bed"
    done < ${mappings_file}
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
    SURVIVOR bedtovcf ${bed_file} ${type} '${strain}-${caller}-${type}.vcf'
    rm ${bed_file}
    """

}


// Transform channel, group by caller (from file name)
vcf_files.map{ file ->
    def name_parts = file.getName().toString().tokenize('.').get(0).tokenize('-')
    def strain_caller = name_parts.get(0) + '-' + name_parts.get(1)
    return tuple(strain_caller, file)
}
.groupTuple()
.set{ grouped_vcfs }


// Join by-type VCF Files into single per-caller VCF File
process join_mapped_vcfs {

    input:
        set strain_caller, file(mapped_vcf) from grouped_vcfs

    output:
        file '*.vcf' into final_vcfs

    script:


    """
    for FILE in ${mapped_vcf}
    do
        if [[ -f ${strain_caller}-joint.vcf ]]
        then
            bcftools view --no-header \$FILE >> '${strain_caller}-joint.vcf'
        else
            bcftools view \$FILE > '${strain_caller}-joint.vcf'
        fi
    done
    rm -rf ${mapped_vcf}
    """
}

final_vcfs.flatMap{ file ->
    def strain = file.getName().toString().tokenize('-').get(0)
    def max_dist_arr = params.max_dist.toString().split(',')
    def min_size_arr = params.min_size.toString().split(',')
    def min_callers_arr = params.min_callers.toString().split(',')

    def tuples = []

    for(max_dist in max_dist_arr){
        for(min_callers in min_callers_arr){
            for(min_size in min_size_arr){
                tuples.add(tuple(strain, max_dist, min_callers, min_size, file))
            }
        }
    }
    return tuples
}
.groupTuple(by: [0, 1, 2, 3])
.set{ grouped_final_vcfs }


// Mege VCF Files using SURVIVOR
process merge_mapped_vcfs {

    publishDir file("${params.output_dir}"), mode: "move", pattern: "*.vcf"

    input:
        set strain, max_dist, min_callers, min_size, file(vcf_to_merge) from grouped_final_vcfs
        val filter_hets from params.filter_hets

    output:
        file "${strain}-survivor*.vcf" into merged_vcf


    """
    echo "${vcf_to_merge}" | tr ' ' '\n' > '${strain}-merged-inputlist.txt'
    SUB=""
    if [ "${filter_hets}" == "1" ]; then
        SUB="_nohets"
    fi
    SURVIVOR merge '${strain}-merged-inputlist.txt' ${max_dist} ${min_callers} ${params.same_type} 1 0 ${params.min_size} "${strain}-survivor_${max_dist}_${min_callers}_${min_size}\$SUB.vcf"
    rm -rf ${vcf_to_merge}
    """ 

}