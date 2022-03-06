#!/usr/bin/env nextflow

nextflow.enable.dsl = 2


process bed_from_full_graph {

    input: 
        val strain
        file folder
    
    output:
        file "*.bed"

    script:
    
    simple_name = strain + '-minigraph'

    """
    awk -F"[\t:]" 'BEGIN {OFS = "\t"} {if(\$6!="."&&(\$3-\$2)<\$7)print \$1,\$2,\$3,\$7}' ${folder}/${strain}.bed > "${strain}-minigraph__INS.bed"
    awk -F"[\t:]" 'BEGIN {OFS = "\t"} {if(\$6!="."&&(\$3-\$2)>\$7)print \$1,\$2,\$3,\$7}' ${folder}/${strain}.bed > "${strain}-minigraph__DEL.bed"
    """

}


process rename_tuples {

    input: 
        tuple val(simple_name), val(type), file(bed_file)
    
    output:
        tuple val(type), file("*.bed")

    script:
    
    """
    cat ${bed_file} > "${simple_name}.${bed_file.name}.bed"
    """
}

process intersect_all_minigraph {

    publishDir file(params.out_dir + '/minigraph_all'), mode: "copy"

    input: 
        tuple val(type), file(bed_files)
        val overlap_fraction

    output:
        file "*.bed"

    script:
    def b6nj = bed_files.findAll{ it.name.contains("C57BL_6NJ") }.join(" ")
    def others = bed_files.findAll{ !it.name.contains("C57BL_6NJ") }.join(" ")

    """
    bedtools intersect -wo -filenames -f ${overlap_fraction} -a ${b6nj} -b ${others} > minigraph-${type}-all.bed
    """

}