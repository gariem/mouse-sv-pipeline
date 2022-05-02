#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process bed_from_full_graph {

    input: 
        file graph_bed
    
    output:
        file "*.bed"

    script:
    
    strain = graph_bed.name.replace(".bed", "")

    """
    awk -F"[\t:]" 'BEGIN {OFS = "\t"} {if(\$6!="."&&(\$3-\$2)<\$7)print \$1,\$2,\$3,\$7-(\$3-\$2)}' ${graph_bed} > "${strain}-minigraph__INS.bed"
    awk -F"[\t:]" 'BEGIN {OFS = "\t"} {if(\$6!="."&&(\$3-\$2)>\$7)print \$1,\$2,\$3,\$7-(\$3-\$2)}' ${graph_bed} > "${strain}-minigraph__DEL.bed"
    awk -F"[\t:]" 'BEGIN {OFS = "\t"} {if(\$6!="."&&(\$3-\$2)==\$7)print \$1,\$2,\$3,\$7-(\$3-\$2)}' ${graph_bed} > "${strain}-minigraph__INV.bed"
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

    output:
        file "*.bed"

    script:

    overlap_fraction = 1.0

    def b6nj = bed_files.findAll{ it.name.contains("C57BL_6NJ") }.join(" ")
    def others = bed_files.findAll{ !it.name.contains("C57BL_6NJ") }.join(" ")

    """
    bedtools intersect -wo -filenames -f ${overlap_fraction} -a ${b6nj} -b ${others} > minigraph-${type}-all.bed
    """

}

process get_features_across_strains {

    input:
        file bed_file

    output:
       tuple file("*.bed"), val(type)

    script:

    type = bed_file.name.contains("INS") ? "INS" : "DEL"
    strains_m1 = params.pattern.count(",")
    
"""
#!python

import pandas as pd

cols = ['Chr1', 'Pos1', 'End1', 'Size1', 'Src', 'Chr2', 'Pos2', 'End2', 'Size2', 'Overlap']
group_cols = ['Chr1', 'Pos1', 'End1', 'Size1', 'Chr2', 'Pos2', 'End2', 'Size2', 'Overlap']
data = pd.read_csv('${bed_file}', sep='\t', names=cols, header=None)

grouped_data = data.groupby(group_cols, as_index=False).count()
repeated = grouped_data[grouped_data["Src"]>=${strains_m1}]
repeated = repeated.drop(['Chr2', 'Pos2', 'End2', 'Size2', 'Overlap', 'Src'], axis=1)

repeated.to_csv("repeated-${type}.bed", sep="\t", header=False, index=False)
"""

}
