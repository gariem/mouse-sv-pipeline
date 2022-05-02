#!/usr/bin/env nextflow

nextflow.enable.dsl = 2


process combine_features { 

    input:
        tuple val(strain), val(type), file(input_files), file(validation_file)

    output:
        tuple val(strain), val(type), val(strain), file("*.bed")

    script:

    if(input_files[0].name.contains("minigraph")){
        minigraph_file = input_files[0]
        pbsv_file = input_files[1]
    }else{
        minigraph_file = input_files[1]
        pbsv_file = input_files[0]
    }

    """
    cat ${minigraph_file} > ${strain}-${type}-combined.bed
    bedtools intersect -a ${pbsv_file} -b ${minigraph_file} -f 0.8 -r -v >> ${strain}-${type}-combined.bed
    """

}