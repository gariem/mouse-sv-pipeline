#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

params.validated_src_files = './data/validated/original/{DBA_2J,C57BL_6NJ}*.bed'
params.validated_output = './data/validated/simple'


process prepare_validated_files {

    publishDir file(params.validated_output), mode: "move", overwrite: true
    
    input:
        tuple val(strain), val(type), file(bed_file)
    output:
        file '*bed'


    """
    cat ${bed_file} | awk '!x[\$0]++' > ${strain}.${type}.validated.bed
    """

}

workflow {

    Channel.fromPath(params.validated_src_files).set{validated}

    validated.map { file ->
        def identifier = file.name.replace("mm39.bed", "")
        def strain = identifier.tokenize('.').get(0)
        def type = identifier.contains('H1') || identifier.contains('H2') ? 'DEL' : 'INS'
        return tuple(strain, type, file)
    }
    .groupTuple(by: [0,1], size: 2)
    .set{ validated_grouped }

    prepare_validated_files(validated_grouped).view()

}