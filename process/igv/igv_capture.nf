#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process radomize_bed_file {
    maxForks 4

    input:
        file bed_file

    output:
        file "*.data"

    script:
    
    """
    shuf -n ${params.random_sample} ${bed_file} > ${bed_file}.data
    """
}

process take_screenshots { 
    
    //save as: {out_dir}/strain/simple_name/chr_pos_end_slopX.png
    publishDir file(params.out_dir), mode: "copy",  saveAs: {
                    filename -> filename.split('_-_')[1].tokenize('-').get(0) + '/' + filename.split('_-_')[1].replace(".png", "") + '/' + filename.replace(filename.split('_-_')[1], "").replace("_-_","") + ".png"
                }    

    input:
        tuple val(strain), val(dir), val(simple_name), file(bed_file)
        file igv_workdir

    output:
        file "${dir}/*.png"

    script:
    type = bed_file.name.toUpperCase().contains("INS") ? "INS_" : bed_file.name.toUpperCase().contains("DEL") ? "DEL_" : ""

    template 'igv_screenshot.sh'

}