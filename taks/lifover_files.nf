#!/usr/bin/env nextflow

params.input_files = "./data/input/validated_mm10.tar.gz"
params.remap_api = '/home/egarcia/appdir/bin/remap_api.pl'

params.remap_src = 'GCF_000001635.20'
params.remap_dest = 'GCF_000001635.27'

params.alias_src = 'mm10'
params.alias_dest = 'mm39'

params.output_dir = './data/validated'


process untar {

    input:
        file 'input_file.tar.gz' from file(params.input_files)

    output:
        file '**.bed' into bed_files

    """
    tar -zxvf input_file.tar.gz
    """

}

process uplift_files {

    maxForks 2
    publishDir file(params.output_dir), mode: "copy"

    input:
        file input_bed from bed_files.flatten()
    
    output:
        file '*.bed'

    script:
    new_name = input_bed.getName().replace(params.alias_src, params.alias_dest)

    """
    perl ${params.remap_api} --mode asm-asm --from ${params.remap_src} --dest ${params.remap_dest} --annotation ${input_bed}  --annot_out ${new_name}
    """
}
