#!/usr/bin/env nextflow

nextflow.enable.dsl = 2


process intersect_files { 

    input: 
        file bed_file_a
        file bed_file_b
        val window
        val options
        val outname

    output:
        file "${outname}"

    script:

    """
    awk -F'\\t' 'BEGIN {OFS = FS} {print \$1,\$2-${window},\$3+${window},\$4}' ${bed_file_a} | awk '!x[\$0]++' > FILE_A
    awk -F'\\t' 'BEGIN {OFS = FS} {print \$1,\$2-${window},\$3+${window},\$4}' ${bed_file_b} | awk '!x[\$0]++' > FILE_B

    bedtools intersect -a FILE_A -b FILE_B ${options} > "${outname}"

    """
}