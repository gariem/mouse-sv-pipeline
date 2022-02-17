#!/usr/bin/env nextflow

nextflow.enable.dsl = 2


process intersect_files { 

    input: 
        tuple val(identifier), val(simple_name), file(b_file), file(a_file)
        val options
        each window

    output:
        file "*.bed"

    script:

    """
    awk -F'\\t' 'BEGIN {OFS = FS} {print \$1,\$2-${window},\$3+${window},\$4}' ${a_file} > FILE_A
    awk -F'\\t' 'BEGIN {OFS = FS} {print \$1,\$2-${window},\$3+${window},\$4}' ${b_file} > FILE_B

    bedtools intersect -a FILE_A -b FILE_B ${options} > intersect_wa


    INTERSECTED="\$(cat intersect_wa | uniq -u | wc -l)"
    TOTAL="\$(cat FILE_A | uniq -u | wc -l)"

    mv intersect_wa "${simple_name}__${identifier.tokenize('.').get(1)}.\${INTERSECTED}_of_\${TOTAL}.bed"
    """
}

process retrieve_validated_features {

    input:
        tuple val(key), file(bed_file)
        val suffix
    
    output:
        file "*.bed"

    script:

    """
    cat ${bed_file} | awk '!x[\$0]++' > ${key}${suffix}.bed
    """    
}