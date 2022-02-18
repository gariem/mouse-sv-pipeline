#!/usr/bin/env nextflow

nextflow.enable.dsl = 2


process intersect_features { 

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

process calc_intersect_stats {

    input:
        tuple val(simple_name), file(bed_file)
        val suffix
    output:
        file "*.data"
    
    script:
    data_file = simple_name + '_-_' + suffix
    """
    ls ${bed_file} > tmp
    sed -i "s/${simple_name}__//g" tmp

    echo "ID=${simple_name}" > data

    GLOBAL_TOTAL=0
    GLOBAL_INTERSECT=0

    while read -r stat_line
    do
        TYPE="\$(echo \$stat_line | cut -d '.' -f 1)"
        DATA="\$(echo \$stat_line | cut -d '.' -f 2)"
        INTERSECT="\$(echo \$DATA | cut -d '_' -f 1)"
        TOTAL="\$(echo \$DATA | cut -d '_' -f 3)"
        
        GLOBAL_TOTAL=\$((\$GLOBAL_TOTAL + \$TOTAL))
        GLOBAL_INTERSECT=\$((\$GLOBAL_INTERSECT + \$INTERSECT))

        echo "\${TYPE}_TOT=\${TOTAL}" >> data
        echo "\${TYPE}_INT=\${INTERSECT}"  >> data
    done < tmp

    SCORE="\$(printf '%.2f\n' \$(echo "scale=2; \${GLOBAL_INTERSECT}/\${GLOBAL_TOTAL}" | bc -l))" 
    echo "SCORE=\${SCORE}"  >> data

    mv data ${data_file}_\${SCORE}.data
    """
}

process save_intersect_stats {

    publishDir file(params.out_dir), mode: "copy", overwrite: true,  saveAs: {filename -> filename.split('-')[0] + '/' + filename.split('_-_')[0] + '/' + filename.split('_-_')[1] }

    input:
        file stat_file
    
    output:
        file "*"

    """
    ln -s ${stat_file} ${stat_file.name.replace(".data", "")}
    """

}