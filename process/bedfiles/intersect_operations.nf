#!/usr/bin/env nextflow

nextflow.enable.dsl = 2


process intersect_features { 

    input: 
        tuple val(source), val(type), file(a_file), file(b_file)
        val options
        each window

    output:
        tuple val(source), val(type), file('*.bed')

    script:

    """
    awk -F'\\t' 'BEGIN {OFS = FS} {print \$1,\$2-${window},\$3+${window},\$4}' ${a_file} > FILE_A
    awk -F'\\t' 'BEGIN {OFS = FS} {print \$1,\$2-${window},\$3+${window},\$4}' ${b_file} > FILE_B

    bedtools intersect -a FILE_A -b FILE_B ${options} > intersect_wa

    INTERSECTED="\$(cat intersect_wa | uniq -u | wc -l)"
    TOTAL_A="\$(cat FILE_A | uniq -u | wc -l)"
    TOTAL_B="\$(cat FILE_B | uniq -u | wc -l)"

    mv intersect_wa "${type}.\${INTERSECTED}_of_\${TOTAL_A}_of_\${TOTAL_B}.bed"
    """
}

process retrieve_validated_features {

    input:
        file(bed_file)
    
    output:
        tuple val(strain), val(type), file('*.bed')

    script:

    strain = bed_file.name.tokenize('.').get(0)
    type = bed_file.name.tokenize('.').get(1)

    """
    cat ${bed_file} | awk '!x[\$0]++' > ${strain}.${type}.validated.bed
    """    
}

process calculate_scores {
    input:
        tuple val(source), val(type), file(previous), file(validated)

    output:
        tuple val(source), file('*_*_*')

    script:

    previous_numbers = previous.name.tokenize('.').get(1).split('_of_')
    previous_intersected = previous_numbers[0]
    previous_total = previous_numbers[1]
    new_total = previous_numbers[2]

    validated_numbers = validated.name.tokenize('.').get(1).split('_of_')
    validated_intersected = validated_numbers[0]
    validated_total = validated_numbers[1]

    """
    echo "${type}.VT=${validated_total}" >> data
    echo "${type}.VI=${validated_intersected}"  >> data
    echo "${type}.PT=${previous_total}" >> data
    echo "${type}.PI=${previous_intersected}"  >> data
    echo "${type}.TOT=${new_total}"  >> data

    SCORE="\$(printf '%.2f\n' \$(echo "scale=2; ${validated_intersected}/${validated_total}" | bc -l))" 
    DIFF="\$(printf '%.2f\n' \$(echo "scale=2; ${new_total}/${previous_total}" | bc -l))" 

    echo "${type}.SCORE=\${SCORE}" >> data
    echo "${type}.DIFF=\${DIFF}" >> data

    mv data "${type}_\${SCORE}_\${DIFF}"
    """
}

process calc_intersect_stats {

    input:
        tuple val(simple_name), file(bed_file)
        val suffix
    output:
        tuple val(simple_name), file('*_*') 
    
    script:
    
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

    mv data ${suffix}_\${SCORE}
    """
}

process save_intersect_stats {

    publishDir file(params.out_dir), mode: "copy", overwrite: true,  saveAs: {filename -> filename.split('-')[1] + '/' + filename.split('-')[2] + '/' + filename.split('-')[0] + '.stats' }

    input:
        tuple val(source_name), file(stat_file)
    
    output:
        file "*"

    """
    ln -s ${stat_file} ${stat_file.name.split('_')[0]}-${source_name}
    """

}

process retrieve_preivous_features {

    input:
        file(bed_file)
    
    output: 
        tuple val(strain), val(type), file('*.bed')
    
    script:

    parts = bed_file.name.tokenize('.')
    strain = parts.get(0)
    type = parts.get(1)

    """
    cat ${bed_file} | awk '!x[\$0]++' > ${strain}.${type}.bed
    """


}