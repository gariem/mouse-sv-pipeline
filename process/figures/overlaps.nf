#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process split_data {

    input:
        tuple val(strain), val(type), file(input_files)
        val start
        val end

    output:
        tuple val(strain), val(type), val(range), file("*.bed")
    
    script:

    range = start.toString() + '_' + end.toString()

    if(input_files[0].name.contains("minigraph")){
        minigraph_file = input_files[0]
        pbsv_file = input_files[1]
    }else{
        minigraph_file = input_files[1]
        pbsv_file = input_files[0]
    }


    """
    cat ${minigraph_file} | awk -F'\t' 'BEGIN {OFS = FS} {print \$1,\$2,\$3,\$4<0?\$4*-1:\$4}' | awk '\$4 >= ${start} && \$4 <= ${end}' > "${strain}.${type}.${range}.minigraph.bed"
    cat ${pbsv_file} | awk -F'\t' 'BEGIN {OFS = FS} {print \$1,\$2,\$3,\$4<0?\$4*-1:\$4}' | awk '\$4 >= ${start} && \$4 <= ${end}' > "${strain}.${type}.${range}.pbsv.bed"
    """
}

process calc_overlaps {

    input:
        tuple val(strain), val(type), val(range), file(input_files)

    output:
        file("*.csv")

    script: 

    if(input_files[0].name.contains("minigraph")){
        minigraph_file = input_files[0]
        pbsv_file = input_files[1]
    }else{
        minigraph_file = input_files[1]
        pbsv_file = input_files[0]
    }

    """
    cat ${minigraph_file} | wc -l
    cat ${pbsv_file} | wc -l

    PBSV_COUNT="\$(cat ${pbsv_file} | uniq -u | wc -l)"
    MNG_COUNT="\$(cat ${minigraph_file} | uniq -u | wc -l)"

    A_IN_B="\$(bedtools intersect -a ${minigraph_file} -b ${pbsv_file} | uniq -u | wc -l )"
    B_IN_A="\$(bedtools intersect -a ${pbsv_file} -b ${minigraph_file} | uniq -u | wc -l )"

    # echo "STRAIN,TYPE,RANGE,PBSV_COUNT,MNG_COUNT,A_IN_B,B_IN_A" > data.csv
    echo "${strain},${type},${range},\$PBSV_COUNT,\$MNG_COUNT,\$A_IN_B,\$B_IN_A" > "${strain}.${type}.${range}.data.csv"

    """
}

process draw_overlaps {

    input:
        file csv_file
    
    output:
        file "*.data"

    """
    cat ${csv_file} > data.data
    """
}

