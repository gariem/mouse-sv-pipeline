#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process split_data {

    input:
        tuple val(strain), val(type), file(input_files), file(validation_file)
        val start
        val end

    output:
        tuple val(strain), val(type), val(range), file("*.pbsv.bed"), file("*.minigraph.bed"), file("*.validation.bed")
    
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
    cat ${pbsv_file} | awk -F'\t' 'BEGIN {OFS = FS} {print \$1,\$2,\$3,\$4<0?\$4*-1:\$4}' | awk '\$4 >= ${start} && \$4 <= ${end}' > "${strain}.${type}.${range}.pbsv.bed"
    cat ${minigraph_file} | awk -F'\t' 'BEGIN {OFS = FS} {print \$1,\$2,\$3,\$4<0?\$4*-1:\$4}' | awk '\$4 >= ${start} && \$4 <= ${end}' > "${strain}.${type}.${range}.minigraph.bed"
    cat ${validation_file} | awk -F'\t' 'BEGIN {OFS = FS} {print \$1,\$2,\$3,\$4<0?\$4*-1:\$4}' | awk '\$4 >= ${start} && \$4 <= ${end}' > "${strain}.${type}.${range}.validation.bed"
    """
}

process calc_overlaps {

    input:
        tuple val(strain), val(type), val(range), file(pbsv_file), file(minigraph_file), file(validation_file)
        val window

    output:
        file("*.csv")

    script: 

    """
    awk -F'\\t' 'BEGIN {OFS = FS} {print \$1,\$2-${window},\$3+${window},\$4}' ${pbsv_file} > FILE_A
    awk -F'\\t' 'BEGIN {OFS = FS} {print \$1,\$2-${window},\$3+${window},\$4}' ${minigraph_file} > FILE_B

    awk -F'\\t' 'BEGIN {OFS = FS} {print \$1,\$2-30,\$3+30,\$4}' ${pbsv_file} > PBSV_VALIDATION
    awk -F'\\t' 'BEGIN {OFS = FS} {print \$1,\$2-30,\$3+30,\$4}' ${minigraph_file} > MG_VALIDATION
    awk -F'\\t' 'BEGIN {OFS = FS} {print \$1,\$2-30,\$3+30,\$4}' ${validation_file} > VALIDATION

    PBSV_COUNT="\$(cat FILE_A | uniq -u | wc -l)"
    MNG_COUNT="\$(cat FILE_B | uniq -u | wc -l)"
    
    VAL_COUNT="\$(cat VALIDATION | uniq -u | wc -l)"

    bedtools intersect -a FILE_A -b FILE_B -wa > INTER

    A_IN_B="\$(cat INTER | uniq -u | wc -l )"

    VAL_IN_PBSV="\$(bedtools intersect -a VALIDATION -b PBSV_VALIDATION -wa | uniq -u | wc -l )"
    VAL_IN_MINI="\$(bedtools intersect -a VALIDATION -b MG_VALIDATION -wa | uniq -u | wc -l )"
    VAL_IN_INTER="\$(bedtools intersect -a VALIDATION -b INTER -wa | uniq -u | wc -l )"

    echo "${window},${strain},${type},${range},\$PBSV_COUNT,\$MNG_COUNT,\$A_IN_B,\$VAL_COUNT,\$VAL_IN_PBSV,\$VAL_IN_MINI,\$VAL_IN_INTER" > "${strain}.${type}.${range}_w${window}.data.csv"

    """
}

process draw_overlaps {

    input:
        file csv_file
    
    output:
        file "*.csv"

    """
    echo "window,strain,type,range,pbsv,minigraph,intersect,t_val,pbsv_val,mini_val,inter_val" >> pbsv_minigraph_overlap.csv
    cat ${csv_file} >> pbsv_minigraph_overlap.csv
     ## 
    """
}

