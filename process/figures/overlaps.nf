#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process split_triple_data {

    input:
        tuple val(strain), val(type), file(input_files), file(validation_file)

    output:
        tuple val(strain), val(type), file("*.pbsv.bed"), file("*.minigraph.bed"), file("*.validation.bed")
    
    script:

    if(input_files[0].name.contains("minigraph")){
        minigraph_file = input_files[0]
        pbsv_file = input_files[1]
    }else{
        minigraph_file = input_files[1]
        pbsv_file = input_files[0]
    }


    """
    start=0
    end=50
    range="\${start}_\${end}"
    cat ${pbsv_file} | awk -F'\t' 'BEGIN {OFS = FS} {print \$1,\$2,\$3,\$4<0?\$4*-1:\$4}' | awk '\$4 >= 0 && \$4 <= 50' > "${strain}.${type}.\${range}.pbsv.bed"
    cat ${minigraph_file} | awk -F'\t' 'BEGIN {OFS = FS} {print \$1,\$2,\$3,\$4<0?\$4*-1:\$4}' | awk '\$4 >= 0 && \$4 <= 50' > "${strain}.${type}.\${range}.minigraph.bed"
    cat ${validation_file} | awk -F'\t' 'BEGIN {OFS = FS} {print \$1,\$2,\$3,\$4<0?\$4*-1:\$4}' | awk '\$4 >= 0 && \$4 <= 50' > "${strain}.${type}.\${range}.validation.bed"

    start=50
    end=100
    range="\${start}_\${end}"
    cat ${pbsv_file} | awk -F'\t' 'BEGIN {OFS = FS} {print \$1,\$2,\$3,\$4<0?\$4*-1:\$4}' | awk '\$4 >= 50 && \$4 <= 100' > "${strain}.${type}.\${range}.pbsv.bed"
    cat ${minigraph_file} | awk -F'\t' 'BEGIN {OFS = FS} {print \$1,\$2,\$3,\$4<0?\$4*-1:\$4}' | awk '\$4 >= 50 && \$4 <= 100' > "${strain}.${type}.\${range}.minigraph.bed"
    cat ${validation_file} | awk -F'\t' 'BEGIN {OFS = FS} {print \$1,\$2,\$3,\$4<0?\$4*-1:\$4}' | awk '\$4 >= 50 && \$4 <= 100' > "${strain}.${type}.\${range}.validation.bed"
    

    start=100
    end=1000
    range="\${start}_\${end}"
    cat ${pbsv_file} | awk -F'\t' 'BEGIN {OFS = FS} {print \$1,\$2,\$3,\$4<0?\$4*-1:\$4}' | awk '\$4 >= 100 && \$4 <= 1000' > "${strain}.${type}.\${range}.pbsv.bed"
    cat ${minigraph_file} | awk -F'\t' 'BEGIN {OFS = FS} {print \$1,\$2,\$3,\$4<0?\$4*-1:\$4}' | awk '\$4 >= 100 && \$4 <= 1000' > "${strain}.${type}.\${range}.minigraph.bed"
    cat ${validation_file} | awk -F'\t' 'BEGIN {OFS = FS} {print \$1,\$2,\$3,\$4<0?\$4*-1:\$4}' | awk '\$4 >= 100 && \$4 <= 1000' > "${strain}.${type}.\${range}.validation.bed"

    start=1000
    end=10000
    range="\${start}_\${end}"
    cat ${pbsv_file} | awk -F'\t' 'BEGIN {OFS = FS} {print \$1,\$2,\$3,\$4<0?\$4*-1:\$4}' | awk '\$4 >= 1000 && \$4 <= 10000' > "${strain}.${type}.\${range}.pbsv.bed"
    cat ${minigraph_file} | awk -F'\t' 'BEGIN {OFS = FS} {print \$1,\$2,\$3,\$4<0?\$4*-1:\$4}' | awk '\$4 >= 1000 && \$4 <= 10000' > "${strain}.${type}.\${range}.minigraph.bed"
    cat ${validation_file} | awk -F'\t' 'BEGIN {OFS = FS} {print \$1,\$2,\$3,\$4<0?\$4*-1:\$4}' | awk '\$4 >= 1000 && \$4 <= 10000' > "${strain}.${type}.\${range}.validation.bed"

    start=10000
    end=10000000
    range="\${start}_\${end}"
    cat ${pbsv_file} | awk -F'\t' 'BEGIN {OFS = FS} {print \$1,\$2,\$3,\$4<0?\$4*-1:\$4}' | awk '\$4 >= 10000 && \$4 <= 1000000' > "${strain}.${type}.\${range}.pbsv.bed"
    cat ${minigraph_file} | awk -F'\t' 'BEGIN {OFS = FS} {print \$1,\$2,\$3,\$4<0?\$4*-1:\$4}' | awk '\$4 >= 10000 && \$4 <= 1000000' > "${strain}.${type}.\${range}.minigraph.bed"
    cat ${validation_file} | awk -F'\t' 'BEGIN {OFS = FS} {print \$1,\$2,\$3,\$4<0?\$4*-1:\$4}' | awk '\$4 >= 10000 && \$4 <= 1000000' > "${strain}.${type}.\${range}.validation.bed"
    
    """
}

process split_survivor_data {
    input:
        tuple val(strain), val(type), file(validation_file), file(survivor_file)

    output:
        tuple val(strain), val(type), file("*.survivor.bed"), file("*.validation.bed")
    
    script:

    """
    start=0
    end=50
    range="\${start}_\${end}"
    cat ${survivor_file} | awk -F'\t' 'BEGIN {OFS = FS} {print \$1,\$2,\$3,\$4<0?\$4*-1:\$4}' | awk '\$4 >= 0 && \$4 <= 50' > "${strain}.${type}.\${range}.survivor.bed"
    cat ${validation_file} | awk -F'\t' 'BEGIN {OFS = FS} {print \$1,\$2,\$3,\$4<0?\$4*-1:\$4}' | awk '\$4 >= 0 && \$4 <= 50' > "${strain}.${type}.\${range}.validation.bed"

    start=50
    end=100
    range="\${start}_\${end}"
    cat ${survivor_file} | awk -F'\t' 'BEGIN {OFS = FS} {print \$1,\$2,\$3,\$4<0?\$4*-1:\$4}' | awk '\$4 >= 50 && \$4 <= 100' > "${strain}.${type}.\${range}.survivor.bed"
    cat ${validation_file} | awk -F'\t' 'BEGIN {OFS = FS} {print \$1,\$2,\$3,\$4<0?\$4*-1:\$4}' | awk '\$4 >= 50 && \$4 <= 100' > "${strain}.${type}.\${range}.validation.bed"
    
    start=100
    end=1000
    range="\${start}_\${end}"
    cat ${survivor_file} | awk -F'\t' 'BEGIN {OFS = FS} {print \$1,\$2,\$3,\$4<0?\$4*-1:\$4}' | awk '\$4 >= 100 && \$4 <= 1000' > "${strain}.${type}.\${range}.survivor.bed"
    cat ${validation_file} | awk -F'\t' 'BEGIN {OFS = FS} {print \$1,\$2,\$3,\$4<0?\$4*-1:\$4}' | awk '\$4 >= 100 && \$4 <= 1000' > "${strain}.${type}.\${range}.validation.bed"

    start=1000
    end=10000
    range="\${start}_\${end}"
    cat ${survivor_file} | awk -F'\t' 'BEGIN {OFS = FS} {print \$1,\$2,\$3,\$4<0?\$4*-1:\$4}' | awk '\$4 >= 1000 && \$4 <= 10000' > "${strain}.${type}.\${range}.survivor.bed"
    cat ${validation_file} | awk -F'\t' 'BEGIN {OFS = FS} {print \$1,\$2,\$3,\$4<0?\$4*-1:\$4}' | awk '\$4 >= 1000 && \$4 <= 10000' > "${strain}.${type}.\${range}.validation.bed"

    start=10000
    end=10000000
    range="\${start}_\${end}"
    cat ${survivor_file} | awk -F'\t' 'BEGIN {OFS = FS} {print \$1,\$2,\$3,\$4<0?\$4*-1:\$4}' | awk '\$4 >= 10000 && \$4 <= 1000000' > "${strain}.${type}.\${range}.survivor.bed"
    cat ${validation_file} | awk -F'\t' 'BEGIN {OFS = FS} {print \$1,\$2,\$3,\$4<0?\$4*-1:\$4}' | awk '\$4 >= 10000 && \$4 <= 1000000' > "${strain}.${type}.\${range}.validation.bed"
    
    """
}

process calc_overlaps {

    input:
        tuple val(strain), val(type), val(range), file(input_files)
        val window

    output:
        file("*.csv")

    script: 

    for(int i = 0;i<=2;i++) {
        name = input_files.get(i).getName()
        if(name.contains('.pbsv.')){
            pbsv_file = input_files.get(i)
            continue
        }
        if(name.contains('.minigraph.')){
            minigraph_file = input_files.get(i)
            continue
        }
        if(name.contains('.validation.')){
            validation_file = input_files.get(i)
            continue
        }
    }

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

process calc_survivor_scores {

    input:
        tuple val(strain), val(type), val(range), file(input_files)
        val window

    output:
        tuple file("*.csv"), file("*.bed")

    script: 

    for(int i = 0;i<=1;i++) {
        name = input_files.get(i).getName()
        if(name.contains('.survivor.')){
            survivor_file = input_files.get(i)
            continue
        }
        if(name.contains('.validation.')){
            validation_file = input_files.get(i)
            continue
        }
    }

    """
    awk -F'\\t' 'BEGIN {OFS = FS} {print \$1,\$2-30,\$3+30,\$4}' ${survivor_file} > SURVIIVOR
    awk -F'\\t' 'BEGIN {OFS = FS} {print \$1,\$2-30,\$3+30,\$4}' ${validation_file} > VALIDATION

    SURV_COUNT="\$(cat SURVIIVOR | uniq -u | wc -l)"    
    VAL_COUNT="\$(cat VALIDATION | uniq -u | wc -l)"

    VAL_IN_SUVR="\$(bedtools intersect -a VALIDATION -b SURVIIVOR -wa | uniq -u | wc -l )"

    bedtools intersect -a VALIDATION -b SURVIIVOR -v > "${strain}-${range}-${type}-missed.bed"

    echo "${window},${strain},${type},${range},\$SURV_COUNT,\$VAL_COUNT,\$VAL_IN_SUVR" > "${strain}.${type}.${range}_w${window}.data.csv"

    """
}

process summarize_survivor {
    input:
        file csv_file
    
    output:
        file "*.csv"

    script:    
    """
    echo "window,strain,type,range,survivor,validated,survivor" >> pbsv_minigraph_merged_score.csv
    cat ${csv_file} >> pbsv_minigraph_merged_score.csv
     ## 
    """

}

process summarize_overlaps {

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

