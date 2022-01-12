#!/usr/bin/env nextflow

params.input_files = './data/merged/DBA_2J-merged_100_1_20.vcf'
params.out_dir = './data/reports/raw'
params.previous_dir = './data/previous'
params.validated_dir = './data/validated'

params.merge_mappings = "./data/input/merge_mappings.txt"
params.metrics_mappings = "./data/input/metrics_mappings.txt"

params.filter_hets = '0'

params.intersect_window_a = '0'
params.intersect_window_b = '20'

Channel.fromPath(params.input_files).set{input_files}

process general_metrics {

    input:
        file vcf_file from input_files
        file mappings_file from file(params.merge_mappings)
        val filter_hets from params.filter_hets

    output:
        file '*.bed' into unflattened_bed_files

    script:

    vcf_name = vcf_file.getName().tokenize(".").get(0)
    strain = vcf_name.tokenize("-").get(0) 
    data_file = file("${params.out_dir}/${vcf_name}.data")

    filter_flag = (filter_hets == 1 || vcf_name.contains("_nohets")) ? 1 : 0

    """
    mkdir -p ${file(params.out_dir)}

    TOTAL_CALLS="\$(bcftools view --no-header ${vcf_file} | awk -F'\\t' 'BEGIN {OFS = FS} \$1 ~/^[0-9]*\$|^X\$/{print}' | wc -l)"
    SUP1="\$(bcftools query -i"SUPP='1'" -f'%CHROM\\t%POS\\t%END\\n' ${vcf_file} | awk -F'\\t' 'BEGIN {OFS = FS} \$1 ~/^[0-9]*\$|^X\$/{print}' | wc -l)"
    SUP2="\$(bcftools query -i"SUPP='2'" -f'%CHROM\\t%POS\\t%END\\n' ${vcf_file} | awk -F'\\t' 'BEGIN {OFS = FS} \$1 ~/^[0-9]*\$|^X\$/{print}' | wc -l)"
    SUP3="\$(bcftools query -i"SUPP='3'" -f'%CHROM\\t%POS\\t%END\\n' ${vcf_file} | awk -F'\\t' 'BEGIN {OFS = FS} \$1 ~/^[0-9]*\$|^X\$/{print}' | wc -l)"

    HOM_COUNT="\$(bcftools query -i"GT='HOM'" -f'%CHROM\\t%POS\\t%END\\n' ${vcf_file} | awk -F'\\t' 'BEGIN {OFS = FS} \$1 ~/^[0-9]*\$|^X\$/{print}' | wc -l)"
    HET_COUNT="\$(bcftools query -i"GT='HET'" -f'%CHROM\\t%POS\\t%END\\n' ${vcf_file} | awk -F'\\t' 'BEGIN {OFS = FS} \$1 ~/^[0-9]*\$|^X\$/{print}' | wc -l)"

    echo "STRAIN=${strain}" > ${data_file}
    echo "SRC_FILE=${vcf_file.getName()}" >> ${data_file}
    echo "FILTER=${filter_flag}" >> ${data_file}
    echo "TOTAL=\$TOTAL_CALLS" >> ${data_file}

    echo "S1=\$SUP1" >> ${data_file}
    echo "S2=\$SUP2" >> ${data_file}
    echo "S3=\$SUP3" >> ${data_file}

    echo "HOM=\$HOM_COUNT" >> ${data_file}
    echo "HET=\$HET_COUNT" >> ${data_file}

    echo "W_A=${params.intersect_window_a}" >> ${data_file}
    echo "W_B=${params.intersect_window_a}" >> ${data_file}

    for ((i=1;i<=19;i++)); 
    do 
        if [ "${filter_hets}" == "1" ]; then
            QUERY="CHROM='\$i' && GT='HOM'"
        else
            QUERY="CHROM='\$i'"
        fi

        CHROM_COUNT="\$(bcftools query -i"\$QUERY" -f'%CHROM\\t%POS\\t%END\\n' ${vcf_file} | wc -l)"
        echo "chr\$i=\$CHROM_COUNT" >> ${data_file}
    done

    if [ "${filter_hets}" == "1" ]; then
        CHROM_COUNT="\$(bcftools query -i"CHROM='X' && GT='HOM'" -f'%CHROM\\t%POS\\t%END\\n' ${vcf_file} | wc -l)"
    else
        CHROM_COUNT="\$(bcftools query -i"CHROM='X'" -f'%CHROM\\t%POS\\t%END\\n' ${vcf_file} | wc -l)"
    fi
    
    echo "chrX=\$CHROM_COUNT" >> ${data_file}

    while read -r line
    do
        TYPE="\$(echo \$line | cut -d ':' -f 1)"

        if [ "${filter_hets}" == "1" ]; then
            QUERY="SVTYPE='\$TYPE' && GT='HOM'"
        else
            QUERY="SVTYPE='\$TYPE' "
        fi

        bcftools query -i"\$QUERY" -f'%CHROM\\t%POS0\\t%END0\\t%SVLEN\\n' ${vcf_file} | \
            awk -F'\\t' 'BEGIN {OFS = FS} \$1 ~/^[0-9]*\$|^X\$/{print \$1,\$2,\$3,\$4}' >> "${vcf_name}.\$TYPE.bed"
    done < ${mappings_file}
    """
}

process validated_metrics {

    input:
        file bed_file from unflattened_bed_files.flatten()
        file mappings_file from file(params.metrics_mappings)

    output:
        file '*.bed' optional true into compare_beds

    script:
    
    validated_dir = file(params.validated_dir)

    name_arr = bed_file.getName().tokenize(".")
    name = name_arr.get(0)
    type = name_arr.get(1)

    strain = name_arr.get(0).tokenize("-").get(0)

    previous = file("${params.previous_dir}/${strain}.${type}.mm39.bed")
    data_file = file("${params.out_dir}/${name_arr.get(0)}.data")

    """
    TYPE_COUNT="\$(cat ${bed_file} | wc -l)"
    PREV_TYPE_COUNT="\$(cat ${previous} | wc -l)"

    echo "${type}=\$TYPE_COUNT" >> ${data_file}
    echo "${type}0=\$PREV_TYPE_COUNT" >> ${data_file}

    echo "\$(cat ${mappings_file} | grep ${type} || true;)" > compare_to.txt

    if [[ -s compare_to.txt ]]
    then
        while read -r line
        do
            VALIDATED_TYPE="\$(echo \$line | cut -d '=' -f 1)"
            VALIDATED_FILE=${validated_dir}/${strain}.\$VALIDATED_TYPE.mm39.bed

            if [[ -e \$VALIDATED_FILE ]]
            then
                ln -s \$VALIDATED_FILE "${type}_\$VALIDATED_TYPE.A.${strain}.\$VALIDATED_TYPE.bed"
                ln -s ${bed_file} "${type}_\$VALIDATED_TYPE.B.${name}.${type}.bed"
            fi
        done < compare_to.txt
    fi
    """
}

compare_beds
        .flatten()
        .map { file ->
            def key = file.name.toString().tokenize('.').get(0)
            return tuple(key, file)
        }
        .groupTuple()
        .set{ grouped_beds }


process intersect_files {
    
    input:
        set key, file(bed_files) from grouped_beds
        val window_a from params.intersect_window_a
        val window_b from params.intersect_window_b
    output:
        file '*.dummy' into final_dummies

    script:

    if (bed_files[0].getName().contains('.A.')){
        file_a = bed_files[0]
        file_b = bed_files[1]
    } else {
        file_a = bed_files[1]
        file_b = bed_files[0]
    }
    
    process = file_a.getName().tokenize(".").get(0)
    source_file = file_b.getName().tokenize(".").get(2)
    data_file=file("${params.out_dir}/${source_file}.data")

    """
    A_TOTAL="\$(cat ${file_a} | wc -l)"
    A_NAME="\$(echo ${file_a} | cut -d '.' -f 4)"

    awk -F'\\t' 'BEGIN {OFS = FS} {print \$1,\$2-${window_a},\$3+${window_a},\$4}' ${file_a} > FILE_A
    awk -F'\\t' 'BEGIN {OFS = FS} {print \$1,\$2-${window_b},\$3+${window_b},\$4}' ${file_b} > FILE_B

    bedtools intersect -a FILE_A -b FILE_B -wa > intersect_wa

    A_INTERSECTED="\$(cat intersect_wa | uniq -u | wc -l)"

    echo "\${A_NAME}_TOT=\$A_TOTAL" >> ${data_file}
    echo "\${A_NAME}_INT=\$A_INTERSECTED" >> ${data_file}

    rm -rf FILE_A FILE_B intersect_wa

    touch "${source_file}.${process}.dummy"
    """
}

process generate_csv {
    
    input:
        file dummies from final_dummies.collect()
        val filter_hets from params.filter_hets

    script:

    data_file_name = "${dummies[0].getName().tokenize(".").get(0)}.data"
    data_file = file("${params.out_dir}/${data_file_name}")
    
    if(filter_hets == 1){
        csv_file_name = "${dummies[0].getName().tokenize(".").get(0)}-${params.intersect_window_a}_${params.intersect_window_b}_nohets.raw.csv"
    }else{
        csv_file_name = "${dummies[0].getName().tokenize(".").get(0)}-${params.intersect_window_a}_${params.intersect_window_b}.raw.csv"
    }
    
    csv_file = file("${params.out_dir}/${csv_file_name}")

    """
    cut -d '=' -f 1 < ${data_file} | awk -v RS= -v OFS=, '{\$1 = \$1} 1' > ${csv_file}
    cut -d '=' -f 2 < ${data_file} | awk -v RS= -v OFS=, '{\$1 = \$1} 1' >> ${csv_file}
    """
}