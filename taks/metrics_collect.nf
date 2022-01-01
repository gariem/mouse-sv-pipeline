#!/usr/bin/env nextflow

params.input_files = './data/merged/DBA_2J-merged_100_1_20.vcf'
params.out_dir = './data/reports/raw'
params.previous_dir = './data/previous'
params.validated_dir = './data/validated'

params.merge_mappings = "./data/input/merge_mappings.txt"
params.metrics_mappings = "./data/input/metrics_mappings.txt"

params.intersect_window_a = '0'
params.intersect_window_b = '20'

params.bcftools = '/home/egarcia/appdir/bcftools/bin/bcftools'

Channel.fromPath(params.input_files).set{input_files}

process general_metrics {

    //publishDir file(params.out_dir), mode: "move", pattern: "*.data"

    input:
        file vcf_file from input_files
        file mappings_file from file(params.merge_mappings)

    output:
        file '*.bed' into unflattened_bed_files

    script:

    vcf_name = vcf_file.getName().tokenize(".").get(0)
    strain = vcf_name.tokenize("-").get(0) 
    data_file = file("${params.out_dir}/${vcf_name}.data")

    """
    mkdir -p ${file(params.out_dir)}
    TOTAL_CALLS="\$(${params.bcftools} view --no-header ${vcf_file} | wc -l)"
    SUP1="\$(${params.bcftools} query -i"SUPP='1'" -f'%CHROM\\t%POS\\t%END\\n' ${vcf_file} | wc -l)"
    SUP2="\$(${params.bcftools} query -i"SUPP='2'" -f'%CHROM\\t%POS\\t%END\\n' ${vcf_file} | wc -l)"

    echo "SRC_FILE=${vcf_file.getName()}" > ${data_file}
    echo "STRAIN=${strain}" >> ${data_file}
    echo "TOTAL=\$TOTAL_CALLS" >> ${data_file}
    echo "SUP1=\$SUP1" >> ${data_file}
    echo "SUP2=\$SUP2" >> ${data_file}
    echo "WINDOW_A=${params.intersect_window_a}" >> ${data_file}
    echo "WINDOW_B=${params.intersect_window_a}" >> ${data_file}

    for ((i=1;i<=19;i++)); 
    do 
        CHROM_COUNT="\$(${params.bcftools} query -i"CHROM='\$i'" -f'%CHROM\\t%POS\\t%END\\n' ${vcf_file} | wc -l)"
        echo "chr\$i=\$CHROM_COUNT" >> ${data_file}
    done

    CHROM_COUNT="\$(${params.bcftools} query -i"CHROM='X'" -f'%CHROM\\t%POS\\t%END\\n' ${vcf_file} | wc -l)"
    echo "chrX=\$CHROM_COUNT" >> ${data_file}

    while read -r line
    do
        TYPE="\$(echo \$line | cut -d ':' -f 1)"
        ${params.bcftools} query -i"SVTYPE='\$TYPE'" -f'%CHROM\\t%POS0\\t%END0\\t%SVLEN\\n' ${vcf_file} | \
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

    data_file = file("${params.out_dir}/${name_arr.get(0)}.data")

    """
    TYPE_COUNT="\$(cat ${bed_file} | wc -l)"
    echo "${type}=\$TYPE_COUNT" >> ${data_file}

    echo "\$(cat ${mappings_file} | grep ${type} || true;)" > compare_to.txt

    if [[ -s compare_to.txt ]]
    then
        while read -r line
        do
            VALIDATED_TYPE="\$(echo \$line | cut -d '=' -f 1)"
            VALIDATED_FILE=${validated_dir}/${strain}.\$VALIDATED_TYPE.mm39.bed

            if [[ -e \$VALIDATED_FILE ]]
            then
                echo "File exists: \$VALIDATED_FILE -> ${type}_\$VALIDATED_TYPE.${strain}.\$VALIDATED_TYPE.bed"
                cat \$VALIDATED_FILE | uniq -u > "${type}_\$VALIDATED_TYPE.A.${strain}.\$VALIDATED_TYPE.bed"
                cp ${bed_file} "${type}_\$VALIDATED_TYPE.B.${name}.${type}.bed"
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
    data_file=file("${params.out_dir}/${file_b.getName().tokenize(".").get(2)}.data")

    """
    A_TOTAL="\$(cat ${file_a} | wc -l)"
    A_NAME="\$(echo ${file_a} | cut -d '.' -f 4)"

    DATA_FILE="\$(echo ${file_b} | cut -d '.' -f 3).data"

    awk -F'\\t' 'BEGIN {OFS = FS} {print \$1,\$2-${window_a},\$3+${window_a},\$4}' ${file_a} > FILE_A
    awk -F'\\t' 'BEGIN {OFS = FS} {print \$1,\$2-${window_b},\$3+${window_b},\$4}' ${file_b} > FILE_B

    bedtools intersect -a FILE_A -b FILE_B -wa > intersect_wa

    A_INTERSECTED="\$(cat intersect_wa | uniq -u | wc -l)"

    echo "\${A_NAME}_TOTAL=\$A_TOTAL" >> ${data_file}
    echo "\${A_NAME}_INTERSECTED=\$A_INTERSECTED" >> ${data_file}

    touch "${file_b.getName().tokenize(".").get(2)}.${process}.dummy"
    """
}

process generate_csv {
    echo true
    
    input:
        file dummies from final_dummies.collect()

    script:

    data_file_name = "${dummies[0].getName().tokenize(".").get(0)}.data"
    data_file = file("${params.out_dir}/${data_file_name}")
    csv_file_name = "${dummies[0].getName().tokenize(".").get(0)}-${params.intersect_window_a}_${params.intersect_window_b}.raw.csv"
    csv_file = file("${params.out_dir}/${csv_file_name}")
    
    """
    cut -d '=' -f 1 < ${data_file} | awk -v RS= -v OFS=, '{\$1 = \$1} 1' > ${csv_file}
    cut -d '=' -f 2 < ${data_file} | awk -v RS= -v OFS=, '{\$1 = \$1} 1' >> ${csv_file}
    """
}