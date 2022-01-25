#!/usr/bin/env nextflow

params.strain = "*"
params.input_dir = "./data/input"
params.input_files = "${params.input_dir}/${params.strain}*.vcf"
params.out_dir = './data/reports/raw'
params.previous_dir = './data/previous'
params.validated_dir = './data/validated'

params.merge_mappings = "./data/input/merge_mappings.txt"
params.metrics_mappings = "./data/input/metrics_mappings.txt"

params.filter_hets = '0'
params.intersect_window = '20'

Channel.fromPath(params.input_files).set{input_files}


input_files.flatMap{ file ->
    def key = file.name.toString().tokenize('.').get(0)
    
    def filter_hets_arr = params.filter_hets.toString().split(',')
    def tuples = []

    for(filter_hets in filter_hets_arr){
        tuples.add(tuple(key,filter_hets,file))
    }

    return tuples
}
.groupTuple(by: [0, 1])
.set{ grouped_inputs }


process general_metrics {

    input:
        set key, filter, file(vcf_files) from grouped_inputs
        
        file merge_mappings_file from file(params.merge_mappings)
        file metrics_mappings_file from file(params.metrics_mappings)

    output:
        file "${vcf_name}.*.*" into prepared_files

    script:

    vcf_file = vcf_files
    vcf_name = vcf_file.getName().tokenize(".").get(0)
    strain = vcf_name.tokenize("-").get(0) 
    data_file = "${vcf_name}.data"

    filter_flag = (filter == '0' || vcf_name.contains("_all")) ? 0 : 1
    filter_name = filter == '1' ? "hom" : "all"

    prev_dir = file(params.previous_dir)
    validated_dir = file(params.validated_dir)

    """
    echo "Collecting general metrics for: ${vcf_file.getName()}"

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

    for ((i=1;i<=19;i++)); 
    do 
        if [ "${filter}" == "1" ]; then
            QUERY="CHROM='\$i' && GT='HOM'"
        else
            QUERY="CHROM='\$i'"
        fi

        CHROM_COUNT="\$(bcftools query -i"\$QUERY" -f'%CHROM\\t%POS\\t%END\\n' ${vcf_file} | wc -l)"
        echo "chr\$i=\$CHROM_COUNT" >> ${data_file}
    done

    if [ "${filter}" == "1" ]; then
        CHROM_COUNT="\$(bcftools query -i"CHROM='X' && GT='HOM'" -f'%CHROM\\t%POS\\t%END\\n' ${vcf_file} | wc -l)"
    else
        CHROM_COUNT="\$(bcftools query -i"CHROM='X'" -f'%CHROM\\t%POS\\t%END\\n' ${vcf_file} | wc -l)"
    fi
    
    echo "chrX=\$CHROM_COUNT" >> ${data_file}

    echo "Generating BED files by type"

    while read -r mappings_line
    do
        TYPE="\$(echo \$mappings_line | cut -d ':' -f 1)"

        if [ "${filter}" == "1" ]; then
            QUERY="SVTYPE='\$TYPE' && GT='HOM'"
        else
            QUERY="SVTYPE='\$TYPE' "
        fi

        bcftools query -i"\$QUERY" -f'%CHROM\\t%POS0\\t%END0\\t%SVLEN\\n' ${vcf_file} | \
            awk -F'\\t' 'BEGIN {OFS = FS} \$1 ~/^[0-9]*\$|^X\$/{print \$1,\$2,\$3,\$4}' > "Temp.${vcf_name}.\$TYPE.bed"

        PREV_FILE="${prev_dir}/${strain}.\$TYPE.mm39.bed"

        TYPE_COUNT="\$(cat Temp.${vcf_name}.\$TYPE.bed | wc -l)"
        PREV_TYPE_COUNT="\$(cat \$PREV_FILE | wc -l)"

        echo "\${TYPE}=\$TYPE_COUNT" >> ${data_file}
        echo "\${TYPE}0=\$PREV_TYPE_COUNT" >> ${data_file}

        echo "\$(cat ${metrics_mappings_file} | grep \${TYPE} || true;)" > compare_to.txt

        if [[ -s compare_to.txt ]]
        then
            while read -r compare_line
            do
                VALIDATED_TYPE="\$(echo \$compare_line | cut -d '=' -f 1)"
                VALIDATED_FILE=${validated_dir}/${strain}.\$VALIDATED_TYPE.mm39.bed

                if [[ -e \$VALIDATED_FILE ]]
                then
                    ln -s \$VALIDATED_FILE "${vcf_name}.${filter_name}.\${TYPE}.\${VALIDATED_TYPE}.A.${strain}.\${VALIDATED_TYPE}.bed"
                    ln -s "Temp.${vcf_name}.\$TYPE.bed" "${vcf_name}.${filter_name}.\${TYPE}.\$VALIDATED_TYPE.B.${vcf_name}.\${TYPE}.bed"
                    ln -s ${data_file} "${vcf_name}.${filter_name}.\${TYPE}.\${VALIDATED_TYPE}.data"
                fi
            done < compare_to.txt
        fi

    done < ${merge_mappings_file}
    """
}

prepared_files.flatten().flatMap{ file ->
    def parts = file.name.toString().tokenize('.')
    key = parts.get(0) + '.' + parts.get(1) + '.' + parts.get(2) + '.' + parts.get(3)

    def intersect_w_arr = params.intersect_window.toString().split(',')
    def tuples = []

    for(intersect_w in intersect_w_arr){
        tuples.add(tuple(key,intersect_w,file))
    }

    return tuples
}
.groupTuple(by: [0, 1], size: 3)
.set{ grouped_files }

process intersect_files {

    input:
        set key, window, file(input_files) from grouped_files

    output:
        file '*.metrics' into metrics_files

    script:
 

    for(int i = 0;i<=2;i++) {
        vcfName = input_files.get(i).getName()
        if(vcfName.contains('.A.') && vcfName.endsWith('.bed')){
            file_a = input_files.get(i)
            continue
        }
        if(vcfName.contains('.B.') && vcfName.endsWith('.bed')){
            file_b = input_files.get(i)
            continue
        }
        if(vcfName.endsWith('.data')){
            src_data_file = input_files.get(i)
            continue
        }
    }
    
    data_file=key + "." + window + ".metrics"

    """
    awk -F'\\t' 'BEGIN {OFS = FS} {print \$1,\$2-${window},\$3+${window},\$4}' ${file_a} | awk '!x[\$0]++' > FILE_A
    awk -F'\\t' 'BEGIN {OFS = FS} {print \$1,\$2-${window},\$3+${window},\$4}' ${file_b} | awk '!x[\$0]++' > FILE_B

    A_TOTAL="\$(cat FILE_A | wc -l)"
    A_NAME="\$(echo ${file_a} | cut -d '.' -f 4)"

    bedtools intersect -a FILE_A -b FILE_B -wa > intersect_wa

    A_INTERSECTED="\$(cat intersect_wa | uniq -u | wc -l)"

    cat ${src_data_file} > ${data_file}

    echo "W_AB=${window}" >> ${data_file}
    echo "T_\${A_NAME}=\$A_TOTAL" >> ${data_file}
    echo "I_\${A_NAME}=\$A_INTERSECTED" >> ${data_file}

    rm -rf FILE_A FILE_B intersect_wa
    """
}

metrics_files.map{ file ->
    def parts = file.name.toString().tokenize('.')
    key = parts.get(0) + '.' + parts.get(1) + '.' + parts.get(4)
    return tuple(key, file)
}
.groupTuple(size: 4)
.set{ metrics_grouped }

process join_metrics() {

    publishDir file("${params.out_dir}"), mode: "move", pattern: "*.csv"

    input:
        set key, file(metrics) from metrics_grouped
    
    output:
        file '*.csv' into final_metrics

    """
    cat ${metrics} | awk '!x[\$0]++' > tmp
    #cat ${metrics} | sort | uniq > tmp

    cut -d '=' -f 1 < tmp | awk -v RS= -v OFS=, '{\$1 = \$1} 1' > "${key}.csv"
    cut -d '=' -f 2 < tmp | awk -v RS= -v OFS=, '{\$1 = \$1} 1' >> "${key}.csv"

    rm tmp
    """
}