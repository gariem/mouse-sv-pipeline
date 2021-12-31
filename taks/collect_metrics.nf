#!/usr/bin/env nextflow

params.input_files = './data/merged/DBA_2J-merged_100_1_30.vcf'
params.previous_dir = './data/previous'
params.validated_dir = './data/validated'

params.mappings = "./data/input/merge_mappings.txt"

params.out_dir = './data/reports/raw'

params.bcftools = '/home/egarcia/appdir/bcftools/bin/bcftools'

Channel.fromPath(params.input_files).set{input_files}

process general_metrics {

    publishDir file(params.out_dir), mode: "copy", pattern: "*.data"

    input:
        file vcf_file from input_files
        file mappings_file from file(params.mappings)

    output:
        file '*.bed' into bed_files

    script:

    vcf_name = vcf_file.getName().tokenize(".").get(0)
    strain = vcf_name.tokenize("-").get(0) 
    data_file = file("${params.out_dir}/${vcf_name}.data")

    """
    TOTAL_CALLS="\$(${params.bcftools} view --no-header ${vcf_file} | wc -l)"
    SINGLE_CALLER="\$(${params.bcftools} query -i"SUPP='1'" -f'%CHROM\\t%POS\\t%END\\n' ${vcf_file} | wc -l)"
    DOUBLE_CALLER="\$(${params.bcftools} query -i"SUPP='2'" -f'%CHROM\\t%POS\\t%END\\n' ${vcf_file} | wc -l)"

    echo "STRAIN=${strain}" > ${data_file}
    echo "TOTAL=\$TOTAL_CALLS" >> ${data_file}
    echo "SUP1=\$SINGLE_CALLER" >> ${data_file}
    echo "SUP2=\$DOUBLE_CALLER" >> ${data_file}

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


process type_metrics {

    input:
        file bed_file from bed_files.flatten()
    
    script:
    
    name_arr = bed_file.getName().tokenize(".")
    data_file = file("${params.out_dir}/${name_arr.get(0)}.data")
    type = name_arr.get(1)
    """
    TYPE_COUNT="\$(cat ${bed_file} | wc -l)"
    echo "${type}=\$TYPE_COUNT" >> ${data_file}
    """
}