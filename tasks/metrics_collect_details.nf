#!/usr/bin/env nextflow

nextflow.enable.dsl = 2


params.file= "DBA_2J-mapped_200_2_20_all.vcf"
params.input_dir = "./data/merged"
params.input_files = "${params.input_dir}/${params.file}"
params.out_dir = './data/reports/raw'
params.previous_dir = './data/previous'
params.validated_dir = './data/validated'

params.merge_mappings = "./data/input/merge_mappings.txt"
params.metrics_mappings = "./data/input/metrics_mappings.txt"
params.genome_sizes = "./data/input/genome_size.txt"

params.size_distribution_script='./templates/size_distribution.py'


process generate_new_data_file {

    input:
        file vcf_file

    output:
        file '*.new_data.csv'

    script:
    file_id = vcf_file.name.toString().tokenize('.').get(0)

    """
    bcftools query -f'%CHROM,%SVTYPE,%SVLEN\\n' ${vcf_file} | awk -F',' 'BEGIN {OFS = FS} {abs=\$3;sub("^-", "",abs); print \$1,\$2,abs}' > "${file_id}.new_data.csv"
    """
}

process generate_previous_data_file {

    echo true

    input:
        file vcf_file
        file previous_dir

    output:
        file '*.previous_data.csv'

    script:

    file_id = vcf_file.name.toString().tokenize('.').get(0)
    strain = file_id.tokenize('-').get(0)
    
    """
    for file in ${previous_dir}/${strain}*; do
        TYPE="\$(echo \$file | cut -d '.' -f 2)"
        cat \$file | awk -v type=\$TYPE '{print \$1","type","\$4}' >> "${file_id}.previous_data.csv"
    done

    """
}

process calculate_size_distribution {

    input:
        tuple val(file_key), file(data_files)
        file size_distribution_script
    
    output:
        file "${file_key}.sizes.png"

    script:
    
    if(data_files.get(0).getName().contains("new_data")){
        new_data = data_files.get(0)
        previous_data = data_files.get(1)
    }else{
        new_data = data_files.get(1)
        previous_data = data_files.get(0)
    }

    """
    python ${size_distribution_script} -i1 ${new_data} -i2 ${previous_data} -o "${file_key}.sizes.png"
    """
}

workflow {

    Channel.fromPath(params.input_files).multiMap{file ->
        new_data: file
        previous_data: file
    }.set {input}

    new_data = generate_new_data_file(input.new_data)
    previous_data = generate_previous_data_file(input.previous_data, file(params.previous_dir))

    new_data.concat(previous_data).map{ file -> 
        def key = file.name.toString().tokenize('.').get(0)
        return tuple(key,file)
    }
    .groupTuple(size: 2)
    .set{ grouped_data }

    calculate_size_distribution(grouped_data, file(params.size_distribution_script))

    
}