params.strain = "DBA_2J"

params.input_dir = "./data/input"
params.output_dir = "./data/merged"

params.max_dist = '100'
params.min_callers = '1'
params.same_type = '1'
params.min_size = '30'

params.survivor = '/home/egarcia/workspace/github/SURVIVOR/Debug/SURVIVOR'

Channel.fromFilePairs("${params.input_dir}/${params.strain}*-{sniffles,pbsv}.vcf").set{samples_ch}

process merge_with_combi{

    publishDir file(params.output_dir), mode: "copy"

    input:
        set strain, file(vcf_files) from samples_ch

    output:
        file "${strain}-unmapped_*.vcf"
    script:
    
    if(vcf_files.get(0).getName().contains("pbsv")){
        pbsv_file = vcf_files.get(0)
        sniffles_file = vcf_files.get(1)
    }else{
        pbsv_file = vcf_files.get(1)
        sniffles_file = vcf_files.get(0)
    }

    """

    echo "${pbsv_file}" > '${params.strain}-merged-inputlist.txt'
    echo "${sniffles_file}" >> '${params.strain}-merged-inputlist.txt'

    ${params.survivor} merge '${params.strain}-merged-inputlist.txt' ${params.max_dist} ${params.min_callers} ${params.same_type} 1 0 ${params.min_size} "${params.strain}-unmapped_${params.max_dist}_${params.min_callers}_${params.min_size}.vcf"

    """
}