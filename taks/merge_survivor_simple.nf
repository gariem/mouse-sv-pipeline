params.strain = "DBA_2J"

params.input_dir = "./data/input"
params.output_dir = "./data/merged"

params.max_dist = '100'
params.min_callers = '1'
params.same_type = '1'
params.min_size = '30'

params.survivor = '/home/egarcia/workspace/github/SURVIVOR/Debug/SURVIVOR'

Channel.fromPath("${params.input_dir}/${params.strain}-*.vcf").set{vcf_files_ch}


// Transform channel, group by caller (from file name)
vcf_files_ch.map{ file ->
    def strain = file.getName().tokenize(".").get(0).tokenize('-').get(0)
    return tuple(strain, file)
}
.groupTuple()
.set{ grouped_vcfs }


process merge_with_survivor{

    publishDir file(params.output_dir), mode: "copy"

    input:
        set strain, file(vcf_files) from grouped_vcfs

    output:
        file "${strain}-unmapped_*.vcf"
    script:
    
    for(int i = 0;i<=2;i++) {
        if(vcf_files.get(i).getName().contains("pbsv")){
            pbsv_file = vcf_files.get(i)
            continue
        }
        if(vcf_files.get(i).getName().contains("cutesv")){
            cutesv_file = vcf_files.get(i)
            continue
        }
        if(vcf_files.get(i).getName().contains("sniffles")){
            sniffles_file = vcf_files.get(i)
            continue
        }
    }

    """

    echo "${pbsv_file}" > '${params.strain}-merged-inputlist.txt'
    echo "${sniffles_file}" >> '${params.strain}-merged-inputlist.txt'
    echo "${cutesv_file}" >> '${params.strain}-merged-inputlist.txt'

    ${params.survivor} merge '${params.strain}-merged-inputlist.txt' ${params.max_dist} ${params.min_callers} ${params.same_type} 1 0 ${params.min_size} "${params.strain}-unmapped_${params.max_dist}_${params.min_callers}_${params.min_size}.vcf"

    """
}