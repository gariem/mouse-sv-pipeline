params.strain = "C57BL_6NJ"

params.input_dir = "./data/input"
params.output_dir = "./data/merged"

params.min_coverage = '3'

params.combi_exe = './templates/combiSV2.1.pl'

Channel.fromPath("${params.input_dir}/${params.strain}-*vcf").set{vcf_files_ch}

vcf_files_ch.map{ file ->
    def strain = file.getName().tokenize(".").get(0).tokenize('-').get(0)
    return tuple(strain, file)
}
.groupTuple()
.set{ grouped_vcfs }

process merge_with_combi{

    publishDir file(params.output_dir), mode: "copy"

    input:
        set strain, file(vcf_files) from grouped_vcfs
        file combisv from file(params.combi_exe)
        val min_coverage from params.min_coverage

    output:
        file "${strain}-combisv_c*.vcf"

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
    perl ${combisv} -pbsv ${pbsv_file} -sniffles  ${sniffles_file} -cutesv ${cutesv_file} -o "${strain}-combisv_c${min_coverage}.vcf" -c ${min_coverage}
    """
}