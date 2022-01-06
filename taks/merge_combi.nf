params.strain = "C57BL_6NJ"

params.input_dir = "./data/input"
params.output_dir = "./data/merged"

params.min_coverage = '3'

params.combi_exe = './templates/combiSV2.1.pl'

Channel.fromFilePairs("${params.input_dir}/${params.strain}*-{sniffles,pbsv}.vcf").set{samples_ch}

process merge_with_combi{

    publishDir file(params.output_dir), mode: "copy"

    input:
        set strain, file(vcf_files) from samples_ch
        file combisv from file(params.combi_exe)
        val min_coverage from params.min_coverage

    output:
        file "${strain}-combisv_c*.vcf"
    script:
    
    if(vcf_files.get(0).getName().contains("pbsv")){
        pbsv_file = vcf_files.get(0)
        sniffles_file = vcf_files.get(1)
    }else{
        pbsv_file = vcf_files.get(1)
        sniffles_file = vcf_files.get(0)
    }

    """
    perl ${combisv} -pbsv ${pbsv_file} -sniffles  ${sniffles_file} -o "${strain}-combisv_c${min_coverage}.vcf" -c ${min_coverage}
    """
}