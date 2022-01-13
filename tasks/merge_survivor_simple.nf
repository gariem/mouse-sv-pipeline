params.strain = "*"

params.input_dir = "./data/input"
params.output_dir = "./data/merged"

params.max_dist = '100'
params.min_callers = '1'
params.same_type = '1'
params.min_size = '30'

Channel.fromPath("${params.input_dir}/${params.strain}-*.vcf").set{vcf_files_ch}

// Transform channel, group by caller (from file name)
vcf_files_ch.flatMap{ file ->
    def strain = file.getName().tokenize(".").get(0).tokenize('-').get(0)
    
    def max_dist_arr = params.max_dist.toString().split(',')
    def min_size_arr = params.min_size.toString().split(',')
    def min_callers_arr = params.min_callers.toString().split(',')

    def tuples = []

    for(max_dist in max_dist_arr){
        for(min_callers in min_callers_arr){
            for(min_size in min_size_arr){
                tuples.add(tuple(strain, max_dist, min_callers, min_size, file))
            }
        }
    }
    return tuples
}
.groupTuple(by: [0, 1, 2, 3])
.set{ grouped_vcfs }


process merge_with_survivor{

    publishDir file(params.output_dir), mode: "move"

    input:
        set strain, max_dist, min_callers, min_size, file(vcf_files) from grouped_vcfs

    output:
        file "${strain}-survivor_*.vcf"
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

    echo "${pbsv_file}" > '${strain}-merged-inputlist.txt'
    echo "${sniffles_file}" >> '${strain}-merged-inputlist.txt'
    echo "${cutesv_file}" >> '${strain}-merged-inputlist.txt'

    SURVIVOR merge '${strain}-merged-inputlist.txt' ${max_dist} ${min_callers} ${params.same_type} 1 0 ${min_size} "${strain}-survivor_${max_dist}_${min_callers}_${min_size}_all.vcf"

    rm -rf '${strain}-merged-inputlist.txt'
    """
}