#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

params.previous_dir = './data/previous'

// params.input = "./data/input/calls/{A_J,DBA_2J,C57BL_6NJ}-*.vcf"
// params.validated_files = './data/validated/simple/{A_J,DBA_2J,C57BL_6NJ}*.bed'
// params.previous_files = './data/previous/{A_J,DBA_2J,C57BL_6NJ}*.bed'

// params.input = "./data/input/calls/{DBA_2J,C57BL_6NJ}-*.vcf"
// params.validated_files = './data/validated/simple/{DBA_2J,C57BL_6NJ}*.bed'
// params.previous_files = './data/previous/{DBA_2J,C57BL_6NJ}*.bed'

params.input = "./data/input/calls/C57BL_6NJ-*.vcf"
params.validated_files = './data/validated/simple/C57BL_6NJ*.bed'
params.previous_files = './data/previous/C57BL_6NJ*.bed'

params.igv_workdir = '/media/egarcia/DataBank/mouse/igv_workfiles'

params.out_dir = './data/analysis/local/details'

params.dysgu_probs='20'
params.read_depth='30,40,50'
params.sv_types='INS,DEL,INV,DUP'

params.filter_hets=true
params.min_score=0.85

params.max_diff_b6n=20
params.min_score_b6n=0.2

params.screenshots_missed=false
params.screenshots_random=false
params.random_sample=20
params.create_figures=true
params.big_inversions=true

params.size_distribution_script='./templates/size_distribution.py'


include { 
    bed_from_vcf; 
    clean_regions; 
    filter_prob_dysgu; 
    clean_hets;
    big_inversions;
    filter_read_depth;
    filter_diff_allele_depth
} from './bedfiles/vcf_operations'

include { 
    intersect_features as find_intersected; 
    intersect_features as find_missed; 
    intersect_features as find_previous; 
    retrieve_validated_features; 
    retrieve_preivous_features;
    calc_intersect_stats;
    save_intersect_stats;
    calculate_scores
} from './bedfiles/intersect_operations'

include { 
    size_data_from_vcf; 
    size_data_from_previous; 
    calculate_size_distribution 
} from './figures/size_distribution'

include {
    take_screenshots as screenshot_missed;
    take_screenshots as screenshot_random;
    radomize_bed_file
} from './igv/igv_capture'


workflow {

    Channel.fromPath(params.input).set{vcf_channel}

    Channel.fromPath(params.validated_files).map{ file ->
        def strain = file.name.tokenize('.').get(0)
        def type = file.name.tokenize('.').get(1)
        return tuple(strain, type, file)
    }.set{validated_features}
    
    Channel.fromPath(params.previous_files).map{ file ->
        def strain = file.name.tokenize('.').get(0)
        def type = file.name.tokenize('.').get(1)
        return tuple(strain, type, file)
    }.set{previous_features}

    clean_regions(vcf_channel).multiMap{file -> 
        all: file
        homs_in: file
    }.set{vcf_files}

    if(params.filter_hets){      
        homs_only = clean_hets(vcf_files.homs_in)
        files = vcf_files.all.concat(homs_only)
    }else{
        files = vcf_files.all
    }

    files.branch {
        dysgu: it.name.contains("dysgu")
        pbsv: it.name.contains("pbsv")
        other: true
    }.set {callers}

    // Create new files from dysgu with specific probability thresholds
    probabilities = Channel.from(params.dysgu_probs).splitCsv().flatten()
    dysgu_files = filter_prob_dysgu(callers.dysgu, probabilities)

    // Create new files from pbsv with specific read depth

    callers.pbsv.multiMap{file ->
        dp_filter: file
        diff_ad_filter: file
    }.set{pbsv}

    read_depths = Channel.from(params.read_depth).splitCsv().flatten()
    filtered_dp = filter_read_depth(pbsv.dp_filter, read_depths)

    diffs_values = Channel.from([-16,3], [0,100], [-100,0])
    filtered_diff_ad = filter_diff_allele_depth(pbsv.diff_ad_filter, diffs_values)

    pbsv_files = filtered_dp.concat(filtered_diff_ad)

    // Join new dsygu filter vcfs with other vcfs 
    //TODO: join merged files too
    dysgu_files.concat(pbsv_files).concat(callers.other).multiMap { file ->
        intersect: file
        high_score: file
    }.set { all_vcfs }

    sv_types = Channel.from(params.sv_types).splitCsv().flatten()
    
    bed_from_vcf(all_vcfs.intersect, sv_types).multiMap{ file ->
        validate: file
        previous: file
    }.set{new_features}

    validated_features.combine(new_features.validate, by: [0, 1]).map { tuple_element ->
        return tuple(tuple_element[3], tuple_element[1], tuple_element[2], tuple_element[4])
    }.multiMap { file ->
        intersect: file
        outersect: file
    }.set {new_and_validated}

    previous_features.combine(new_features.previous, by: [0, 1]).map { tuple_element ->
        return tuple(tuple_element[3], tuple_element[1], tuple_element[2], tuple_element[4])
    }.set {new_and_previous}

    previous_intersected = find_previous(new_and_previous, '-wa', 30)
    validated_intersected = find_intersected(new_and_validated.intersect, '-wa', 30)

    previous_intersected.combine(validated_intersected, by: [0, 1]).set {data}
    scores = calculate_scores(data)
    
    // // filter intersect scores higher than param.min_score
    scores.filter {
        (it[0].contains("C57BL_6NJ") && Float.parseFloat(it[1].name.split('_')[1]) >= params.min_score_b6n && Float.parseFloat(it[1].name.split('_')[2]) <= params.max_diff_b6n)
    }.multiMap { file ->
        outersect: file
        save: file
    }.set {src_high_score_files}

    // src_high_score_files.save.view()
    
    // // save high scores to data/analysis/strain/simple_name
    save_intersect_stats(src_high_score_files.save).view()

    // // map files with high scores back to feature beds 
    // // outersect is used to map with original feature beds
    // // figures is used to map with original vcf file
    // src_high_score_files.outersect.map { file ->
    //     def simple_name = file.name.split('_-_')[0]
    //     def strain = simple_name.tokenize('-').get(0)
    //     return tuple(strain, simple_name, file)
    // }.multiMap{tuple ->
    //     outersect: tuple
    //     vcf_candidates: tuple
    // }.set {high_score_tuples}

    // high_score_tuples.outersect.combine(new_and_validated.outersect, by: 1).map { tuple_element ->
    //     return tuple(tuple_element[3], tuple_element[0], tuple_element[4], tuple_element[5])
    // }.multiMap{tuple ->
    //     out: tuple
    //     in: tuple
    // }.set { new_and_validated_high }

    // // new_and_validated_out.view() =>
    // // C57BL_6NJ.DEL, C57BL_6NJ-cutesv.all, C57BL_6NJ-cutesv.all__DEL.B.bed, C57BL_6NJ.DEL.A.bed

    // // *******************************
    // // prepare high score vcfs
    // // this section is very important
    // // *******************************
    // if(params.create_figures || params.screenshots_missed || params.screenshots_random || params.big_inversions) {
        
    //     all_vcfs.high_score.map{ file ->
    //         def simple_name = file.name.replace(".vcf", "")
    //         def strain = simple_name.tokenize('-').get(0)
    //         return tuple(strain, simple_name, file)
    //     }.set{ all_vcfs_candidates }

    //     high_score_tuples.vcf_candidates.combine(all_vcfs_candidates, by:1).map {tuple_element ->
    //         return tuple_element[4]
    //     }.multiMap{file ->
    //         fig_pacbio: file
    //         fig_ilumina: file
    //         big_invs: file
    //     }.set {high_score_vcfs}
    // }

    // // Create figures from high scoring VCF files
    // if(params.create_figures) {

    //     size_data_from_vcf(high_score_vcfs.fig_pacbio).map{file ->
    //         def simple_name = file.name.replace(".pacbio.csv", "")
    //         return tuple(simple_name, file)
    //     }set {pacbio_data}

    //     size_data_from_previous(high_score_vcfs.fig_ilumina).map{file ->
    //         def simple_name = file.name.replace(".ilumina.csv", "")
    //         return tuple(simple_name, file)
    //     }set {ilumina_data}

    //     calculate_size_distribution(pacbio_data.combine(ilumina_data, by: 0), file(params.size_distribution_script))
    // }

    // // Take screenshots from missed validated features
    // if(params.screenshots_missed){

    //     find_missed(new_and_validated_high.out, '-v', 30).map{file ->
    //         def simple_name = file.name.split('__')[0]
    //         def strain = simple_name.tokenize('-').get(0)
    //         return tuple(strain, 'captures/missed', simple_name, file)
    //     }.set{out_data}

    //     screenshot_missed(out_data, file(params.igv_workdir))
    // }

    // if(params.screenshots_random) {
    //     radomize_bed_file(new_and_validated_high.in.map{ tuple_element ->
    //         return tuple_element[2]
    //     }).map{ file ->
    //         def simple_name = file.name.split('__')[0]
    //         def strain = simple_name.tokenize('-').get(0)
    //         return tuple(strain, 'captures/random', simple_name, file)
    //     }.set{random_data}
        
    //     screenshot_random(random_data, file(params.igv_workdir))
    // }

    // if(params.big_inversions){
    //     big_inversions(high_score_vcfs.big_invs).map{file ->
    //         def simple_name = file.name.replace(".bed", "")
    //         def strain = simple_name.tokenize('-').get(0)
    //         return tuple(strain, simple_name, file)
    //     }.set{src_biginvs}
    // }
    
}