#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

params.previous_dir = './data/previous'

params.pattern = "{A_J,DBA_2J,C57BL_6NJ,C3H_HeJ,AKR_J,C57BL_6JEve}"
// params.pattern = "{A_J,DBA_2J,C57BL_6NJ,C57BL_6JEve}"
// params.pattern = "{DBA_2J,C57BL_6NJ}"
// params.pattern = "C57BL_6NJ"
// params.pattern = "DBA_2J"

params.input = "./data/input/calls/${params.pattern}-*.vcf"
params.validated_files = "./data/validated/simple/${params.pattern}*.bed"
params.previous_files = "./data/previous/${params.pattern}*.bed"

params.graph_files = "./data/input/minigraph/15plusEve/${params.pattern}.bed"
params.full_graph_files = './data/input/minigraph/full'

params.igv_workdir = '/media/egarcia/DataBank/mouse/igv_workfiles'

params.out_dir = './data/analysis/local/details'

params.dysgu_probs='20'
params.read_depth='30,40,50'
params.sv_types='INS,DEL,INV,DUP'

params.filter_hets=true
params.min_score=0.75
params.max_diff=10

params.max_diff_b6n=60
params.min_score_b6n=0.50

params.screenshots_missed=false
params.screenshots_random=false
params.random_sample=20
params.create_figures=true
params.big_inversions=true

params.size_distribution_script='./templates/size_distribution.py'

include { 
    bed_from_vcf; 
    bed_from_vcf as bed_from_survivor;
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
    size_data_from_bedfiles; 
    filtered_size_data_from_bedfiles;
    calculate_size_distribution;
    calculate_size_distribution_filtered
} from './figures/size_distribution'

include {
    split_triple_data;
    split_survivor_data;
    calc_overlaps as overlaps_0;
    calc_overlaps as overlaps_5;
    calc_overlaps as overlaps_10;
    calc_survivor_scores;
    summarize_overlaps;
    summarize_survivor;
} from './figures/overlaps'

include {
    take_screenshots as screenshot_missed;
    take_screenshots as screenshot_random;
    radomize_bed_file
} from './igv/igv_capture'

include {
    bed_from_full_graph;
    rename_tuples;
    intersect_all_minigraph;
    get_features_across_strains
} from './graph/minigraph'

include {
    survivor_vcf_from_bed;
    join_survivor_vcfs;
    merge_survivor_vcfs
} from './merge/survivor'


workflow {

    Channel.fromPath(params.input).set{vcf_channel}

    Channel.fromPath(params.graph_files).set{graph_channel}

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
    // TODO: join merged files too
    dysgu_files.concat(pbsv_files).concat(callers.other).multiMap { file ->
        intersect: file
        high_score: file
    }.set { all_vcfs }

    // Bed files must end with __TYPE.bed
    bed_from_full_graph(graph_channel).flatten().map{ file ->
        def simple_name = file.name.split('__')[0]
        def strain = simple_name.tokenize('-').get(0)
        def type = file.name.split('__')[1].replace(".bed", "")
        return tuple(strain, type, simple_name, file)
    }.multiMap{file ->
        concat: file
        compare: file
    }.set {minigraph_beds}

    sv_types = Channel.from(params.sv_types).splitCsv().flatten()

    bed_from_vcf(all_vcfs.intersect, sv_types).concat(minigraph_beds.concat).multiMap{ file ->
        validate: file
        previous: file
        intercaller: file
    }.set{new_features}

    validated_features.combine(new_features.validate, by: [0, 1]).map { tuple_element ->
        return tuple(tuple_element[3], tuple_element[1], tuple_element[2], tuple_element[4])
    }.multiMap { file ->
        intersect: file
        outersect: file
        intercaller: file
    }.set {new_and_validated}

    previous_features.combine(new_features.previous, by: [0, 1]).map { tuple_element ->
        return tuple(tuple_element[3], tuple_element[1], tuple_element[2], tuple_element[4])
    }.multiMap { file ->
        intersect: file
        figures: file
    }.set {new_and_previous}

    previous_intersected = find_previous(new_and_previous.intersect, '-wa')
    validated_intersected = find_intersected(new_and_validated.intersect, '-wa')

    previous_intersected.combine(validated_intersected, by: [0, 1]).set {data}
    scores = calculate_scores(data)

    scores.filter {
        (   
            (it[1].contains("C57BL_6") && Float.parseFloat(it[2].name.split('_')[1]) >= params.min_score_b6n && Float.parseFloat(it[2].name.split('_')[2]) <= params.max_diff_b6n) 
            || (Float.parseFloat(it[2].name.split('_')[1]) >= params.min_score && Float.parseFloat(it[2].name.split('_')[2]) <= params.max_diff)
        )
            && 
        ( 
            it[1].contains("minigraph") || it[1].contains("pbsv.all.ad_")
        ) 
        
    }.map{tuple_element ->
        return tuple(tuple_element[1], tuple_element[2])
    }.groupTuple(by: 0, size: 2).multiMap { file ->
        outersect: file
        save: file
        minigraph: file
        intercaller: file
    }.set {src_high_score_files}

    // save high scores to data/analysis/strain/simple_name
    save_intersect_stats(src_high_score_files.save)

    // Prepate tuples to compare between callers (minigraph vs pbsv) and validated data
    src_high_score_files.intercaller.filter {
        it[0].contains("minigraph") || it[0].contains("pbsv")
    }.set {src_interstrain_highscores}


    src_interstrain_highscores.combine(new_and_validated.intercaller.map{ tuple_element -> 
        def strain = tuple_element[0].tokenize("-").get(0)
        return tuple(tuple_element[0], strain, tuple_element[1], tuple_element[2], tuple_element[3])
    }.filter{it[2] == 'INS' || it[2] == 'DEL'}, by: 0).map{tuple_element ->
        return tuple(tuple_element[2], tuple_element[3], tuple_element[5], tuple_element[4])
    }.groupTuple(by: [0, 1, 3], size: 2).set {src_intercaller_data}

    survivor_vcf_from_bed(src_intercaller_data).flatMap { item ->
        def tuples = []
        for (file in item[2]){
            def caller = file.name.tokenize('-').get(1)
            tuples.add(tuple(item[0], caller, item[1], file))
        }
        return tuples
    }.groupTuple(by: [0, 1]).set{survivor_vcfs}    

    survivor_types = Channel.from('INS', 'DEL').splitCsv().flatten()
    survivor_features = bed_from_survivor(merge_survivor_vcfs(join_survivor_vcfs(survivor_vcfs).groupTuple(by: 0)), survivor_types)

    validated_features.combine(survivor_features, by: [0, 1]).map { tuple_element ->
        return tuple(tuple_element[0], tuple_element[1], tuple_element[2], tuple_element[4])
    }.set {survivor_and_validated}

    split_survivor_data(survivor_and_validated).flatMap{tuple_element->
        tuples = []
        for (file in tuple_element[2]){
            def range = file.name.tokenize('.').get(2)
            tuples.add(tuple(tuple_element[0], tuple_element[1], range, file))
        }
        for (file in tuple_element[3]){
            def range = file.name.tokenize('.').get(2)
            tuples.add(tuple(tuple_element[0], tuple_element[1], range, file))
        }
        return tuples
    }.groupTuple(by: [0, 1, 2]).set {split_survivor_data}

    summarize_survivor(calc_survivor_scores(split_survivor_data, 10).collect()).view()

    split_triple_data(src_intercaller_data).flatMap{tuple_element->
        tuples = []
        for (file in tuple_element[2]){
            def range = file.name.tokenize('.').get(2)
            tuples.add(tuple(tuple_element[0], tuple_element[1], range, file))
        }
        for (file in tuple_element[3]){
            def range = file.name.tokenize('.').get(2)
            tuples.add(tuple(tuple_element[0], tuple_element[1], range, file))
        }   
        for (file in tuple_element[4]){
            def range = file.name.tokenize('.').get(2)
            tuples.add(tuple(tuple_element[0], tuple_element[1], range, file))
        }   
        return tuples
    }.groupTuple(by: [0, 1, 2]).set {split_data}

    summarize_overlaps(overlaps_0(split_data, 0).concat(overlaps_5(split_data, 5)).concat(overlaps_10(split_data, 10)).collect())

    // Transform high score tuples 
    src_high_score_files.minigraph.filter {
        it[0].contains("minigraph")
    }.flatMap{ tuple_element ->
        def tuples = []
        tuple_element[1].each { elem -> tuples.add(tuple(tuple_element[0], elem.name.tokenize('_').get(0), elem)) }
        return tuples
    }.set {src_minigraph_highscores}

    // Combine high score minigraph with original bedfiles

    src_minigraph_highscores.combine(
        minigraph_beds.compare.map{ tuple_element ->
            return tuple(tuple_element[2], tuple_element[1], tuple_element[3])
        }, by: [0, 1]
    ).map{tuple_element ->
        return tuple(tuple_element[1], tuple_element[3])
    }.groupTuple(by: 0).set {high_score_minigraph_data}

    repeated = get_features_across_strains(intersect_all_minigraph(high_score_minigraph_data))
    
    // map files with high scores back to feature beds 
    // outersect is used to map with original feature beds
    // figures is used to map with original vcf file

    src_high_score_files.outersect.combine(new_and_validated.outersect, by: 0).map { tuple_element ->
        return tuple(tuple_element[0], tuple_element[2], tuple_element[3], tuple_element[4])
    }.multiMap{tuple ->
        out: tuple
        in: tuple
        figs: tuple
    }.set { new_and_validated_high }

    // Take screenshots from missed validated features
    if(params.screenshots_missed){

        find_missed(new_and_validated_high.out, '-v').map{tuple_element ->
            def simple_name = tuple_element[0]
            def strain = simple_name.tokenize('-').get(0)
            return tuple(strain, 'captures/missed', simple_name, tuple_element[2])
        }.set{out_data}

        screenshot_missed(out_data, file(params.igv_workdir))
    }

    // Take screenshots from random predicted features

    if(params.screenshots_random) {

        radomize_bed_file(new_and_validated_high.in.map{ tuple_element ->
            return tuple_element[3]
        }).map{ file ->
            def simple_name = file.name.split('__')[0]
            def strain = simple_name.tokenize('-').get(0)
            return tuple(strain, 'captures/random', simple_name, file)
        }.set{random_data}
        
        screenshot_random(random_data, file(params.igv_workdir))
    }


    // Create figures from high scoring files
    if(params.create_figures) {

        new_and_validated_high.figs.map {tuple_element ->
            return tuple(tuple_element[0], tuple_element[1])
        }.combine(new_and_previous.figures, by: [0, 1])
        .multiMap{ file -> 
            all: file
            filtered: file
        }.set{bed_size_data}

        size_data = size_data_from_bedfiles(bed_size_data.all).groupTuple(by: 0).map{tuple_element ->
            return tuple(tuple_element[0], tuple_element[2], tuple_element[3])
        }

        calculate_size_distribution(size_data, file(params.size_distribution_script))
        
        filtered_size_data = filtered_size_data_from_bedfiles( bed_size_data.filtered.combine(repeated, by: 1)).groupTuple(by: 0).map{tuple_element ->
            return tuple(tuple_element[0], tuple_element[2], tuple_element[3])
        }

        // filtered_size_data.view()
        // size_data.view()

        calculate_size_distribution_filtered(filtered_size_data, file(params.size_distribution_script))

    }
    
}