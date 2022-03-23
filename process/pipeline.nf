#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

params.previous_dir = './data/previous'

// params.pattern = "{A_J,DBA_2J,C57BL_6NJ,C3H_HeJ,AKR_J,C57BL_6JEve}"
params.pattern = "{A_J,DBA_2J,C57BL_6NJ,C57BL_6JEve}"
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
    split_data as split_0_100;
    split_data as split_100_1K;
    split_data as split_1K_100K;
    split_data as split_30_100;
    calc_overlaps;
    draw_overlaps;
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
    }
    .set {minigraph_beds}

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

    // raw_scores.filter {
    //     (it[1].contains("C57BL_6") && Float.parseFloat(it[2].name.split('_')[1]) >= params.min_score_b6n && Float.parseFloat(it[2].name.split('_')[2]) <= params.max_diff_b6n) ||
    //     (Float.parseFloat(it[2].name.split('_')[1]) >= params.min_score && Float.parseFloat(it[2].name.split('_')[2]) <= params.max_diff)
    // }.multiMap{file -> 
    //     raw: file
    //     calc: file
    // }.set{scores}

    // scores.calc.map{tuple_element ->
    //     def name = tuple_element[2].name
    //     def parts = name.tokenize("_")
    //     def type = parts.get(0)
    //     def score = parts.get(1)
    //     def diff = parts.get(2)
    //     return tuple(tuple_element[0], tuple_element[1], type, score, diff)
    // }.groupTuple(by: [0, 2]).map{ items ->
    //     def tuples = []
    //     for(int i=0; i<items[1].size(); i++) {
    //         tuples.add(tuple(items[1].get(i), items[3].get(i), items[4].get(i)))
    //     }
    //     tuples.sort { a,b -> b[1] <=> a[1] ?: a[2] <=> b[2] }
    //     return tuple(items[0], items[2], tuples.findAll{t -> t[0].contains("pbsv.all.ad_")}) // Modify here to find other scores for pbsv or other callers
    // }.view()

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

    // src_high_score_files.save.view()

    // save high scores to data/analysis/strain/simple_name
    save_intersect_stats(src_high_score_files.save)

    // Prepate tuples to compare between callers (minigraph vs pbsv)
    src_high_score_files.intercaller.filter {
        it[0].contains("minigraph") || it[0].contains("pbsv")
    }.set {src_interstrain_highscores}

    src_interstrain_highscores.combine(new_features.intercaller.map{ tuple_element -> 
        return tuple(tuple_element[2],tuple_element[0], tuple_element[1], tuple_element[3])
    }.filter{it[2] == 'INS' || it[2] == 'DEL'}, by: 0).map{tuple_element ->
        return tuple(tuple_element[2], tuple_element[3], tuple_element[4])
    }.groupTuple(by: [0, 1], size: 2).set {src_intercaller_data}

    // src_intercaller_data.view()
    split_data = split_0_100(src_intercaller_data, 0, 100).concat(split_100_1K(src_intercaller_data, 100, 1000)).concat(split_1K_100K(src_intercaller_data, 1000, 100000)).concat(split_30_100(src_intercaller_data, 30, 100))
    draw_overlaps(calc_overlaps(split_data).collect()).view()

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