#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

params.previous_dir = './data/previous'
params.input = "./data/input/calls/{DBA_2J,C57BL_6NJ}-*.vcf"

params.validated_files = './data/validated/{DBA_2J,C57BL_6NJ}*.bed'

params.out_dir = './data/analysis/details'

params.dysgu_probs='50,65,70'
params.sv_types='INS,DEL,INV,DUP'

params.filter_hets=true
params.min_score=0.8
params.take_screenshots=true


params.size_distribution_script='./templates/size_distribution.py'

include { bed_from_vcf; clean_regions; filter_prob_dysgu; clean_hets } from './bedfiles/vcf_operations'
include { intersect_features; retrieve_validated_features; save_intersect_stats } from './bedfiles/intersect_operations'
include { generate_len_csv_from_vcf; generate_len_csv_from_previous_tsv; calculate_size_distribution } from './figures/size_distribution'


workflow {

    Channel.fromPath(params.input).set{vcf_channel}

    clean_regions(vcf_channel).multiMap{file -> 
        all: file
        homs_in: file
    }.set{vcf_files}

    files = vcf_files.all

    if(params.filter_hets){      
        homs_only = clean_hets(vcf_files.homs_in)
        files = vcf_files.all.concat(homs_only)
    }

    files.multiMap{file ->
        source: file
        merge: file
    }.set {prepared_files}

    // use [[prepared_files.source]] to merge and later inject the results back in the pipeline for further processing

    prepared_files.source.branch {
        dysgu: it.name.contains("dysgu")
        other: true
    }.set {callers}

    // Create new files from dysgu with specific probability thresholds
    probabilities = Channel.from(params.dysgu_probs).splitCsv().flatten()
    dysgu_files = filter_prob_dysgu(callers.dysgu, probabilities)

    // Join new dsygu filter vcfs with other vcfs 
    //TODO: join merged files too
    all_vcfs=dysgu_files.concat(callers.other)

    sv_types = Channel.from(params.sv_types).splitCsv().flatten()
    src_new_features = bed_from_vcf(all_vcfs, sv_types, '.B')

    Channel.fromPath(params.validated_files).set{validated}

    validated.map { file ->
        def identifier = file.name.replace("mm39.bed", "")
        def strain = identifier.tokenize('.').get(0)
        def type = identifier.contains('H1') || identifier.contains('H2') ? 'DEL' : 'INS'
        def key = strain + '.' + type
        return tuple(key, file)
    }
    .groupTuple(size: 2)
    .set{ validated_grouped }

    src_validated_features = retrieve_validated_features(validated_grouped, '.A')

    src_new_features.map{ file ->
        def parts = file.name.split('__')
        def strain =parts[0].split('-')[0]
        def type = parts[1].tokenize('.').get(0)
        def key = strain + '.' + type

        return tuple(key, parts[0], file)

    }.set{ new_features }

    src_validated_features.map{ file ->
        def key = file.name.replace(".A.bed", "")
        return tuple(key, file)
    }.set{ validated_features }

    intersect_results = intersect_features(new_features.combine(validated_features, by: 0), '-wa', 30)

    intersect_results.map { file ->
        def simple_name = file.name.split('__')[0]
        return tuple(simple_name, file)
    }.groupTuple()
    .set{ grouped_intersect_result }

    save_intersect_stats(grouped_intersect_result, 'validated').view()

    // if(params.take_screenshots){

    // }

    // dysgu_files.multiMap{file ->
    //     pacbio: file
    //     ilumina: file
    // }.set {figure_sizes}

    // pacbio_data = generate_len_csv_from_vcf(figure_sizes.pacbio)
    // ilumina_data = generate_len_csv_from_previous_tsv(figure_sizes.ilumina, file(params.previous_dir))

    // figures  = calculate_size_distribution(pacbio_data, ilumina_data, file(params.size_distribution_script))

    
}