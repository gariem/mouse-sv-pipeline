#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

params.previous_dir = './data/previous'
params.input = "./data/input/DBA_2J-*.vcf"

params.validated_files = './data/validated/DBA_2J*.bed'

params.out_dir = './data/analysis/details'

params.filter_hets=1
params.dysgu_probs='50,65,70'


params.size_distribution_script='./templates/size_distribution.py'

include { bed_from_vcf; clean_regions; filter_prob_dysgu; clean_hets } from './bedfiles/vcf_process'
include { intersect_files; retrieve_validated_features } from './bedfiles/intersect_bed'
include { generate_len_csv_from_vcf; generate_len_csv_from_previous_tsv; calculate_size_distribution } from './figures/size_distribution'

workflow {

    Channel.fromPath(params.input).set{vcf_channel}

    clean_regions(vcf_channel).multiMap{file -> 
        all: file
        homs_in: file
    }.set{vcf_files}

    files = vcf_files.all

    if(params.filter_hets==1){      
        homs_only = clean_hets(vcf_files.homs_in)
        files = vcf_files.all.concat(homs_only)
    }

    files.multiMap{file ->
        source: file
        merge: file
    }.set {prepared_files}

    prepared_files.source.branch {
        dysgu: it.name.contains("dysgu")
        other: true
    }.set {callers}

    probs = Channel.from(params.dysgu_probs).splitCsv().flatten()
    
    dysgu_files = filter_prob_dysgu(callers.dysgu, probs)

    sv_types = Channel.from('INS,DEL,INV,DUP').splitCsv().flatten()

    new_features = bed_from_vcf(dysgu_files, sv_types, '.B')

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

    validated_features = retrieve_validated_features(validated_grouped, '.A')

    new_features.map{ file ->
        def parts = file.name.split('__')
        def strain =parts[0].split('-')[0]
        def type = parts[1].tokenize('.').get(0)
        def key = strain + '.' + type

        return tuple(key, parts[0], file)

    }.set{ x_features }

    validated_features.map{ file ->
        def key = file.name.replace(".A.bed", "")
        return tuple(key, file)
    }.set{ x_validated }

    intersect_files(x_features.combine(x_validated, by: 0), '-wa', 30).view()

    // dysgu_files.multiMap{file ->
    //     pacbio: file
    //     ilumina: file
    // }.set {figure_sizes}

    // pacbio_data = generate_len_csv_from_vcf(figure_sizes.pacbio)
    // ilumina_data = generate_len_csv_from_previous_tsv(figure_sizes.ilumina, file(params.previous_dir))

    // figures  = calculate_size_distribution(pacbio_data, ilumina_data, file(params.size_distribution_script))

    
}