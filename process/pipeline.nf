#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { bed_from_vcf; clean_regions; filter_prob_dysgu; clean_hets } from './bedfiles/vcf_process'
include { intersect_files } from './bedfiles/intersect_bed'
include { generate_len_csv_from_vcf; generate_len_csv_from_previous_tsv; calculate_size_distribution } from './figures/size_distribution'


params.previous_dir = './data/previous'
params.input = "./data/input/C57BL_6NJ-*.vcf"

params.filter_hets=1
params.dysgu_probs='60,70,80'


params.size_distribution_script='./templates/size_distribution.py'

workflow {

    Channel.fromPath(params.input).set{input}

    clean_regions(input).multiMap{file -> 
        all: file
        homs: file
    }.set{input_files}

    files = input_files.all

    if(params.filter_hets==1){      
        homs_only = clean_hets(input_files.homs)
        files = input_files.all.concat(homs_only)
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
    
    dysgu_files.multiMap{file ->
        pacbio: file
        ilumina: file
    }.set {figure_sizes}

    pacbio_data = generate_len_csv_from_vcf(figure_sizes.pacbio)
    ilumina_data = generate_len_csv_from_previous_tsv(figure_sizes.ilumina, file(params.previous_dir))

    figures  = calculate_size_distribution(pacbio_data, ilumina_data, file(params.size_distribution_script))

    figures.view()

}