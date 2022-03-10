#!/usr/bin/env nextflow

nextflow.enable.dsl = 2


process size_data_from_bedfiles {

    input:
        tuple val(simple_name), val(type), file(ilumina), file(pacbio)
        
    output:
        tuple val(simple_name), val(type), file("*.ilumina.csv"), file("*.pacbio.csv")

    script:

    """
    awk '{abs=\$4; sub("^-", "",abs); print \$1",${type},"abs}' ${ilumina} > "${simple_name}.${type}.ilumina.csv"
    awk '{abs=\$4; sub("^-", "",abs); print \$1",${type},"abs}' ${pacbio} > "${simple_name}.${type}.pacbio.csv"

    """
}

process filtered_size_data_from_bedfiles {

    input:
        tuple val(type), val(simple_name), file(ilumina), file(pacbio), file(repeated)
        
    output:
        tuple val(simple_name), val(type), file("*.ilumina.csv"), file("*.pacbio.csv")

    script:

    """
    bedtools intersect -a ${pacbio} -b ${repeated} -v > filtered_pacbio.bed
    awk '{abs=\$4; sub("^-", "",abs); print \$1",${type},"abs}' ${ilumina} > "${simple_name}.${type}.ilumina.csv"
    awk '{abs=\$4; sub("^-", "",abs); print \$1",${type},"abs}' filtered_pacbio.bed > "${simple_name}.${type}.pacbio.csv"

    """
}

process calculate_size_distribution {

    publishDir file(params.out_dir), mode: "copy",  saveAs: {filename -> filename.tokenize('-').get(0) + '/' + filename.replace(".sizes.png","")+'/figure_sizes.png' }
    
    input:
        tuple val(simple_name), file(ilumina_data), file(pacbio_data)
        file size_distribution_script
    
    output:
        file "*.png"

    script:

    strain = simple_name.tokenize('-').get(0)

    """

    cat ${pacbio_data} > pacbio.csv
    cat ${ilumina_data} > ilumina.csv

    python ${size_distribution_script} -s ${strain} -i1 pacbio.csv -i2 ilumina.csv -o "${simple_name}.sizes.png"
    rm pacbio.csv
    rm ilumina.csv
    """
}


process calculate_size_distribution_filtered {

    publishDir file(params.out_dir), mode: "copy",  saveAs: {filename -> filename.tokenize('-').get(0) + '/' + filename.replace(".sizes.png","")+'/figure_sizes_filtered.png' }
    
    input:
        tuple val(simple_name), file(ilumina_data), file(pacbio_data)
        file size_distribution_script
    
    output:
        file "*.png"

    script:

    strain = simple_name.tokenize('-').get(0)

    """

    cat ${pacbio_data} > pacbio.csv
    cat ${ilumina_data} > ilumina.csv

    python ${size_distribution_script} -s ${strain} -i1 pacbio.csv -i2 ilumina.csv -o "${simple_name}.sizes.png"
    rm pacbio.csv
    rm ilumina.csv
    """
}
