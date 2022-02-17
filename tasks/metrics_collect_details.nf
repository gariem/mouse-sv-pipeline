#!/usr/bin/env nextflow

nextflow.enable.dsl = 2


params.file= "C57BL_6NJ-dysgu.vcf"
params.input_dir = "./data/input"
params.input_files = "${params.input_dir}/${params.file}"
params.out_dir = './data/reports/details'
params.previous_dir = './data/previous'
params.validated_dir = './data/validated'

params.merge_mappings = "./data/input/merge_mappings.txt"
params.metrics_mappings_file = "./data/input/metrics_mappings.txt"
params.genome_sizes = "./data/input/genome_size.txt"

params.all_sv_types='INS,DEL,INV,DUP'
params.compare_sv_types='INS,DEL'

params.reference_path='/home/egarcia/data/mouse/reference/Mus_musculus.GRCm39.dna.toplevel.fa'
params.new_alignment='/home/egarcia/data/mouse/minimap2/C57BL_6NJ.sorted.bam' 
params.old_alignment='/media/egarcia/DataBank/mouse/aligments_old/C57BL_6NJ.bam'


params.intersect_window = '30'

params.size_distribution_script='./templates/size_distribution.py'


process generate_new_data_file {

    input:
        file vcf_file

    output:
        file '*.new_data.csv'

    script:
    file_id = vcf_file.name.toString().tokenize('.').get(0)

    """
    bcftools query -f'%CHROM,%SVTYPE,%SVLEN\\n' ${vcf_file} | awk -F',' 'BEGIN {OFS = FS} \$1 ~/^[0-9]*\$|^X\$/{abs=\$3;sub("^-", "",abs); print \$1,\$2,abs}' > "${file_id}.new_data.csv"
    """
}

process generate_previous_data_file {

    input:
        file vcf_file
        file previous_dir

    output:
        file '*.previous_data.csv'

    script:

    file_id = vcf_file.name.toString().tokenize('.').get(0)
    strain = file_id.tokenize('-').get(0)
    
    """
    for file in ${previous_dir}/${strain}*; do
        TYPE="\$(echo \$file | cut -d '.' -f 2)"
        cat \$file | awk -v type=\$TYPE '{print \$1","type","\$4}' >> "${file_id}.previous_data.csv"
    done

    """
}

process calculate_size_distribution {

    publishDir file(params.out_dir), mode: "copy",  saveAs: {filename -> filename.replace(".sizes.png","")+'/size_distribution.png' }

    input:
        tuple val(file_key), file(data_files)
        file size_distribution_script
    
    output:
        file "${file_key}.sizes.png"

    script:
    
    if(data_files.get(0).getName().contains("new_data")){
        new_data = data_files.get(0)
        previous_data = data_files.get(1)
    }else{
        new_data = data_files.get(1)
        previous_data = data_files.get(0)
    }

    """
    python ${size_distribution_script} -i1 ${new_data} -i2 ${previous_data} -o "${file_key}.sizes.png"
    """
}

process generate_bed_files {

    publishDir file(params.out_dir), mode: "copy",  saveAs: {
                            filename -> filename.tokenize('.').get(0) + '/' + filename.replace(filename.tokenize('.').get(0)+ '.', "")
                            }

    input: 
        file vcf_file
        each type
    
    output:
        file '*.bed'

    script:
    
    file_id = vcf_file.name.toString().tokenize('.').get(0)
    strain = file_id.tokenize('-').get(0)

    """
    bcftools query -i"SVTYPE='${type}'" -f'%CHROM\\t%POS0\\t%END0\\t%SVLEN\\n' ${vcf_file} | \
            awk -F'\\t' 'BEGIN {OFS = FS} \$1 ~/^[0-9]*\$|^X\$/{print \$1,\$2,\$3,\$4}' > "${file_id}.${strain}.${type}.bed"
    """
}

process intersect_out_with_validated { 

    input: 
        file bed_file
        file validated_dir
        file metrics_mappings_file
        val window

    output:
        file '*.bed'

    script:

    file_id = bed_file.name.toString().tokenize('.').get(0)
    strain = file_id.tokenize('-').get(0)

    type = bed_file.name.tokenize('.').get(2)

    """
    echo "\$(cat ${metrics_mappings_file} | grep ${type} || true;)" > compare_to.txt

    awk -F'\\t' 'BEGIN {OFS = FS} {print \$1,\$2-${window},\$3+${window},\$4}' ${bed_file} | awk '!x[\$0]++' > FILE_B

    while read -r compare_line
    do
        VALIDATED_TYPE="\$(echo \$compare_line | cut -d '=' -f 1)"
        VALIDATED_FILE=${validated_dir}/${strain}.\$VALIDATED_TYPE.mm39.bed

        awk -F'\\t' 'BEGIN {OFS = FS} {print \$1,\$2-${window},\$3+${window},\$4}' \${VALIDATED_FILE} | awk '!x[\$0]++' > FILE_A

        bedtools intersect -a FILE_A -b FILE_B -v > "${file_id}.${strain}.${type}_\${VALIDATED_TYPE}.out.bed"
    done < compare_to.txt
    """
}

process take_screenshots { 

    publishDir file(params.out_dir), mode: "copy",  saveAs: {
                            filename -> filename.tokenize('.').get(1) + '/' + filename.replace(filename.tokenize('.').get(1)+ '.', "")
                            }

    input:
        file bed_file
        file reference
        file reference_index
        file new_alignment
        file new_alignment_index
        file old_alignment
        file new_alignment_index

    output:
        file 'snapshots/*.png'

    script:

    file_id = bed_file.name.toString().tokenize('.').get(0)

    """
    echo "new" > snapshots.txt
    echo "genome ${reference}" >> snapshots.txt
    echo "snapshotDirectory ./snapshots" >> snapshots.txt
    echo "load ${new_alignment}" >> snapshots.txt
    echo "load ${old_alignment}" >> snapshots.txt
    echo "load ${bed_file}" >> snapshots.txt
    
    bedToIgv -path ./snapshots -slop 50 -i ${bed_file} >> snapshots.txt
    bedToIgv -path ./snapshots -slop 200 -clps -i ${bed_file} >> snapshots.txt

    echo "exit" >> snapshots.txt

    sed -i -e 's/.png/.${file_id}.png/g' snapshots.txt

    echo "IGV.Bounds=0,0,1440,810" > prefs.properties
    echo "SAM.SHOW_ALL_BASES=false" >> prefs.properties
    echo "SAM.SHOW_SOFT_CLIPPED=true" >> prefs.properties
    echo "SAM.SHOW_JUNCTION_TRACK=false" >> prefs.properties
    echo "SAM.SHOW_JUNCTION_FLANKINGREGIONS=false" >> prefs.properties
    echo "DETAILS_BEHAVIOR=CLICK" >> prefs.properties

    xvfb-run --auto-servernum -s "-screen 0 1440x810x24" java -Xmx4000m --module-path=/IGV_Linux_2.12.2/lib --module=org.igv/org.broad.igv.ui.Main -b snapshots.txt -o prefs.properties

    """

}

workflow {

    Channel.fromPath(params.input_files).multiMap{file ->
        new_data: file
        previous_data: file
        final_data: file
        comparative_data: file
    }.set {input}

    new_data = generate_new_data_file(input.new_data)
    previous_data = generate_previous_data_file(input.previous_data, file(params.previous_dir))

    new_data.concat(previous_data).map{ file -> 
        def key = file.name.toString().tokenize('.').get(0)
        return tuple(key,file)
    }
    .groupTuple(size: 2)
    .set{ grouped_data }

    calculate_size_distribution(grouped_data, file(params.size_distribution_script))

    all_sv_types = Channel.value(params.all_sv_types.split(",")).flatten()
    compare_sv_types = Channel.value(params.compare_sv_types.split(",")).flatten()

    final_beds = generate_bed_files(input.final_data, all_sv_types).filter { 
                    params.compare_sv_types.split(',').contains(it.name.toString().tokenize('.').get(2)) 
                }
    
    out_features = intersect_out_with_validated(final_beds, file(params.validated_dir), file(params.metrics_mappings_file), params.intersect_window).flatten()

    // take_screenshots(out_features, file(params.reference_path), file(params.reference_path + '.fai'),
    //                     file(params.new_alignment), file(params.new_alignment + '.bai'), 
    //                     file(params.old_alignment), file(params.old_alignment + '.bai')
    //                 )


    
}