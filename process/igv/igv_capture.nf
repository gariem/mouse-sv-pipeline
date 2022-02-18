#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process take_screenshots { 
    
    //save as: {out_dir}/strain/simple_name/chr_pos_end_slopX.png
    publishDir file(params.out_dir), mode: "copy",  saveAs: {
                    filename -> filename.split('_-_')[1].tokenize('-').get(0) + '/' + filename.split('_-_')[1].replace(".png", "") + '/' + filename.replace(filename.split('_-_')[1], "").replace("_-_","") + ".png"
                }    

    input:
        file bed_file
        file igv_workdir

    output:
        file 'snapshots/*.png'

    script:

    simple_name = bed_file.name.split('__')[0]
    strain = simple_name.tokenize('-').get(0)

    """
    echo "new" > snapshots.txt
    echo "genome ${igv_workdir}/Mus_musculus.GRCm39.dna.toplevel.fa" >> snapshots.txt
    echo "snapshotDirectory ./snapshots" >> snapshots.txt
    echo "load ${igv_workdir}/${strain}.pacbio.bam" >> snapshots.txt
    echo "load ${igv_workdir}/${strain}.ilumina.bam" >> snapshots.txt
    echo "load ${bed_file}" >> snapshots.txt
    
    bedToIgv -slop 50 -i ${bed_file} | grep -v snapshotDirectory >> snapshots.txt
    bedToIgv -slop 500 -clps -i ${bed_file} | grep -v snapshotDirectory >> snapshots.txt

    echo "exit" >> snapshots.txt

    sed -i -e 's/.png/_-_${simple_name}.png/g' snapshots.txt

    echo "IGV.Bounds=0,0,1440,810" > prefs.properties
    echo "SAM.SHOW_ALL_BASES=false" >> prefs.properties
    echo "SAM.SHOW_SOFT_CLIPPED=true" >> prefs.properties
    echo "SAM.SHOW_JUNCTION_TRACK=false" >> prefs.properties
    echo "SAM.SHOW_JUNCTION_FLANKINGREGIONS=false" >> prefs.properties
    echo "DETAILS_BEHAVIOR=CLICK" >> prefs.properties

    echo "processing ${bed_file}"
    xvfb-run --auto-servernum -s "-screen 0 1024x768x24" java -Xmx4000m --module-path=/IGV_Linux_2.12.2/lib --module=org.igv/org.broad.igv.ui.Main -b snapshots.txt -o prefs.properties


    """

}