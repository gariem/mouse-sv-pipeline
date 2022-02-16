#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process take_screenshots { 

    input:
        file bed_file
        file reference
        file reference_index
        file new_alignment
        file new_alignment_index
        file old_alignment
        file new_alignment_index
        val outname

    output:
        file "${outname}"

    script:


    """
    echo "new" > snapshots.txt
    echo "genome ${reference}" >> snapshots.txt
    echo "snapshotDirectory ./snapshots" >> snapshots.txt
    echo "load ${new_alignment}" >> snapshots.txt
    echo "load ${old_alignment}" >> snapshots.txt
    echo "load ${bed_file}" >> snapshots.txt
    
    bedToIgv -path ./snapshots -slop 50 -i ${bed_file} >> snapshots.txt
    bedToIgv -path ./snapshots -slop 500 -clps -i ${bed_file} >> snapshots.txt

    echo "exit" >> snapshots.txt

    sed -i -e 's/.png/.${outname}.png/g' snapshots.txt

    echo "IGV.Bounds=0,0,1440,810" > prefs.properties
    echo "SAM.SHOW_ALL_BASES=false" >> prefs.properties
    echo "SAM.SHOW_SOFT_CLIPPED=true" >> prefs.properties
    echo "SAM.SHOW_JUNCTION_TRACK=false" >> prefs.properties
    echo "SAM.SHOW_JUNCTION_FLANKINGREGIONS=false" >> prefs.properties
    echo "DETAILS_BEHAVIOR=CLICK" >> prefs.properties

    xvfb-run --auto-servernum -s "-screen 0 1440x810x24" java -Xmx4000m --module-path=/IGV_Linux_2.12.2/lib --module=org.igv/org.broad.igv.ui.Main -b snapshots.txt -o prefs.properties

    """

}