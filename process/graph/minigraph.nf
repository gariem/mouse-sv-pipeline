#!/usr/bin/env nextflow

nextflow.enable.dsl = 2


process bed_from_full_graph {
    cache false

    input: 
        val strain
        file folder
    
    output:
        file "*.bed"

    script:
    
    simple_name = strain + '-minigraph'

    """
    awk -F"[\t:]" 'BEGIN {OFS = "\t"} {if(\$6!="."&&(\$3-\$2)<\$7)print \$1,\$2,\$3,\$3-\$2}' ${folder}/${strain}.bed > "${strain}-minigraph__INS.bed"
    awk -F"[\t:]" 'BEGIN {OFS = "\t"} {if(\$6!="."&&(\$3-\$2)>\$7)print \$1,\$2,\$3,\$3-\$2}' ${folder}/${strain}.bed > "${strain}-minigraph__DEL.bed"
    """

}


process bed_from_simple_graph {
    cache false

    input: 
        val strain
        file folder
    
    output:
        file "*.bed"

    script:
    
    simple_name = strain + '-minigraph'

    """
    awk -F"[\t:]" 'BEGIN {OFS = "\t"} {if(\$6!="."&&(\$3-\$2)<\$7)print \$1,\$2,\$3,\$3-\$2}' ${folder}/${strain}.bed > "${strain}-minigraph__INS.bed"
    awk -F"[\t:]" 'BEGIN {OFS = "\t"} {if(\$6!="."&&(\$3-\$2)>\$7)print \$1,\$2,\$3,\$3-\$2}' ${folder}/${strain}.bed > "${strain}-minigraph__DEL.bed"
    """

}