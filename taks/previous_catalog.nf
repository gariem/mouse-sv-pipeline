#!/usr/bin/env nextflow

params.input_catalog = './data/input/18strains.REL-1302-SV-GRCm38.sdp.tab.gz'
params.strain_defs = './data/input/strain_defs.txt'


process uncompress {

    input:
        file compressed_catalog from file(params.input_catalog)
    
    output:
        file '18strains.REL-1302-SV-GRCm38.sdp.tab' into catalog_data

    """
    gzip -d --force ${compressed_catalog}
    """
}


process split_columns {
    
    input:
        file data from catalog_data
        file strain_defs from file(params.strain_defs)
    
    output:
        file '*.txt' into strain_data

    """
    grep "#CHROM" ${data} |  tr "\\t" "\\n" > columns

    COL=1

    while read -r line
    do
        
        STRAIN_DEF="\$(grep \$line ${strain_defs} | head -1 || true )"

        if [ ! -z "\$STRAIN_DEF" ]
        then
            SRAIN_NAME="\$(echo \$STRAIN_DEF | cut -d '=' -f 2)"
            cut -f \$COL < ${data} > \$SRAIN_NAME.txt
        fi

        COL=\$(( \$COL + 1 ))
    done < columns
    """

}


process parse_strain_data {

    input: 
        file data from strain_data.flatten()
    
    output:
        file '*.bed' into bed_files

    script:

    strain = data.getName().tokenize(".").get(0)

    """
    grep ":" $data | grep ";" | grep -v "#" | \
        awk -F'\\t' 'BEGIN {OFS = FS} {split(\$1,arr,";"); split(arr[1],pos,":"); split(pos[2],bps,"-"); \
        len=bps[2]-bps[1]; split(arr[2],type,"|"); print pos[1],bps[1],bps[2],len,type[1]}' > ${strain}.bed
    """
}

bed_files.view()