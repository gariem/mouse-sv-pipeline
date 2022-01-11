#!/usr/bin/env nextflow

params.issues_url = "https://ftp.ncbi.nlm.nih.gov/pub/grc/mouse/GRC/Issue_Mapping/GRCm39_issues.gff3"
params.output_bed = "./data/issues/GRCm39_issues.bed"

final_name = file(params.output_bed).getName()
final_dir = file(params.output_bed).getParent()

process download_issues {
    
    output:
        file "GRCm39_issues.gff3" into gff_issues

    """
    wget -qO- ${params.issues_url} > GRCm39_issues.gff3
    """
}

process issues_to_bed {
    
    publishDir final_dir, mode: "copy"

    input:
        file "GRCm39_issues.gff3" from gff_issues
    
    output:
        file "${final_name}"

    """
    awk  -F'\t' 'BEGIN {OFS = FS} \
    {split(\$9,info,";"); split(info[3],chr,"="); split(info[1],key,"="); \
    split(info[2],issue,"="); if (chr[2]!="Un" && \$4!="") print chr[2],\$4,\$5,issue[2]" ("key[2]")"}' "GRCm39_issues.gff3" > ${final_name}
    """

}