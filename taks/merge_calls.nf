#!/usr/bin/env nextflow

params.strain = "DBA_2J"
params.caller = "pbsv"
params.calls = "./data/input/DBA_2J-pbsv.vcf"
params.mapping = "./data/input/mappings.txt"
params.output_bed = "./data/merged/file.bed"

final_name = file(params.output_bed).getName()
final_dir = file(params.output_bed).getParent()


sv_types_ch = Channel.of("DEL","INS")
mappings_ch = Channel.of('SVTYPE="DEL" || SVTYPE="DEL/INV"', 'SVTYPE="INS"')

calls = file(params.calls)

process vcf_to_bed {

    publishDir final_dir, mode: "copy"

    input:
        val type from sv_types_ch
        val mapping from mappings_ch
        file 'DBA_2J-pbsv.vcf' from calls
    
    output:
        file "${params.strain}-${params.caller}-${type}.bed" into file
   

    """
    /home/egarcia/appdir/bcftools/bin/bcftools query -i '${mapping}' -f'%CHROM\t%POS0\t%END0\t%SVLEN\n' DBA_2J-pbsv.vcf | \
    awk -F'\t' 'BEGIN {OFS = FS} \$1 ~/^[0-9]*\$|^X\$/{print \$1,\$2,\$3,\$4}' >> "${params.strain}-${params.caller}-${type}.bed"
    """
}
