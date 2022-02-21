for file in data/input/calls/*.vcf; do
	OUT="$(echo $file | cut -d '/' -f 4 | cut -d '.' -f 1)"
	echo "$OUT"
	bcftools query -i"SVTYPE='INV' && SVLEN>1000000" -f'%CHROM\t%POS0\t%END0\t%SVLEN\n' $file | \
		awk -F'\t' 'BEGIN {OFS = FS} {print $1,$2,$3,$4}' > data/analysis/biginv/$OUT.bed

done
