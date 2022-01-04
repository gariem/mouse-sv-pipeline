while read -r STRAIN
do

    nextflow run taks/merge_survivor_complex.nf --strain $STRAIN --max_dist 50 --min_callers 1 --min_size 20 --filter_hets 0
    nextflow run taks/merge_survivor_complex.nf --strain $STRAIN --max_dist 100 --min_callers 1 --min_size 20 --filter_hets 0
    nextflow run taks/merge_survivor_complex.nf --strain $STRAIN --max_dist 50 --min_callers 1 --min_size 20 --filter_hets 1
    nextflow run taks/merge_survivor_complex.nf --strain $STRAIN --max_dist 100 --min_callers 1 --min_size 20 --filter_hets 1

    nextflow run taks/metrics_collect.nf --intersect_window_a 20 --intersect_window_b 20 --filter_hets 0 --input_files ./data/input/$STRAIN-sniffles.vcf
    nextflow run taks/metrics_collect.nf --intersect_window_a 30 --intersect_window_b 30 --filter_hets 0 --input_files ./data/input/$STRAIN-sniffles.vcf
    nextflow run taks/metrics_collect.nf --intersect_window_a 50 --intersect_window_b 50 --filter_hets 0 --input_files ./data/input/$STRAIN-sniffles.vcf

    nextflow run taks/metrics_collect.nf --intersect_window_a 20 --intersect_window_b 20 --filter_hets 1 --input_files ./data/input/$STRAIN-sniffles.vcf
    nextflow run taks/metrics_collect.nf --intersect_window_a 30 --intersect_window_b 30 --filter_hets 1 --input_files ./data/input/$STRAIN-sniffles.vcf
    nextflow run taks/metrics_collect.nf --intersect_window_a 50 --intersect_window_b 50 --filter_hets 1 --input_files ./data/input/$STRAIN-sniffles.vcf

    nextflow run taks/metrics_collect.nf --intersect_window_a 20 --intersect_window_b 20 --filter_hets 0 --input_files ./data/input/$STRAIN-pbsv.vcf
    nextflow run taks/metrics_collect.nf --intersect_window_a 30 --intersect_window_b 30 --filter_hets 0 --input_files ./data/input/$STRAIN-pbsv.vcf
    nextflow run taks/metrics_collect.nf --intersect_window_a 50 --intersect_window_b 50 --filter_hets 0 --input_files ./data/input/$STRAIN-pbsv.vcf

    nextflow run taks/metrics_collect.nf --intersect_window_a 20 --intersect_window_b 20 --filter_hets 1 --input_files ./data/input/$STRAIN-pbsv.vcf
    nextflow run taks/metrics_collect.nf --intersect_window_a 30 --intersect_window_b 30 --filter_hets 1 --input_files ./data/input/$STRAIN-pbsv.vcf
    nextflow run taks/metrics_collect.nf --intersect_window_a 50 --intersect_window_b 50 --filter_hets 1 --input_files ./data/input/$STRAIN-pbsv.vcf

    nextflow run taks/metrics_collect.nf --intersect_window_a 20 --intersect_window_b 20 --input_files ./data/merged/$STRAIN-survivor_50_1_20.vcf
    nextflow run taks/metrics_collect.nf --intersect_window_a 30 --intersect_window_b 30 --input_files ./data/merged/$STRAIN-survivor_50_1_20.vcf
    nextflow run taks/metrics_collect.nf --intersect_window_a 50 --intersect_window_b 50 --input_files ./data/merged/$STRAIN-survivor_50_1_20.vcf

    nextflow run taks/metrics_collect.nf --intersect_window_a 20 --intersect_window_b 20 --input_files ./data/merged/$STRAIN-survivor_100_1_20.vcf
    nextflow run taks/metrics_collect.nf --intersect_window_a 30 --intersect_window_b 30 --input_files ./data/merged/$STRAIN-survivor_100_1_20.vcf
    nextflow run taks/metrics_collect.nf --intersect_window_a 50 --intersect_window_b 50 --input_files ./data/merged/$STRAIN-survivor_100_1_20.vcf

    nextflow run taks/metrics_collect.nf --intersect_window_a 20 --intersect_window_b 20 --input_files ./data/merged/$STRAIN-survivor_50_1_20_nohets.vcf
    nextflow run taks/metrics_collect.nf --intersect_window_a 30 --intersect_window_b 30 --input_files ./data/merged/$STRAIN-survivor_50_1_20_nohets.vcf
    nextflow run taks/metrics_collect.nf --intersect_window_a 50 --intersect_window_b 50 --input_files ./data/merged/$STRAIN-survivor_50_1_20_nohets.vcf

    nextflow run taks/metrics_collect.nf --intersect_window_a 20 --intersect_window_b 20 --input_files ./data/merged/$STRAIN-survivor_100_1_20_nohets.vcf
    nextflow run taks/metrics_collect.nf --intersect_window_a 30 --intersect_window_b 30 --input_files ./data/merged/$STRAIN-survivor_100_1_20_nohets.vcf
    nextflow run taks/metrics_collect.nf --intersect_window_a 50 --intersect_window_b 50 --input_files ./data/merged/$STRAIN-survivor_100_1_20_nohets.vcf
done < strains.txt

nextflow run taks/metrics_consolidate.nf