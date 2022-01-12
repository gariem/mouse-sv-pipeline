PARAMS_FILE="run-params$1.txt"
STRAIN_VALUES="$(cat $PARAMS_FILE | grep 'strains=' | cut -d '=' -f 2 )"
SURVIVOR_MAX_DIST_VALUES="$(cat $PARAMS_FILE | grep 'survivor.max_dist=' | cut -d '=' -f 2 )"
SURVIVOR_MIN_SIZE_VALUES="$(cat $PARAMS_FILE | grep 'survivor.min_size=' | cut -d '=' -f 2 )"
SURVIVOR_MIN_CALLERS_VALUES="$(cat $PARAMS_FILE | grep 'survivor.min_callers=' | cut -d '=' -f 2 )"
COMBI_MIN_COVERAGE_VALUES="$(cat $PARAMS_FILE | grep 'combi.min_coverage=' | cut -d '=' -f 2 )"
INTERSECT_WINDOW_VALUES="$(cat $PARAMS_FILE | grep 'intersect.window=' | cut -d '=' -f 2 )"
FILTER_VALUES="$(cat $PARAMS_FILE | grep 'filter_hets=' | cut -d '=' -f 2 )"

mkdir -p ./work/commands

for STRAIN in $(echo $STRAIN_VALUES | sed "s/,/ /g")
do
    CMD_FILE=./work/commands/$STRAIN.sh
    echo "#Commands for Strain $STRAIN" > $CMD_FILE

    echo "Processing strain: $STRAIN"
    for MAX_DIST in $(echo $SURVIVOR_MAX_DIST_VALUES | sed "s/,/ /g")
    do 
        for MIN_SIZE in $(echo $SURVIVOR_MIN_SIZE_VALUES | sed "s/,/ /g")
        do 
            for MIN_CALLERS in $(echo $SURVIVOR_MIN_CALLERS_VALUES | sed "s/,/ /g")
            do
                for FILTER in $(echo $FILTER_VALUES | sed "s/,/ /g")
                do 
                    echo "nextflow run tasks/merge_survivor_mapped.nf --strain $STRAIN --max_dist $MAX_DIST --min_callers $MIN_CALLERS --min_size $MIN_SIZE --filter_hets $FILTER" >> $CMD_FILE
                    echo "nextflow run tasks/merge_survivor_simple.nf --strain $STRAIN --max_dist $MAX_DIST --min_callers $MIN_CALLERS --min_size $MIN_SIZE --filter_hets $FILTER" >> $CMD_FILE
                done
            done
        done
    done

    for MIN_COVERAGE in $(echo $COMBI_MIN_COVERAGE_VALUES | sed "s/,/ /g")
    do 
        echo "nextflow run tasks/merge_combi.nf --strain $STRAIN --min_coverage $MIN_COVERAGE" >> $CMD_FILE
    done

    for WINDOW in $(echo $INTERSECT_WINDOW_VALUES | sed "s/,/ /g")
    do
        for FILTER in $(echo $FILTER_VALUES | sed "s/,/ /g")
        do 
            echo "nextflow run tasks/metrics_collect.nf --intersect_window_a $WINDOW --intersect_window_b $WINDOW --filter_hets $FILTER --input_files ./data/input/$STRAIN-sniffles.vcf" >> $CMD_FILE
            echo "nextflow run tasks/metrics_collect.nf --intersect_window_a $WINDOW --intersect_window_b $WINDOW --filter_hets $FILTER --input_files ./data/input/$STRAIN-pbsv.vcf" >> $CMD_FILE
            echo "nextflow run tasks/metrics_collect.nf --intersect_window_a $WINDOW --intersect_window_b $WINDOW --filter_hets $FILTER --input_files ./data/input/$STRAIN-cutesv.vcf" >> $CMD_FILE

            for MIN_COVERAGE in $(echo $COMBI_MIN_COVERAGE_VALUES | sed "s/,/ /g")
            do
                echo "nextflow run tasks/metrics_collect.nf --intersect_window_a $WINDOW --intersect_window_b $WINDOW --filter_hets $FILTER --input_files ./data/merged/$STRAIN-combisv_c$MIN_COVERAGE.vcf" >> $CMD_FILE
            done
        done

        for MAX_DIST in $(echo $SURVIVOR_MAX_DIST_VALUES | sed "s/,/ /g")
        do 
            for MIN_SIZE in $(echo $SURVIVOR_MIN_SIZE_VALUES | sed "s/,/ /g")
            do 
                for MIN_CALLERS in $(echo $SURVIVOR_MIN_CALLERS_VALUES | sed "s/,/ /g")
                do
                    echo "nextflow run tasks/metrics_collect.nf --intersect_window_a $WINDOW --intersect_window_b $WINDOW --input_files ./data/merged/$STRAIN-survivor_${MAX_DIST}_${MIN_CALLERS}_${MIN_SIZE}.vcf" >> $CMD_FILE
                    echo "nextflow run tasks/metrics_collect.nf --intersect_window_a $WINDOW --intersect_window_b $WINDOW --input_files ./data/merged/$STRAIN-survivor_${MAX_DIST}_${MIN_CALLERS}_${MIN_SIZE}_nohets.vcf" >> $CMD_FILE

                    for FILTER in $(echo $FILTER_VALUES | sed "s/,/ /g")
                    do 
                        echo "nextflow run tasks/metrics_collect.nf --intersect_window_a $WINDOW --intersect_window_b $WINDOW --filter_hets $FILTER --input_files ./data/merged/$STRAIN-unmapped_${MAX_DIST}_${MIN_CALLERS}_${MIN_SIZE}.vcf" >> $CMD_FILE
                    done
                done
            done
        done
        
    done
    echo "echo \"$STRAIN Finished\"" >> $CMD_FILE
    echo "Launchin script from: $CMD_FILE"
    #bash $CMD_FILE &
done

# echo "nextflow run tasks/metrics_consolidate.nf"