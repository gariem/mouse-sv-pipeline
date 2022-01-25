DATA_TAG=$1 #Either 'prod' or 'test'
PARAMS_FILE="run-params-$DATA_TAG.txt"
PROFILE=$2 #Empty or '-profile lsf'
QUIET=$3 #Empty or '-q'

STRAIN_VALUES="$(cat $PARAMS_FILE | grep 'strains=' | cut -d '=' -f 2 )"
SURVIVOR_MAX_DIST_VALUES="$(cat $PARAMS_FILE | grep 'survivor.max_dist=' | cut -d '=' -f 2 )"
SURVIVOR_MIN_SIZE_VALUES="$(cat $PARAMS_FILE | grep 'survivor.min_size=' | cut -d '=' -f 2 )"
SURVIVOR_MIN_CALLERS_VALUES="$(cat $PARAMS_FILE | grep 'survivor.min_callers=' | cut -d '=' -f 2 )"
COMBI_MIN_COVERAGE_VALUES="$(cat $PARAMS_FILE | grep 'combi.min_coverage=' | cut -d '=' -f 2 )"
INTERSECT_WINDOW_VALUES="$(cat $PARAMS_FILE | grep 'intersect.window=' | cut -d '=' -f 2 )"
FILTER_VALUES="$(cat $PARAMS_FILE | grep 'filter_hets=' | cut -d '=' -f 2 )"

echo "========== PARAMS =========="
cat $PARAMS_FILE
echo ""
echo "DATA_TAG=$DATA_TAG"
echo "PROFILE=$PROFILE"
echo "QUIET=$QUIET"
echo ""
echo "============================"
echo ""

echo "Cleaning work, merge, and reports folders"
rm -rf ./work
rm -rf ./data/reports
rm -rf ./data/merged 
echo "  Done -> Clean"

echo "Collecting metrics From source VCFs"
nextflow $QUIET run tasks/metrics_collect_general.nf \
    --intersect_window $INTERSECT_WINDOW_VALUES \
    --filter_hets $FILTER_VALUES \
    --input_dir ./data/input $PROFILE

rm -rf ./work
echo "  Done -> Source Metrics\n"

echo "Merging with CombiSV"
    nextflow $QUIET run tasks/merge_combi.nf \
        --min_coverage $COMBI_MIN_COVERAGE_VALUES $PROFILE

echo "  Done -> Merge [CombiSV]\n"

echo "Collecting metrics for CombiSV"
    nextflow run tasks/metrics_collect_general.nf \
        --intersect_window $INTERSECT_WINDOW_VALUES \
        --filter_hets $FILTER_VALUES \
        --input_dir ./data/merged $PROFILE

rm -rf ./work
echo "  Done -> Metrics [CombiSV]\n"

echo "Merging with mapped method"
nextflow $QUIET run tasks/merge_survivor_mapped.nf \
    --max_dist $SURVIVOR_MAX_DIST_VALUES \
    --min_callers $SURVIVOR_MIN_CALLERS_VALUES \
    --min_size $SURVIVOR_MIN_SIZE_VALUES \
    --filter_hets $FILTER_VALUES $PROFILE

echo "  Done"

echo "Merging with survivor"
    nextflow run tasks/merge_survivor_simple.nf \
        --max_dist $SURVIVOR_MAX_DIST_VALUES \
        --min_callers $SURVIVOR_MIN_CALLERS_VALUES \
        --min_size $SURVIVOR_MIN_SIZE_VALUES $PROFILE

echo "  Done"

for STRAIN in $(echo $STRAIN_VALUES | sed "s/,/ /g")
do
    echo "Collecting metrics for merged $STRAIN [Mapped]"
    nextflow $QUIET run tasks/metrics_collect_general.nf \
        --intersect_window $INTERSECT_WINDOW_VALUES \
        --filter_hets 0 \
        --input_dir ./data/merged \
        --strain "$STRAIN-mapped" $PROFILE

    rm -rf ./work
    echo "  Done ($STRAIN [Mapped])"
    

    echo "Collecting metrics for merged $STRAIN [Survivor]"
    nextflow $QUIET run tasks/metrics_collect_general.nf \
        --intersect_window $INTERSECT_WINDOW_VALUES \
        --filter_hets $FILTER_VALUES \
        --input_dir ./data/merged \
        --strain "$STRAIN-survivor" $PROFILE

    rm -rf ./work
    echo "  Done ($STRAIN [Survivor])"
    
done

nextflow $QUIET run tasks/metrics_consolidate.nf $PROFILE