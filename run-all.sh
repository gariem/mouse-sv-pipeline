DATA_TAG=$1 #Either 'prod' or 'test'
PARAMS_FILE="run-params-$DATA_TAG.txt"
NXF_ENV=$2 #Empty or '-q'
QUIET=$3 #Empty or '-profile lsf'

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
echo "NXF_ENV=$NXF_ENV"
echo "QUIET=$QUIET"
echo ""
echo "============================\n"
echo ""

echo "Cleaning work, merge, and reports folders"
rm -rf ./work
rm -rf ./data/reports
rm -rf ./data/merged 
echo "\tDone -> Clean\n"

echo "Collecting metrics From source VCFs"
nextflow $QUIET run tasks/metrics_collect.nf \
    --intersect_window $INTERSECT_WINDOW_VALUES \
    --filter_hets $FILTER_VALUES \
    --input_dir ./data/input $NXF_ENV

rm -rf ./work
echo "  Done -> Source Metrics\n"

echo "Merging with CombiSV"
    nextflow $QUIET run tasks/merge_combi.nf \
        --min_coverage $COMBI_MIN_COVERAGE_VALUES $NXF_ENV

echo "  Done -> Merge [CombiSV]\n"

echo "Collecting metrics for CombiSV"
    nextflow run tasks/metrics_collect.nf \
        --intersect_window $INTERSECT_WINDOW_VALUES \
        --filter_hets $FILTER_VALUES \
        --input_dir ./data/merged $NXF_ENV

rm -rf ./work
echo "  Done -> Metrics [CombiSV]\n"

echo "Merging with mapped method"
nextflow $QUIET run tasks/merge_survivor_mapped.nf \
    --max_dist $SURVIVOR_MAX_DIST_VALUES \
    --min_callers $SURVIVOR_MIN_CALLERS_VALUES \
    --min_size $SURVIVOR_MIN_SIZE_VALUES \
    --filter_hets $FILTER_VALUES $NXF_ENV

echo "  Done\n"

echo "Merging with survivor"
    nextflow run tasks/merge_survivor_simple.nf \
        --max_dist $SURVIVOR_MAX_DIST_VALUES \
        --min_callers $SURVIVOR_MIN_CALLERS_VALUES \
        --min_size $SURVIVOR_MIN_SIZE_VALUES $NXF_ENV

echo "  Done"

for STRAIN in $(echo $STRAIN_VALUES | sed "s/,/ /g")
do
    echo "Collecting metrics for merged $STRAIN [Mapped]"
    nextflow $QUIET run tasks/metrics_collect.nf \
        --intersect_window $INTERSECT_WINDOW_VALUES \
        --filter_hets 0 \
        --input_dir ./data/merged \
        --strain "$STRAIN-mapped" $NXF_ENV

    rm -rf ./work
    echo "  Done ($STRAIN [Mapped])\n"
    

    echo "Collecting metrics for merged $STRAIN [Survivor]"
    nextflow $QUIET run tasks/metrics_collect.nf \
        --intersect_window $INTERSECT_WINDOW_VALUES \
        --filter_hets $FILTER_VALUES \
        --input_dir ./data/merged \
        --strain "$STRAIN-survivor" $NXF_ENV

    rm -rf ./work
    echo "  Done ($STRAIN [Survivor])\n"
    
done

nextflow run tasks/metrics_consolidate.nf