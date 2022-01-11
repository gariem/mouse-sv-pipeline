# Mouse SV Pipelines


### Merge SVs from VCF Files using SURVIVOR-MAPPED pipeline:

```
nextflow run -with-docker raphsoft/mousesv:1.0 \
  tasks/merge_survivor_mapped.nf \
  --strain DBA_2J \
  --max_dist 200 --min_callers 1 \
  --min_size 20 --filter_hets 1
```

### Merge SVs from VCF Files using SURVIVOR-MAPPED pipeline:
```
nextflow run -with-docker raphsoft/mousesv:1.0 \
  tasks/metrics_collect.nf \
  --intersect_window_a 30 --intersect_window_b 30 \
  --filter_hets 0 
  --input_files ./data/merged/DBA_2J-survivor_200_1_20_nohets.vcf
```