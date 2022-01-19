#!/usr/bin/env nextflow

params.strain = 'DBA_2J'
params.raw_metrics = "./data/reports/raw/*.csv"
params.out_dir = "./data/reports"

Channel.fromPath(params.raw_metrics).set{files}

process create_file_list {
    echo true

    input:
        file csv_files from files.collect()

    output:
        file 'list.txt' into list_file
        file '*.csv' into csvs_files

    """
    for FILE in ${csv_files}
    do
        echo "cp_\$FILE" >> list.txt
        mv \$FILE "cp_\$FILE"
    done
    """
}

process process_list {

    publishDir file(params.out_dir), mode: 'copy'

    input:
        file list from list_file
        file csvs from csvs_files
        val strain from params.strain
    output:
        file "report_data.csv"
    
    """
#!python

import pandas as pd

data = pd.DataFrame()

with open('${list}', "r") as list:
    for file in list:
        file_data = pd.read_csv(file.replace("\\n", ""), low_memory=False)
        data = data.append(file_data)
data['SCORE'] = (data['H1_INT'] + data['H2_INT'] + data['H6_INT'] + data['H7_INT']) / (data['H1_TOT'] + data['H2_TOT'] + data['H6_TOT'] + data['H7_TOT'])
data['DIFF'] = (data['INS0'] + data['DEL0'] + data['INV0'] + data['DUP0']) / (data['INS'] + data['DEL'] + data['INV'] + data['DUP'])
data.to_csv("report_data.csv", index=False, header=True)


    """

}