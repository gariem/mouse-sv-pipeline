#!/usr/bin/env nextflow

params.filter = "*.csv"
params.input_dir = "./data/reports/raw"

params.test_raw_file = "${params.input_dir}/test.csv"

params.out_dir = "./data/reports"

inputDir = file(params.test_raw_file).getParent()

process consolidate_metrics {

    publishDir file(params.out_dir), mode: 'copy'
    
    input:
        val filter from params.filter

    output:
        file "report_data.csv"
    
    script:

    """
#!python

import pandas as pd
import glob
import os

all_files = glob.glob(os.path.join("${inputDir}", "${filter}"))

data = pd.DataFrame()

for file in all_files:
    file_data = pd.read_csv(file.replace("\\n", ""), low_memory=False)
    data = data.append(file_data)

data['INS%'] = (data['I_H6'] + data['I_H7']) / (data['T_H6'] + data['T_H7'])
data['DEL%'] = (data['I_H1'] + data['I_H2']) / (data['T_H1'] + data['T_H2'])
data['TOT%'] = (data['I_H1'] + data['I_H2'] + data['I_H6'] + data['I_H7']) / (data['T_H1'] + data['T_H2'] + data['T_H6'] + data['T_H7'])
data['DIFF'] = (data['INS'] + data['DEL'] + data['INV'] + data['DUP']) / (data['INS0'] + data['DEL0'] + data['INV0'] + data['DUP0']) 

data.to_csv("report_data.csv", index=False, header=True)
    """
}
