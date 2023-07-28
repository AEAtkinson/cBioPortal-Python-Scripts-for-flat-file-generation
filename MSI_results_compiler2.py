#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan  7 14:14:58 2020

@author: aaronatkinson
"""

#!/usr/bin/env python

# USAGE: collect_kraken2_reports.py --dir [directory to look for reports in] --out [results file name] 

import argparse, glob, os, pandas


# build arg parser
parser = argparse.ArgumentParser(description = 'parse and assemble MANTIS reports')
#dir = '/Users/aaronatkinson/Documents/Foundation/FoundationReports/MSIresults/MantisAvatarResults'
parser.add_argument('--out', help = 'output filename')
args = parser.parse_args()

# glob in files
report_list = glob.glob('*.txt')

samples = []
# loop through MANTIS reports
for report in report_list:	
# get sample ID from filename
    sample_name = report.split('.')[-2]
    print(sample_name)
	# read in as pandas dataframe
    sample = pandas.read_csv(report, sep="\t")
    sample.columns = ["Average Metric Value (Abbr)", "Value", "Threshold", "Status"]
    # drop undesired columns
    sample = sample[["Threshold", "Status"]]	
    sample.loc[[0], "Threshold"] = [sample_name]
    sample.rename(columns={"Threshold": "HCI"}, inplace =True)
    sample = sample[:1]
    print(sample)
    csv_file = os.path.splitext(sample_name)[0]+".temp.tsv"
    sample.to_csv(csv_file, sep='\t', index=False)

filenames = glob.glob('*.temp.tsv')
dfs = []
for filename in filenames:
    dfs.append(pandas.read_table(filename))
results = pandas.concat(dfs, ignore_index=True)
results.to_csv('AvatarMsiResults.tsv', sep='\t', index=False)
