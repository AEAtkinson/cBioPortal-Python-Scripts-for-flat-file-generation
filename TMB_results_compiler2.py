#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan  7 14:14:58 2020

@author: aaronatkinson
"""

#!/usr/bin/env python

# USAGE: collect_kraken2_reports.py --dir [directory to look for reports in] --out [results file name] 

import glob
import pandas as pd


#dir = '/Users/aaronatkinson/Documents/Foundation/FoundationReports/MSIresults/MantisAvatarResults'

# glob in files
report_list = glob.glob('*_Tmb.txt')

dfs = []
samples = []
# loop through MANTIS reports
for report in report_list:	
# get sample ID from filename
    sample_name = report
    print(sample_name)
	# read in as pandas dataframe
    sample = pd.read_csv(report, sep="\t", skiprows = 2, header = None)
    sample.insert(0, column='SAMPLE_ID', value=sample_name)
    sample.rename(columns={0: 'TMB'}, inplace=True)
    sample['TMB'] = sample['TMB'].str.replace('tmb ', '')
    print(sample)
    dfs.append(sample)
Concat_table = pd.concat(dfs)
Concat_table = Concat_table.apply(lambda x: x.str.replace('_Illumina_Hg38_Tmb.txt',''))
Concat_table.to_csv('Concat_table.txt', sep='\t', index=False)
