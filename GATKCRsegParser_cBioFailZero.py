#!/usr/bin/env python

# Gather Seg reports from a directory and assemble them into 1 large table


import pandas as pd
import glob, os
import numpy as np

# build arg parser


# glob in files
report_list = glob.glob('*.anno.seg')

samples = []
# loop through kraken2 reports
for report in report_list:
    
    # get sample ID from filename
    sample_name = report.split('.')[0]
    print(sample_name)
    # read in as pandas dataframe
    result = pd.read_csv(report, sep="\t")
    #stripped failed lines - while a good idea this leaves Gaps that cuase Gistic to fail.
    # the hash above led to disjointed seg files and were rejected by GISTIC
    #result = result[result['Sample'].str.contains('Pass')]
    # drop undesired columns and strip white space
    result = result[["Sample", "Chr", "Start", "End", "CR_NumPoints", "GATK_LgTumorGeoMean"]]
    result[['Sample', 'STATUS']] = result['Sample'].str.split('_', 1, expand=True)
    
    print(result)
    #change Sample to reflect sample_name
    result['Sample'] = sample_name
    #rename columns to run Gistic
    result.rename(columns={"Chr": "Chromosome", "Start": "Start Position", "End": "End Position", "CR_NumPoints": "Num markers", "GATK_LgTumorGeoMean": "Seg.CN"}, inplace=True)
    result.loc[result.STATUS == 'Fail', ['Seg.CN']] = 0
    result.drop(['STATUS'], axis=1, inplace = True)
    #print(result)
    #save to file
    csv_file = os.path.splitext(report)[-2]+".Gistic.seg"
    result.to_csv(csv_file, sep='\t', index=False, header=None)
    
    
filenames = glob.glob('*.Gistic.seg')

df = (pd.read_table(f, header=None) for f in filenames)
#df.drop([0])
#print(df)
concatenated_df = pd.concat(df, ignore_index=True)

concatenated_df.to_csv('UroMasterCNV_NIMwFailZero4GISTIC.seg', sep='\t', index=False, header=None)
concatenated_df.columns =['ID', 'chrom', 'loc.start', 'loc.end', 'num.mark', 'seg.mean']
concatenated_df.to_csv('UroMasterCNV_NIMGATKCRsegParser_cBiowFailZero4seg.seg', sep='\t', index=False)


