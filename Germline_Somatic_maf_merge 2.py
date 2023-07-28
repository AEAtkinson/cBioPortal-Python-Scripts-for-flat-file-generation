#!/usr/bin/env python3

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan  3 08:29:37 2019 through 2/24/2019

@author: aaronatkinson
"""
#import packages
import pandas as pd, numpy as np
import glob

df=[]
filenames = glob.glob('*Final3.maf')
for filename in filenames:
    columns = ['Hugo_Symbol', 'Entrez_Gene_Id', 'Center', 'NCBI_Build', 'Chromosome', 'Start_Position', 'End_Position', 'Strand', 'Variant_Classification', 'Variant_Type', 'Reference_Allele', 'Tumor_Seq_Allele1', 'Tumor_Seq_Allele2', 'dbSNP_RS', 'dbSNP_Val_Status', 'Tumor_Sample_Barcode', 'Matched_Norm_Sample_Barcode', 'Match_Norm_Seq_Allele1', 'Match_Norm_Seq_Allele2', 'Tumor_Validation_Allele1', 'Tumor_Validation_Allele2', 'Match_Norm_Validation_Allele1', 'Match_Norm_Validation_Allele2', 'Verification_Status', 'Validation_Status', 'Mutation_Status', 'Sequencing_Phase', 'Sequence_Source', 'Validation_Method', 'Score', 'BAM_File', 'Sequencer', 'Tumor_Sample_UUID', 'Matched_Norm_Sample_UUID', 'HGVSc', 'HGVSp', 'HGVSp_Short', 'Transcript_ID', 'Exon_Number', 't_depth', 't_ref_count', 't_alt_count', 'n_depth', 'n_ref_count', 'n_alt_count', 'all_effects', 'Allele', 'Gene', 'Feature', 'Feature_type', 'Consequence', 'cDNA_position', 'CDS_position', 'Protein_position', 'Amino_acids', 'Codons', 'Existing_variation', 'ALLELE_NUM', 'DISTANCE', 'STRAND_VEP', 'SYMBOL', 'SYMBOL_SOURCE', 'HGNC_ID', 'BIOTYPE', 'CANONICAL', 'CCDS', 'ENSP', 'SWISSPROT', 'TREMBL', 'UNIPARC', 'RefSeq', 'SIFT', 'PolyPhen', 'EXON', 'INTRON', 'DOMAINS', 'AF', 'AFR_AF', 'AMR_AF', 'ASN_AF', 'EAS_AF', 'EUR_AF', 'SAS_AF', 'AA_AF', 'EA_AF', 'CLIN_SIG', 'SOMATIC', 'PUBMED', 'MOTIF_NAME', 'MOTIF_POS', 'HIGH_INF_POS', 'MOTIF_SCORE_CHANGE', 'IMPACT', 'PICK', 'VARIANT_CLASS', 'TSL', 'HGVS_OFFSET', 'PHENO', 'MINIMISED', 'GENE_PHENO', 'FILTER', 'flanking_bps', 'vcf_id', 'vcf_qual', 'gnomAD_AF', 'gnomAD_AFR_AF', 'gnomAD_AMR_AF', 'gnomAD_ASJ_AF', 'gnomAD_EAS_AF', 'gnomAD_FIN_AF', 'gnomAD_NFE_AF', 'gnomAD_OTH_AF', 'gnomAD_SAS_AF', 'vcf_pos']
    df1 = pd.read_csv(filename, sep="\t", usecols=columns)
    print(*df1)
    df.append(df1)
df2 = pd.concat(df)
#Concat_table = pd.concat(df1)
cols = ['t_depth', 't_ref_count', 't_alt_count', 'n_depth', 'n_ref_count', 'n_alt_count']

df2[cols] = df2[cols].fillna(0)
df2[cols] = np.ceil(df2[cols]).astype(int)
#df2[cols] = df2[cols].astype(str)
#df2[cols] = df2[cols].replace(regex={r'.0': ''})
df2 = df2.sort_values(by=['Chromosome', 'Start_Position', 'Tumor_Sample_Barcode', 'Center'], ascending=False)
df2.to_csv('test.maf', sep='\t', index=False)
df2.drop_duplicates(subset=['Tumor_Sample_Barcode', 'Hugo_Symbol', 'Start_Position', 'End_Position'], inplace=True)
df2.to_csv('test2.maf', sep='\t', index=False)
df2 = df2.sort_values(by=['Chromosome', 'Start_Position'])
df2['NCBI_Build'] = 'hg38'
#df2['Matched_Norm_Sample_Barcode'] = df2['Tumor_Sample_Barcode'].str.split('_').str[0]
#df2['Matched_Norm_Sample_Barcode'] = df2['Matched_Norm_Sample_Barcode'].str.split('.').str[0]
#add Mutation_Status as either Germline or Somatic

df2.to_csv('CombinedACMG_MasterFilteredSortedFinal.maf', sep='\t', index=False)
