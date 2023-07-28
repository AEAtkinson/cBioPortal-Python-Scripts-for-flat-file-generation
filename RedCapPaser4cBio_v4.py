#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep  1 10:43:22 2021

@author: aaronatkinson
"""
#import packages
import pandas as pd 
import numpy as np
import glob, os
import re
from pathlib import Path
#import re

##Assign directory to change or use within directory
directory = "TempDataClinicalExtraction"
  
# Parent Directory path
path = '$pwd'
# Path
TempDir = os.path.join(path, directory)
  

##BUILD SAMPLE FILE
print('\n\n|||||||||||||||||||||||||||||||||||||||||||||| BUILDING SAMPLE FILE FOR cBIO ||||||||||||||||||||||||||||||||||||||||||||||\n\n')
t1 = []
FileSample = glob.glob('*_ClinicalMolLinkage_V4.csv')
for f in FileSample:
    dfstart = pd.read_csv(f, skiprows=[1], header=None)
    #dfstart.rename(columns={'ORIENAvatarKey': 'AvatarKey'}, inplace=True)
    t1.append(dfstart)
    dfc = pd.concat(t1, axis=0, ignore_index=True)
    dfc.columns = ['Patient_ID', 'Tumor/Germline', 'SAMPLE_ID', 'WES Batch', 'RNASeq', 'RNA Batch', 'Disease Type', 'SpecimenSiteOfOrigin', 'SpecimenSiteOfCollection', 'Primary/Met', 'SpecimenType', 'SpecimenDerivativeSource', 'Histology/Behavior', 'PreservationMethod', 'Percent Tumor Content', 'Age At Specimen Collection', 'SpecimenSiteOfOriginRollUp', 'SampleAgeInDays']
    dfc.columns = [c.replace(" ", "_") for c in list(dfc.columns)]
    dfc.rename(columns={'Tumor/Germline': 'Tumor_Germline'}, inplace=True)
    dfc = dfc.tail(-1)
    #cheat to remove unstripped header
    dfc = dfc.loc[~dfc['Patient_ID'].isin(['ORIENAvatarKey'])]
    print('\n\n|||||||||||||||||||||||||||||||||||||||||||||| CLINICAL FILE HEADERS dfc |||||||||||||||||||||||||||||||||||||||||||\n\n')
    #Drop germline data, but correct germline only samples first with germline SL number from sample sheet
    dfc.loc[dfc.Patient_ID == 'BCZK5PTSG9', ['SAMPLE_ID']] = 'SL400287'
    dfc.loc[dfc.Patient_ID == 'BCZK5PTSG9', ['Tumor_Germline']] = 'Tumor'
    #dfc.loc[dfc.Patient_ID == '731EJZUG4W', ['SAMPLE_ID']] = 'SL413287'
    dfc.loc[dfc.Patient_ID == 'JLGA08C3K5', ['SAMPLE_ID']] = 'SL365551'
    dfc.loc[dfc.Patient_ID == 'JLGA08C3K5', ['Tumor_Germline']] = 'Tumor'
    dfc.to_csv('FindME1.txt', sep='\t', index=False)
    dfc.loc[dfc.Patient_ID == 'RL11WPPBNV', ['SAMPLE_ID']] = 'SL359266'
    dfc.loc[dfc.Patient_ID == 'X0W1ONC8B2', ['SAMPLE_ID']] = 'SL365542'
    dfc.loc[dfc.Patient_ID == 'XASH2D2L7T', ['SAMPLE_ID']] = 'SL390930'
    dfc.loc[dfc.Patient_ID == 'YF6JV0C9MC', ['SAMPLE_ID']] = 'SL400151'
    dfc.loc[dfc.Patient_ID == 'YSVK3UWBEX', ['SAMPLE_ID']] = 'SL429082'
    dfc = dfc.loc[dfc['Tumor_Germline'].isin(['Tumor'])]
    dfc = dfc.loc[~dfc['SpecimenType'].isin(['RNA'])]
    #drop samples without Sample_ID
    dfc['SAMPLE_ID'].replace('', np.nan, inplace=True)
    dfc.dropna(subset=['SAMPLE_ID'], inplace=True)
    dfc.rename(columns={'Age At Specimen Collection': 'Age_At_Specimen_Collection', 'WES Batch': 'wes_batch'}, inplace=True)
    dfc.drop(['Tumor_Germline', 'RNASeq', 'RNA_Batch'], axis=1, inplace = True)
    dfc.replace({'Age 90 or older': 91, 'Unknown/Not Applicable': ''}, inplace=True)
    #print(dfc.dtypes)
    #dfc['Age_At_Specimen_Collection'] = dfc['Age_At_Specimen_Collection'].astype(float)
    dfc['Age_At_Specimen_Collection'] = pd.to_numeric(dfc.Age_At_Specimen_Collection, errors='coerce')
    dfc['Age_At_Specimen_Collection'] = dfc['Age_At_Specimen_Collection'].fillna(0)
    #dfc.rename(columns={'Age_At_Specimen_Collection': 'AgeAtDiagnosis'}, inplace=True)
    #print(dfc.dtypes)
    dfc = dfc.sort_values(by=['Age_At_Specimen_Collection', 'Patient_ID'])
    #dfc.drop_duplicates(subset=['Patient_ID', 'SAMPLE_ID'], keep='last', inplace=True)
    
    print('\n\n|||||||||||||||||||||||||||||||||||||||||||||| CLINICAL FILE HEADERS dfc |||||||||||||||||||||||||||||||||||||||||||\n\n')
t2 = []
FileSample1 = glob.glob('*_Diagnosis_V4.csv')
for f in FileSample1:
    dfstart2 = pd.read_csv(f, index_col=None, header=0)
    dfstart2.rename(columns={'AvatarKey': 'PatientId'}, inplace=True)
    t2.append(dfstart2)
    dfc1 = pd.concat(t2, axis=0, ignore_index=True)
    dfc1.drop(['AgeAtFirstContactFlag', 'AgeAtDiagnosisFlag', 'AgeAtOtherStagingSystemFlag', 'YearOfOtherStagingSystem', 'AgeAtPerformanceStatusAtDiagnosisFlag'], axis=1, inplace = True)
    #dfc1.rename(columns={'AgeAtDiagnosis': 'Age_At_Specimen_Collection'}, inplace=True)
    print(*dfc1)
    #dfc1.replace({'Age 90 or older': 91, 'Unknown/Not Applicable': ''}, inplace=True)
    dfc1 = dfc1.loc[~dfc1['AgeAtFirstContact'].isin(['Unknown/Not Applicable'])]
    #dfc1['AgeAtDiagnosis'] = dfc1['AgeAtDiagnosis'].str.replace('Age 90 or older', '91')
    dfc1['AgeAtDiagnosis'] = dfc1['AgeAtDiagnosis'].astype(str).str.replace('Age 90 or older', '91')
    #dfc1['Age_At_Specimen_Collection'] = dfc1.Age_At_Specimen_Collection.add(dfc1.AgeAtFirstContact).str.replace('Unknown/Not Applicable', '')
    dfc1['AgeAtDiagnosis'] = dfc1['AgeAtDiagnosis'].str.replace('Unknown/Not Applicable', '')

    #dfc1['Age_At_Specimen_Collection'] = dfc1['Age_At_Specimen_Collection'].astype(float)
    dfc1['AgeAtDiagnosis'] = pd.to_numeric(dfc1.AgeAtDiagnosis, errors='coerce')
    dfc1['AgeAtDiagnosis'] = dfc1['AgeAtDiagnosis'].fillna(0)
    #print(dfc1.dtypes)
    dfc1 = dfc1.sort_values(by=['AgeAtDiagnosis', 'PatientId'])
    def camel_to_snake(str):
       return re.sub(r'(?<!^)(?=[A-Z])', '_', str)
    dfc1.columns = [camel_to_snake(col) for col in dfc1.columns]
    dfc1.rename(columns={'Patient_Id': 'Patient_ID'}, inplace=True)
    #dfc1 = dfc1.loc[~dfc1['Age_At_Specimen_Collection'].isin([''])]
#merge Diagnosis and mol linkage by nearest date
SampleCDE = pd.merge_asof(dfc, dfc1, left_on='Age_At_Specimen_Collection', right_on='Age_At_Diagnosis', by='Patient_ID', direction='nearest')

#SampleCDE.to_csv('Temp1.txt', sep='\t', index=False)

#SampleCDE.to_csv('Temp2.txt', sep='\t', index=False)
SampleCDE.rename(columns={'primary/_met': 'primary_or_met', 't_n_m_classification_yc_t_n_m': 'tnm_classification_yctnm', 't_n_m_classification_a_t_n_m': 'tnm_classification_atnm', 't_n_m_classification_r_t_n_m': 'tnm_classification_rtnm', 't_n_m_classification_yp_t_n_m': 'tnm_classification_yptnm', 't_n_m_edition_number': 'tnm_edition_number'}, inplace=True)
print('\n\n************************************ SL_ID CHECK ***********************\n')

#Cleanup strings and numbers
SampleCDE['Percent_Tumor_Content'] = SampleCDE['Percent_Tumor_Content'].str.replace('> 30', '31')
SampleCDE.replace({'Age 90 or older': 91}, inplace=True)
#Drop duplicates - need to reduce duplicates by match on sample ?Age at diagnosis or AgeAtFirstContact with pd.merge_asof Age At Specimen Collection
#SampleCDE.drop_duplicates(subset=['PATIENT_ID', 'SAMPLE_ID'], keep='first', inplace=True)
SampleCDE.rename(columns = str.upper , inplace = True)
SampleCDE['AGE_AT_PERFORMANCE_STATUS_AT_DIAGNOSIS'] = SampleCDE['AGE_AT_PERFORMANCE_STATUS_AT_DIAGNOSIS'].str.replace('Unknown/Not Applicable', '')
SampleCDE['YEAR_OF_PERFORMANCE_STATUS_AT_DIAGNOSIS'] = SampleCDE['YEAR_OF_PERFORMANCE_STATUS_AT_DIAGNOSIS'].str.replace('Unknown/Not Applicable', '')
SampleCDE['AGE_AT_OTHER_STAGING_SYSTEM'] = SampleCDE['AGE_AT_OTHER_STAGING_SYSTEM'].str.replace('Unknown/Not Applicable', '')
SampleCDE.rename(columns={'SAMPLE_ID': 'SL_ID'}, inplace=True)
#SampleCDE.drop(['AGE_AT_DIAGNOSIS'], axis=1, inplace = True)
SL_SAMPLE = pd.read_csv('SAMPLE_ID_TMB_MSI.txt', sep='\t')
SampleCDE = SampleCDE.merge(SL_SAMPLE, on='SL_ID', how='right')
SampleCDE.rename(columns={'SPECIMENSITEOFORIGIN': 'SPECIMEN_SITE_OF_ORIGIN', 'SPECIMENSITEOFCOLLECTION': 'SPECIMEN_SITE_OF_COLLECTION', 'PRIMARY/MET': 'PRIMARY_OR_MET', 'SPECIMENTYPE': 'SPECIMEN_TYPE', 'SPECIMENDERIVATIVESOURCE': 'SPECIMEN_DERIVATIVE_SOURCE', 'HISTOLOGY/BEHAVIOR': 'HISTOLOGY/_BEHAVIOR', 'PRESERVATIONMETHOD': 'PRESERVATION_METHOD', 'SPECIMENSITEOFORIGINROLLUP': 'SPECIMEN_SITE_OF_ORIGIN_ROLL_UP', 'SAMPLEAGEINDAYS': 'SAMPLE_AGE_IN_DAYS', 'T_N_M_CLASSIFICATION_A_T_N_M': 'TNM_CLASSIFICATION_ATNM', 'T_N_M_CLASSIFICATION_R_T_N_M': 'TNM_CLASSIFICATION_RTNM', 'T_N_M_CLASSIFICATION_YC_T_N_M': 'TNM_CLASSIFICATION_YCTNM', 'T_N_M_CLASSIFICATION_YP_T_N_M': 'TNM_CLASSIFICATION_YPTNM', 'T_N_M_EDITION_NUMBER': 'TNM_EDITION_NUMBER'}, inplace=True)
SampleCDE = SampleCDE[['PATIENT_ID', 'SAMPLE_ID', 'SL_ID', 'WES_BATCH', 'DISEASE_TYPE', 'SPECIMEN_SITE_OF_ORIGIN', 'SPECIMEN_SITE_OF_COLLECTION', 'PRIMARY_OR_MET', 'SPECIMEN_TYPE', 'SPECIMEN_DERIVATIVE_SOURCE', 'HISTOLOGY/_BEHAVIOR', 'PRESERVATION_METHOD', 'PERCENT_TUMOR_CONTENT', 'AGE_AT_SPECIMEN_COLLECTION', 'SPECIMEN_SITE_OF_ORIGIN_ROLL_UP', 'SAMPLE_AGE_IN_DAYS', 'AGE_AT_FIRST_CONTACT', 'YEAR_OF_FIRST_CONTACT', 'AGE_AT_DIAGNOSIS', 'YEAR_OF_DIAGNOSIS', 'PRIMARY_DIAGNOSIS_SITE_CODE', 'PRIMARY_DIAGNOSIS_SITE', 'LATERALITY', 'HISTOLOGY_CODE', 'HISTOLOGY', 'GRADE_CLINICAL', 'GRADE_PATHOLOGICAL', 'GRADE_POST_THERAPY', 'HEM_MALIG_IMMUNOPHENOTYPE', 'HEM_MALIG_PHASE', 'CLIN_T_STAGE', 'CLIN_N_STAGE', 'CLIN_M_STAGE', 'CLIN_GROUP_STAGE', 'PATH_T_STAGE', 'PATH_N_STAGE', 'PATH_M_STAGE', 'PATH_GROUP_STAGE', 'TNM_CLASSIFICATION_ATNM', 'TNM_CLASSIFICATION_RTNM', 'TNM_CLASSIFICATION_YCTNM', 'TNM_CLASSIFICATION_YPTNM', 'TNM_EDITION_NUMBER', 'AGE_AT_OTHER_STAGING_SYSTEM', 'OTHER_STAGING_SYSTEM', 'OTHER_STAGING_VALUE', 'CURRENTLY_SEEN_FOR_PRIMARY_OR_RECURR', 'AGE_AT_PERFORMANCE_STATUS_AT_DIAGNOSIS', 'YEAR_OF_PERFORMANCE_STATUS_AT_DIAGNOSIS', 'PERFORM_STATUS_AT_DIAGNOSIS', 'PERFORM_STATUS_AT_DIAGNOSIS_SCALE', 'MSI_STATUS', 'TMB', 'SEQUENCING_PANEL']]
SampleCDE = SampleCDE.sort_values(by=['SAMPLE_ID', 'SAMPLE_AGE_IN_DAYS'])
SampleCDE.drop_duplicates(subset=['SAMPLE_ID'], keep='first', inplace=True)
#MSI = pd.read_csv('CompiledMsiResults.tsv', sep='\t', header=0)
#MSI.rename(columns={'Status': 'MSI_STATUS'}, inplace=True)
#SampleCDE = SampleCDE.merge(MSI, left_on='SAMPLE_ID', right_on='HCI', how='outer')
#TMBtxt = pd.read_csv('Concat_table.txt', sep='\t', header=0)
#TMBtxt['TMB'] = pd.to_numeric(TMBtxt.TMB, errors='coerce')
#SampleCDE = SampleCDE.merge(TMBtxt, left_on='SAMPLE_ID', right_on='SAMPLE_ID', how='left')
#SampleCDE['TMB'] = SampleCDE['TMB'].fillna(0)
#SampleCDE.drop(['HCI'], axis=1, inplace = True)
print('\n\n************************************ SAMPLE_ID CHECK ***********************\n')
#print(SampleCDE['SAMPLE_ID'])
###Merge in LoH
#LoH = pd.read_csv('uroLoHResults.txt', sep='\t')
##LoH['SAMPLE_ID'] = LoH['SAMPLE_ID'].replace('-', '_')
#SampleCDE = SampleCDE.merge(LoH, on='SAMPLE_ID', how='left')
#SampleCDE.insert(1, column='SAMPLE_ID_2', value='')
#SampleCDE.rename(columns={1: 'SAMPLE_ID_2'}, inplace=True)
#SampleCDE['SAMPLE_ID_2'] = SampleCDE['PATIENT_ID']+"_"+SampleCDE['SAMPLE_ID']
#SampleCDE.rename(columns={'SAMPLE_ID': 'SL_ID', 'SAMPLE_ID_2': 'SAMPLE_ID'}, inplace=True)
#add cBio specific header meta data
SampleCDE3meta = SampleCDE.dtypes
SampleCDE = SampleCDE.append(SampleCDE3meta, ignore_index=True)
cols = SampleCDE.columns.tolist()
SampleCDE.loc[-1] = cols
SampleCDE = SampleCDE.reindex(np.roll(SampleCDE.index, shift=1))
SampleCDE.loc[-1] = ''  # adding a row
SampleCDE = SampleCDE.reindex(np.roll(SampleCDE.index, shift=1))
SampleCDE.to_csv('Temp.txt', sep='\t', index=False)
SampleCDE.to_csv('SampleCDE.txt', sep='\t', index=False)

print('\n\n|||||||||||||||||||||||||||||||||||||||||||||| COMPLETE ---- BUILT CLINICAL FILE FOR cBIO |||||||||||||||||||||||||||||||||||||||||||\n\n')


##BUILD PATIENT FILES
print('\n******************************************* BUILDING PATIENT FILE FOR cBIO *********************************************\n\n')
#FilePatient = glob.glob('*_PatientMaster_V4.csv')
#for f in FilePatient:
#    df = pd.read_csv(f)
#    print(*df)
#    print('\n\n')
t3 = []
FilePatient = glob.glob('*_PatientMaster_V4.csv')
for f in FilePatient:
    dfstart3 = pd.read_csv(f, index_col=None, header=0)
    #dfstart3.rename(columns={'ORIENAvatarKey': 'AvatarKey'}, inplace=True)
    t3.append(dfstart3)
    df = pd.concat(t3, axis=0, ignore_index=True)
    print('\n\n')
#FilePatient1 = glob.glob('*_PatientHistory_V4.csv')
#for f in FilePatient1:
#    df1 = pd.read_csv(f)
#    print(*df1)    
t4 = []
FilePatient1 = glob.glob('*_PatientHistory_V4.csv')
for f in FilePatient1:
    dfstart4 = pd.read_csv(f, index_col=None, header=0)
    #dfstart3.rename(columns={'ORIENAvatarKey': 'AvatarKey'}, inplace=True)
    t4.append(dfstart4)
    df1 = pd.concat(t4, axis=0, ignore_index=True)
    print('\n\n')  
#df1.to_csv('PatientTemp.txt', sep='\t', index=False)
t5 = []
FilePatient2 = glob.glob('*_FamilyHistory_V4.csv')
for f in FilePatient2:
    dfstart5 = pd.read_csv(f, index_col=None, header=0)
    dfstart3.rename(columns={'ORIENAvatarKey': 'AvatarKey'}, inplace=True)
    t5.append(dfstart5)
    df2 = pd.concat(t5, axis=0, ignore_index=True)
    df2b = df2
    df2 = df2[['AvatarKey', 'CancerInFamilyInd']]
    df2.drop_duplicates(subset=['AvatarKey'], keep='first', inplace=True)
    print(*df2)
    print('\n******************************************* Parse Family Hx *********************************************\n\n')    
    df2b_father = df2b.loc[df2b['FamilyRelation'].isin(['Father'])]
    df2b_father['FamilyRelation'] = df2b_father['FamilyRelation'] + '_' + df2b_father['FamilyRelationPrimarySite'] + '(' + df2b_father['FamilyRelationPrimarySiteCode'] +')' + '_Age:' + df2b_father['FamilyRelationAgeGroup']
    df2b_father.rename(columns={'FamilyRelation': 'paternal_hx'}, inplace=True)
    df2b_father = df2b_father[['AvatarKey', 'paternal_hx']]
    df2b_father.drop_duplicates(subset=['AvatarKey'], keep='first', inplace=True)
    
    df2b_mother = df2b.loc[df2b['FamilyRelation'].isin(['Mother'])]
    df2b_mother['FamilyRelation'] = df2b_mother['FamilyRelation'] + '_' + df2b_mother['FamilyRelationPrimarySite'] + '(' + df2b_mother['FamilyRelationPrimarySiteCode'] +')' + '_Age:' + df2b_mother['FamilyRelationAgeGroup']
    df2b_mother.rename(columns={'FamilyRelation': 'maternal_hx'}, inplace=True)
    df2b_mother = df2b_mother[['AvatarKey', 'maternal_hx']]
    df2b_mother.drop_duplicates(subset=['AvatarKey'], keep='first', inplace=True)
    
    df2b_brother = df2b.loc[df2b['FamilyRelation'].isin(['Brother'])]
    df2b_brother['FamilyRelation'] = df2b_brother['FamilyRelation'] + '_' + df2b_brother['FamilyRelationPrimarySite'] + '(' + df2b_brother['FamilyRelationPrimarySiteCode'] +')' + '_Age:' + df2b_brother['FamilyRelationAgeGroup']
    df2b_brother.rename(columns={'FamilyRelation': 'fraternal_hx'}, inplace=True)
    df2b_brother = df2b_brother[['AvatarKey', 'fraternal_hx']]
    df2b_brother.drop_duplicates(subset=['AvatarKey'], keep='first', inplace=True)
    
    df2b_sister = df2b.loc[df2b['FamilyRelation'].isin(['Sister'])]
    df2b_sister['FamilyRelation'] = df2b_sister['FamilyRelation'] + '_' + df2b_sister['FamilyRelationPrimarySite'] + '(' + df2b_sister['FamilyRelationPrimarySiteCode'] +')' + '_Age:' + df2b_sister['FamilyRelationAgeGroup']
    df2b_sister.rename(columns={'FamilyRelation': 'sororal_hx'}, inplace=True)
    df2b_sister = df2b_sister[['AvatarKey', 'sororal_hx']]
    df2b_sister.drop_duplicates(subset=['AvatarKey'], keep='first', inplace=True)
    
    #Combine Indications with Family Hx
    df2 = df2.merge(df2b_father, on='AvatarKey', how='left').merge(df2b_mother, on='AvatarKey', how='left').merge(df2b_brother, on='AvatarKey', how='left').merge(df2b_sister, on='AvatarKey', how='left')
    
print('\n******************************************* IMPORT VITAL STATUS FILE FOR OUTCOMES *********************************************\n\n')
t6 = []
FilePatient3 = glob.glob('*_VitalStatus_V4.csv')
for f in FilePatient3:
    dfstart6 = pd.read_csv(f, index_col=None, header=0)
    #dfstart3.rename(columns={'ORIENAvatarKey': 'AvatarKey'}, inplace=True)
    t6.append(dfstart6)
    df3 = pd.concat(t6, axis=0, ignore_index=True)
    print(*df3)
#FilePatient3 = glob.glob('*_VitalStatus_V4.csv')
#for f in FilePatient3:
#    df3 = pd.read_csv(f)
#    print(*df3)
    df3.drop(['AgeAtLastContactFlag', 'AgeAtDeathFlag', 'VitalStatusConfirmed'], axis=1, inplace = True)
    df3['VitalStatus'] = df3['VitalStatus'].str.replace('Alive', 'LIVING')
    df3['VitalStatus'] = df3['VitalStatus'].str.replace('Dead', 'DECEASED')
    #Should I do this confirm with others:
    df3['VitalStatus'] = df3['VitalStatus'].str.replace('Lost to follow-up', 'LIVING')
    df3.replace({'Age 90 or older': 91}, inplace=True)
    df3['AgeAtDeath'] = df3['AgeAtDeath'].str.replace('Unknown/Not Applicable', '')
    df3['AgeAtDeath'] = df3['AgeAtDeath'].replace(np.nan, '', regex=True)
    df3[['AgeAtDeath']] = df3[['AgeAtDeath']].astype(float, errors= 'ignore')
    df3.rename(columns={'VitalStatus': 'OsStatus'}, inplace=True)
    #Note keep columns as SnakeCase so later universal reformat doesn't add extra characters. 
    df3['OsMonths'] = ''
    df3['MonthsLastFollowup'] = ''
    print(*df3)
    #print(df3.dtypes)
    print('\n\n********************************************* VITAL STATUS **********************************************************************\n\n')
#merge and clean patient files 
PatientCDE = df.merge(df1, on='AvatarKey', how='left').merge(df2, on='AvatarKey', how='left').merge(df3, on='AvatarKey', how='left')
PatientCDE = PatientCDE.sort_values(by=['AvatarKey'])
PatientCDE.drop_duplicates(subset=['AvatarKey'], keep='first', inplace=True)
#PatientCDE.to_csv('PatientTemp.txt', sep='\t', index=False)
#Change CamelCase to SNAKE_CASE for cbio requirement of all uppercase headers
def camel_to_snake(str):
   return re.sub(r'(?<!^)(?=[A-Z])', '_', str).lower()

PatientCDE.columns = [camel_to_snake(col) for col in PatientCDE.columns]
PatientCDE.rename(columns = str.upper , inplace = True)
#fix other columns
PatientCDE.replace({'_N_O_S': '_NOS', 'Age 90 or older': 91}, inplace=True)
PatientCDE.rename(columns={'AVATAR_KEY': 'PATIENT_ID', 'PRE_EXIST_H_G_S_I_L': 'PRE_EXIST_HGSIL', 'PRE_EXIST_H_P_V': 'PRE_EXIST_HPV', 'PRE_EXIST_A_D_H': 'PRE_EXIST_ADH', 'PRE_EXIST_A_L_H': 'PRE_EXIST_ALH', 'PRE_EXIST_A_I_N': 'PRE_EXIST_AIN', 'PRE_EXIST_C_I_N': 'PRE_EXIST_CIN', 'PRE_EXIST_L_I_N': 'PRE_EXIST_LIN', 'PRE_EXIST_P_A_I_N': 'PRE_EXIST_PAIN', 'PRE_EXIST_P_I_N': 'PRE_EXIST_PIN', 'PRE_EXIST_V_A_I_N': 'PRE_EXIST_VAIN', 'PRE_EXIST_V_I_N': 'PRE_EXIST_VIN', 'PRE_EXIST_I_E_N': 'PRE_EXIST_IEN', 'PRE_EXIST_M_G_U_S': 'PRE_EXIST_MGUS', 'PRE_EXIST_M_E_N': 'PRE_EXIST_MEN', 'PRE_EXIST_P_J_S': 'PRE_EXIST_PJS', 'PRE_EXIST_P_A_M': 'PRE_EXIST_P_A_M', 'C_O_P_D_DIAGNOSIS_IND': 'COPD_DIAGNOSIS_IND', 'C_V_A_STROKE_DIAGNOSIS_IND': 'CVA_STROKE_DIAGNOSIS_IND', 'D_V_T_DIAGNOSIS_IND': 'DVT_DIAGNOSIS_IND', 'H_I_V_A_I_D_S_DIAGNOSIS_IND': 'HIV_AIDS_DIAGNOSIS_IND', 'M_I_HEART_FAILURE_DIAGNOSIS_IND': 'MI_HEART_FAILURE_DIAGNOSIS_IND'}, inplace=True)
print('\n\n\n############################################################ PATIENT HEADERS ###############################################\n\n\n')
print(*PatientCDE)
print ('\n\n\n\n\n')
#Grab age at diagnosis and merge by patient on closest "AGE_AT_CLINICAL_RECORD_CREATION"
df4 = pd.read_csv('SampleCDE.txt', sep='\t')
df4 = df4[['PATIENT_ID', 'AGE_AT_DIAGNOSIS']]
df4.to_csv('df4.txt', sep='\t', index=False)
df4 = df4.sort_values(by=['AGE_AT_DIAGNOSIS'])
df4.drop_duplicates(subset=['PATIENT_ID'], keep='first', inplace=True)
df4['AGE_AT_DIAGNOSIS'] = pd.to_numeric(df4.AGE_AT_DIAGNOSIS, errors='coerce')

df4.to_csv('df4.txt', sep='\t', index=False)

PatientCDE['AGE_AT_LAST_UPDATE'] = pd.to_numeric(PatientCDE.AGE_AT_LAST_UPDATE, errors='coerce')


#Grab dates for outcomes, meds other timeline data
PatientCDE['AGE_AT_CLINICAL_RECORD_CREATION'] = pd.to_numeric(PatientCDE.AGE_AT_CLINICAL_RECORD_CREATION, errors='coerce')
PatientCDE = PatientCDE.sort_values(by=['AGE_AT_CLINICAL_RECORD_CREATION'])
PatientCDE = pd.merge(PatientCDE, df4, left_on = 'PATIENT_ID', right_on='PATIENT_ID', how='left')
#PatientCDE = pd.merge_asof(PatientCDE, df4, left_on = 'AGE_AT_CLINICAL_RECORD_CREATION', right_on='AGE_AT_DIAGNOSIS', by='PATIENT_ID', direction='backward')
#PatientCDE.to_csv('PatientCDEmerge.txt', sep='\t', index=False)
PatientCDE['MONTHS_LAST_FOLLOWUP'] = (PatientCDE['AGE_AT_LAST_UPDATE'] - PatientCDE['AGE_AT_DIAGNOSIS'])*12
PatientCDE['AGE_AT_DEATH'] = pd.to_numeric(PatientCDE.AGE_AT_DEATH, errors='coerce')
PatientCDE['OS_MONTHS'] = (PatientCDE['AGE_AT_DEATH'] - PatientCDE['AGE_AT_DIAGNOSIS'])*12
sub_list = ['LIVING']
PatientCDE.loc[PatientCDE['OS_STATUS'].isin(sub_list), 'OS_MONTHS'] = PatientCDE['MONTHS_LAST_FOLLOWUP']

#Patient Filter
#PatientCDE = PatientCDE[PatientCDE['PATIENT_ID'].isin(['03GCDX0K5N', '06B7Y294XL', '07EF826VTZ', '0C6NW2TC07', '0G5X064EB6', '0GFLCECW75', '0MJYHKZ6AK', '0RA5URUMVK', '0RFVUNS0M6', '0XBPD7NGNT', '15B16LU7OM', '1H50O85Z25', '1HJFO6EFWS', '1IQTFPMLQR', '1K0JG7QIX8', '1K0JG7QIX8', '1Q58ISCSYD', '1RQGB4W7T4', '1YUCT7U163', '274XK9OLE5', '282INULRAK', '28MTFJ7TIV', '2C9ADE4KW6', '2HNKJUZDQR', '2JRJWC9O7M', '2L2ZRV1LZY', '2N04GOV5XK', '2R23NCXQE0', '2UN4T04H08', '37SHF80OE0', '3C2HGYMBKU', '3E15DF0JVG', '3TLQLPNZP8', '41PP6CPRV7', '45PTETBCPX', '45VCVC3MKN', '4AFTCY1HP0', '4B7HJMVOOV', '4FLIDTYACP', '4FLIDTYACP', '4FTJM67GLA', '4OGM5KMQX2', '4S1A75B0CZ', '4UOOLFGTJC', '4UY5RRTZDA', '506HNR0L9Y', '56QVAODIUI', '5CE11AIS9P', '5JWOSJEGEQ', '5UUQPS8Y5W', '5VD97K8TUD', '63JFZJ9HYR', '64F1NAL85E', '66LL4MEA0D', '69NEPW50R9', '6A9ORBLX71', '6DRRKFWIL1', '6JD9VGBWZB', '6PX4GRURXV', '6T5GP2FH7P', '6U163EB0I2', '70HLF9IQQQ', '70P236XSX2', '731EJZUG4W', '73BBJCCFCG', '75HWMQ9THU', '75QKEG7QT5', '761ATF11EM', '78SU3ENGQ3', '78SU3ENGQ3', '797KYKKNTB', '7DQNA9XN5I', '7JT8CIRL35', '7O8FXBKCMP', '7P2075M6NA', '7XZV89ZKWL', '80II7FE25J', '80TN2OH9U2', '8A9JJV56BE', '8C9KM403EG', '8I52ZJY44W', '8IPGJJGS93', '8K3F2U4AEG', '8SSG3I43TG', '8SWRWW6ZR8', '8U1UG57FYC', '8XOABM89XD', '8YEVOS4E0I', '8ZVLI42YTE', '913EQAJLBI', '923EC0J63U', '928IO9LH0Z', '960N6W495Z', '988C2SHC0G', '995JAM1PGX', '9BZRRRQBLI', '9DFUAVXQTU', '9FUJS5X3V0', '9H5WMVOS0K', '9K443LZVX9', '9RQJ8JFILJ', '9SKYB6E4B4', '9ZXCC6ULKK', '9ZXCC6ULKK', 'A049JYM87K', 'A3AM9T7VPE', 'ACLOR9FFS5', 'ALL5VBRLOM', 'AP93ZGZX2W', 'AVMUPL2GS0', 'B3D07F2E2S', 'BCZK5PTSG9', 'BK6DKQTSI5', 'BTHE1B9DZ2', 'BU8GCS3MYT', 'BV2N69FYI6', 'BV2N69FYI6', 'BXXKG0HZIK', 'C5A7KUNVCS', 'CADAN6NOJ4', 'CGDPFH8IBW', 'CK9QL0A5ET', 'CM03MT997O', 'CM03MT997O', 'CM03MT997O', 'CTGHR4UZ18', 'CW0ZVIVU0P', 'CZE113UZIR', 'D0WP0B97Q7', 'D2C5F4UMXY', 'D7KAPTXFOG', 'DCPKIAYXDB', 'DQF8CHSKPG', 'DYKRF4KH3H', 'E2H89HF5KZ', 'E4H964SOU9', 'EF7AOZI79O', 'EHOVA8D07G', 'EKSOBYXBVG', 'EKXLPR0080', 'EMVBZPMWST', 'EOZ9LELKBM', 'EQNGL8N9W7', 'EU4WQTA1KL', 'F1XUR9R1RK', 'F2TQDDM3HH', 'F905RYZFIO', 'FNVA6L545Q', 'FOTCBZ3OMJ', 'FPLKDOPOH3', 'FVZV8LMNI5', 'FXGR65TB4V', 'FXSR99KNNT', 'GA7GNCV7KM', 'GF6S6LZVF5', 'GH6KKNEPWX', 'GIJZDVWC2Q', 'GNJVFE3HRM', 'GNXLDBYWN1', 'GOYPZUZNAW', 'GPLOBS1HRX', 'GU9WE974XH', 'GWNOJYMTTB', 'HG7EO5096W', 'HUCN0XE91R', 'HWECMUCW0Z', 'HWIG0FSC6Z', 'HXKS3QWP6F', 'HY2ZUQI0SD', 'I2EM49R0AO', 'I65GB8ZEJ2', 'IMH5IF30JJ', 'IW3OCCZBPJ', 'IW7K4RUJG9', 'IXEPJTR5B9', 'J0D6I815XU', 'J15Q9C8MIE', 'J3C9QS0P78', 'J60EKEW4OF', 'J81JUVHO9L', 'J9NHJFS8IN', 'JBLZPPZSFH', 'JG2Z3K6UC5', 'JHBPEVZ48F', 'JHZSJNIKIP', 'JJGV4LPENY', 'JLGA08C3K5', 'JM3MBDWG5H', 'JNXYDUFEED', 'JPSK03CCDH', 'JPSK03CCDH', 'JQHVTFHS1X', 'JRC4RTZNTK', 'JSWZX87JYV', 'JZKKZ1OTQJ', 'K1DVC6O9C7', 'K60LO6G1FG', 'K8DFO8FTT1', 'K8OFSLEHIX', 'KAN1H1VWQG', 'KC8HPZPX7T', 'KPCHILTX70', 'KSTUVZCJIF', 'KW650F0557', 'L410MYBLDN', 'L6RMSU0BXG', 'L8HCOD7F64', 'LATEJ509E1', 'LLXVLRWEAE', 'LM4ZZPQQFO', 'LS6UIWKFEJ', 'LX5K9JZKLH', 'LXELBL78SE', 'M1D0DA5RWT', 'M9G3Y6A4X7', 'MG1SUQR5K8', 'MIVJ67B9F5', 'MJYLGG6W0K', 'MRZW12GMST', 'MYG3XG6SAQ', 'N4K1N7UUMM', 'N675BI25Y6', 'N8USV2HKJ6', 'NAMDA2GZXA', 'NAR8MDP01R', 'NDV4VM27O6', 'NE17FOVUC5', 'NK65MD354P', 'NOD926GZ41', 'NPUNSPR11V', 'O8RFN0K0CC', 'OA1K7SKWYW', 'OFMS9RN5W8', 'OFOSHGCQE4', 'OJQ5IBJCDB', 'ONHYTARPID', 'OPI88VIT01', 'OQYBJZXP4D', 'ORL6LF2IMD', 'OWNR0K9J5O', 'OXKYPBVO0Q', 'OXKYPBVO0Q', 'P0YOHVOGGY', 'P6F7GJOYOS', 'PCM0C3L6G1', 'PGVXTE18CS', 'PJ316U0AO6', 'PQDMEIJU8A', 'PTN0DBFVOH', 'Q1LFWVL6GK', 'Q1TM0Z8WWH', 'Q3XX8RO9CY', 'QB511DI6QC', 'QD6U65T2FH', 'QDTQI9DG2U', 'QJ80IGAAAA', 'QO27BH1WTV', 'QUEJJTOM01', 'QVEQP2FBPS', 'QVLKW6YJ9S', 'QZ68W571DZ', 'R0JJZPXHE6', 'R9EFJQNQD7', 'RC5DO9WG9A', 'RDNMW7PS5V', 'RL11WPPBNV', 'RN9YUHZDTE', 'RRWTZZSW3F', 'S414XWRR6V', 'S5AEWMYN2T', 'S5R5RW4L05', 'S7GI3EL7N1', 'S9CLR0UMQD', 'S9GIV41A0F', 'SF953NJFYV', 'SQV8WBN78P', 'STZJJQ7LKL', 'SUYGY0MFVB', 'SW0L8B0LS4', 'T3BB8UT9KP', 'T653MJNWUY', 'T8IWENLDFE', 'TEEDQQPZOI', 'TFJ3UKACLL', 'TG4OSS40W2', 'TGCXP2RE55', 'TTMBZYSHU8', 'TTN1Q2AA35', 'TTN1Q2AA35', 'TTUKMQJUXT', 'TX7XF0XMBL', 'TXJ9VXE4QC', 'TYUGO1UOO3', 'U0JP7MO0CU', 'U1KHQ8NI5I', 'U3YN74R7B9', 'U6SEJCBUKI', 'U86PUE8URG', 'UB8PDQ9GI1', 'UD8WP9Y59H', 'UFKCGRWOE6', 'ULH5Q1UQZ9', 'UMZ44YQOYC', 'UMZ44YQOYC', 'UOHLJNY3WH', 'UOWQ0WV182', 'UYBI9N0WBE', 'V3V21CTLZR', 'VE1W0KD4IA', 'VF0HIW4X63', 'VI0WMPQNAD', 'VIMLBEI2ZQ', 'VKOQ16NG0M', 'VKOQ16NG0M', 'VMKJ7V6SRD', 'VMKJ7V6SRD', 'VQS7PT887V', 'VV8VUVF9XA', 'W3W1ROPU95', 'W6XNV1LQB1', 'W7UAIA7T0S', 'W80YCK8SOZ', 'WCYEL3APQW', 'WIYGVFIKGG', 'WK78CI5LDO', 'WMJ41DPJ3Y', 'WPYHWCFES9', 'WR6V5879WH', 'WXJA24FUA7', 'WYMINM1LOF', 'X0W1ONC8B2', 'X1MS489ILB', 'X2AGMD6P59', 'X4HHLWE7SC', 'X4JDSBA0MD', 'X4JDSBA0MD', 'X5A3YUL1OI', 'X5WMJ9MUOV', 'X7566XEPLR', 'X8NG9V0HJ6', 'XASH2D2L7T', 'XG77M6C95Z', 'XLDFSI23XN', 'XTIR5SW5E8', 'Y0S85J0AZR', 'Y1DJA21YHZ', 'Y4HVM00UR9', 'YBIJBXII7J', 'YF6JV0C9MC', 'YP49T16HM1', 'YV85IDK7KQ', 'YY88X645NE', 'Z0075688AQ', 'Z0RE2PT1R9', 'Z1WO6PGN3E', 'Z9478QITNI', 'ZA6WQHDGBM', 'ZBG7MKFOZ1', 'ZC60EAOUMH', 'ZE7RW2IJMR', 'ZHRR08KKVM', 'ZIN73MJOAB', 'ZINL8265JQ', 'ZKG95SXJUU', 'ZKZ02PIDIU', 'ZM2Y4LL5QP', 'ZQRZQA3DZT', 'ZW84KJ7LJA'])]
#PatientCDE.drop_duplicates(subset=['PATIENT_ID', 'AGE_AT_CLINICAL_RECORD_CREATION'], keep='first', inplace=True)
Timeline = PatientCDE[['PATIENT_ID','AGE_AT_DIAGNOSIS', 'AGE_AT_LAST_UPDATE']]
PatientCDE.to_csv('PatientCDE.txt', sep='\t', index=False)
print('\n******************************************* BUILT PATIENT FILE FOR cBIO *********************************************\n\n')
PatientCDE = PatientCDE
print('\n\n|||||||||||||||||||||||||||||||||||||||||||||| BUILDING MEDICATION FILE FOR cBIO |||||||||||||||||||||||||||||||||||||||||||\n\n')
t7 = []
FilePatient4 = glob.glob('*_Medications_V4.csv')
for f in FilePatient4:
    dfstart7 = pd.read_csv(f, index_col=None, header=0)
    #dfstart3.rename(columns={'ORIENAvatarKey': 'AvatarKey'}, inplace=True)
    t7.append(dfstart7)
    dfm = pd.concat(t7, axis=0, ignore_index=True)
    #Patient Filter
    #dfm = dfm[dfm['AvatarKey'].isin(['03GCDX0K5N', '06B7Y294XL', '07EF826VTZ', '0C6NW2TC07', '0G5X064EB6', '0GFLCECW75', '0MJYHKZ6AK', '0RA5URUMVK', '0RFVUNS0M6', '0XBPD7NGNT', '15B16LU7OM', '1H50O85Z25', '1HJFO6EFWS', '1IQTFPMLQR', '1K0JG7QIX8', '1K0JG7QIX8', '1Q58ISCSYD', '1RQGB4W7T4', '1YUCT7U163', '274XK9OLE5', '282INULRAK', '28MTFJ7TIV', '2C9ADE4KW6', '2HNKJUZDQR', '2JRJWC9O7M', '2L2ZRV1LZY', '2N04GOV5XK', '2R23NCXQE0', '2UN4T04H08', '37SHF80OE0', '3C2HGYMBKU', '3TLQLPNZP8', '41PP6CPRV7', '45PTETBCPX', '45VCVC3MKN', '4AFTCY1HP0', '4B7HJMVOOV', '4FLIDTYACP', '4FLIDTYACP', '4FTJM67GLA', '4OGM5KMQX2', '4S1A75B0CZ', '4UOOLFGTJC', '4UY5RRTZDA', '506HNR0L9Y', '56QVAODIUI', '5CE11AIS9P', '5JWOSJEGEQ', '5UUQPS8Y5W', '5VD97K8TUD', '63JFZJ9HYR', '64F1NAL85E', '66LL4MEA0D', '69NEPW50R9', '6A9ORBLX71', '6DRRKFWIL1', '6JD9VGBWZB', '6PX4GRURXV', '6T5GP2FH7P', '6U163EB0I2', '70HLF9IQQQ', '70P236XSX2', '731EJZUG4W', '73BBJCCFCG', '75HWMQ9THU', '75QKEG7QT5', '761ATF11EM', '78SU3ENGQ3', '78SU3ENGQ3', '797KYKKNTB', '7DQNA9XN5I', '7JT8CIRL35', '7O8FXBKCMP', '7P2075M6NA', '7XZV89ZKWL', '80II7FE25J', '80TN2OH9U2', '8A9JJV56BE', '8C9KM403EG', '8I52ZJY44W', '8IPGJJGS93', '8K3F2U4AEG', '8SSG3I43TG', '8SWRWW6ZR8', '8U1UG57FYC', '8XOABM89XD', '8YEVOS4E0I', '8ZVLI42YTE', '913EQAJLBI', '923EC0J63U', '928IO9LH0Z', '960N6W495Z', '988C2SHC0G', '995JAM1PGX', '9BZRRRQBLI', '9DFUAVXQTU', '9FUJS5X3V0', '9H5WMVOS0K', '9K443LZVX9', '9RQJ8JFILJ', '9SKYB6E4B4', '9ZXCC6ULKK', '9ZXCC6ULKK', 'A049JYM87K', 'A3AM9T7VPE', 'ACLOR9FFS5', 'ALL5VBRLOM', 'AP93ZGZX2W', 'AVMUPL2GS0', 'B3D07F2E2S', 'BCZK5PTSG9', 'BK6DKQTSI5', 'BTHE1B9DZ2', 'BU8GCS3MYT', 'BV2N69FYI6', 'BV2N69FYI6', 'BXXKG0HZIK', 'C5A7KUNVCS', 'CADAN6NOJ4', 'CGDPFH8IBW', 'CK9QL0A5ET', 'CM03MT997O', 'CM03MT997O', 'CM03MT997O', 'CTGHR4UZ18', 'CW0ZVIVU0P', 'CZE113UZIR', 'D0WP0B97Q7', 'D2C5F4UMXY', 'D7KAPTXFOG', 'DCPKIAYXDB', 'DQF8CHSKPG', 'DYKRF4KH3H', 'E2H89HF5KZ', 'E4H964SOU9', 'EF7AOZI79O', 'EHOVA8D07G', 'EKSOBYXBVG', 'EKXLPR0080', 'EMVBZPMWST', 'EOZ9LELKBM', 'EQNGL8N9W7', 'EU4WQTA1KL', 'F1XUR9R1RK', 'F2TQDDM3HH', 'F905RYZFIO', 'FNVA6L545Q', 'FOTCBZ3OMJ', 'FPLKDOPOH3', 'FVZV8LMNI5', 'FXGR65TB4V', 'FXSR99KNNT', 'GA7GNCV7KM', 'GF6S6LZVF5', 'GH6KKNEPWX', 'GIJZDVWC2Q', 'GNJVFE3HRM', 'GNXLDBYWN1', 'GOYPZUZNAW', 'GPLOBS1HRX', 'GU9WE974XH', 'GWNOJYMTTB', 'HG7EO5096W', 'HUCN0XE91R', 'HWECMUCW0Z', 'HWIG0FSC6Z', 'HXKS3QWP6F', 'HY2ZUQI0SD', 'I2EM49R0AO', 'I65GB8ZEJ2', 'IMH5IF30JJ', 'IW3OCCZBPJ', 'IW7K4RUJG9', 'IXEPJTR5B9', 'J0D6I815XU', 'J15Q9C8MIE', 'J3C9QS0P78', 'J60EKEW4OF', 'J81JUVHO9L', 'J9NHJFS8IN', 'JBLZPPZSFH', 'JG2Z3K6UC5', 'JHBPEVZ48F', 'JHZSJNIKIP', 'JJGV4LPENY', 'JLGA08C3K5', 'JM3MBDWG5H', 'JNXYDUFEED', 'JPSK03CCDH', 'JPSK03CCDH', 'JQHVTFHS1X', 'JRC4RTZNTK', 'JSWZX87JYV', 'JZKKZ1OTQJ', 'K1DVC6O9C7', 'K60LO6G1FG', 'K8DFO8FTT1', 'K8OFSLEHIX', 'KAN1H1VWQG', 'KC8HPZPX7T', 'KPCHILTX70', 'KSTUVZCJIF', 'KW650F0557', 'L410MYBLDN', 'L6RMSU0BXG', 'L8HCOD7F64', 'LATEJ509E1', 'LLXVLRWEAE', 'LM4ZZPQQFO', 'LS6UIWKFEJ', 'LX5K9JZKLH', 'LXELBL78SE', 'M1D0DA5RWT', 'M9G3Y6A4X7', 'MG1SUQR5K8', 'MIVJ67B9F5', 'MJYLGG6W0K', 'MRZW12GMST', 'MYG3XG6SAQ', 'N4K1N7UUMM', 'N675BI25Y6', 'N8USV2HKJ6', 'NAMDA2GZXA', 'NAR8MDP01R', 'NDV4VM27O6', 'NE17FOVUC5', 'NK65MD354P', 'NOD926GZ41', 'NPUNSPR11V', 'O8RFN0K0CC', 'OA1K7SKWYW', 'OFMS9RN5W8', 'OFOSHGCQE4', 'OJQ5IBJCDB', 'ONHYTARPID', 'OPI88VIT01', 'OQYBJZXP4D', 'ORL6LF2IMD', 'OWNR0K9J5O', 'OXKYPBVO0Q', 'OXKYPBVO0Q', 'P0YOHVOGGY', 'P6F7GJOYOS', 'PCM0C3L6G1', 'PGVXTE18CS', 'PJ316U0AO6', 'PQDMEIJU8A', 'PTN0DBFVOH', 'Q1LFWVL6GK', 'Q1TM0Z8WWH', 'Q3XX8RO9CY', 'QB511DI6QC', 'QD6U65T2FH', 'QDTQI9DG2U', 'QJ80IGAAAA', 'QO27BH1WTV', 'QUEJJTOM01', 'QVEQP2FBPS', 'QVLKW6YJ9S', 'QZ68W571DZ', 'R0JJZPXHE6', 'R9EFJQNQD7', 'RC5DO9WG9A', 'RDNMW7PS5V', 'RL11WPPBNV', 'RN9YUHZDTE', 'RRWTZZSW3F', 'S414XWRR6V', 'S5AEWMYN2T', 'S5R5RW4L05', 'S7GI3EL7N1', 'S9CLR0UMQD', 'S9GIV41A0F', 'SF953NJFYV', 'SQV8WBN78P', 'STZJJQ7LKL', 'SUYGY0MFVB', 'SW0L8B0LS4', 'T3BB8UT9KP', 'T653MJNWUY', 'T8IWENLDFE', 'TEEDQQPZOI', 'TFJ3UKACLL', 'TG4OSS40W2', 'TGCXP2RE55', 'TTMBZYSHU8', 'TTN1Q2AA35', 'TTN1Q2AA35', 'TTUKMQJUXT', 'TX7XF0XMBL', 'TXJ9VXE4QC', 'TYUGO1UOO3', 'U0JP7MO0CU', 'U1KHQ8NI5I', 'U3YN74R7B9', 'U6SEJCBUKI', 'U86PUE8URG', 'UB8PDQ9GI1', 'UD8WP9Y59H', 'UFKCGRWOE6', 'ULH5Q1UQZ9', 'UMZ44YQOYC', 'UMZ44YQOYC', 'UOHLJNY3WH', 'UOWQ0WV182', 'UYBI9N0WBE', 'V3V21CTLZR', 'VE1W0KD4IA', 'VF0HIW4X63', 'VI0WMPQNAD', 'VIMLBEI2ZQ', 'VKOQ16NG0M', 'VKOQ16NG0M', 'VMKJ7V6SRD', 'VMKJ7V6SRD', 'VQS7PT887V', 'VV8VUVF9XA', 'W3W1ROPU95', 'W6XNV1LQB1', 'W7UAIA7T0S', 'W80YCK8SOZ', 'WCYEL3APQW', 'WIYGVFIKGG', 'WK78CI5LDO', 'WMJ41DPJ3Y', 'WPYHWCFES9', 'WR6V5879WH', 'WXJA24FUA7', 'WYMINM1LOF', 'X0W1ONC8B2', 'X1MS489ILB', 'X2AGMD6P59', 'X4HHLWE7SC', 'X4JDSBA0MD', 'X4JDSBA0MD', 'X5A3YUL1OI', 'X5WMJ9MUOV', 'X7566XEPLR', 'X8NG9V0HJ6', 'XASH2D2L7T', 'XG77M6C95Z', 'XLDFSI23XN', 'XTIR5SW5E8', 'Y0S85J0AZR', 'Y1DJA21YHZ', 'Y4HVM00UR9', 'YBIJBXII7J', 'YF6JV0C9MC', 'YP49T16HM1', 'YV85IDK7KQ', 'YY88X645NE', 'Z0075688AQ', 'Z0RE2PT1R9', 'Z1WO6PGN3E', 'Z9478QITNI', 'ZA6WQHDGBM', 'ZBG7MKFOZ1', 'ZC60EAOUMH', 'ZE7RW2IJMR', 'ZHRR08KKVM', 'ZIN73MJOAB', 'ZINL8265JQ', 'ZKG95SXJUU', 'ZKZ02PIDIU', 'ZM2Y4LL5QP', 'ZQRZQA3DZT', 'ZW84KJ7LJA', '3E15DF0JVG'])]
    #dfm2.to_csv('DiscordantPatientListIO.txt', sep='\t', index=False)
#FileSample1 = glob.glob('*_Medications_V4.csv')
#for f in FileSample1:
#    dfm = pd.read_csv(f)
    #cleanup unused columns
    dfm.drop(['MedicationInd', 'MedReasonNoneGiven', 'AgeAtMedStopFlag', 'AgeAtMedStartFlag'], axis=1, inplace = True)
    #Convert to snake_case
    dfm.columns = [camel_to_snake(col) for col in dfm.columns]
    dfm.rename(columns = str.upper , inplace = True)
    dfm.rename(columns={'AVATAR_KEY': 'PATIENT_ID', 'MEDICATION': 'AGENT(S)', 'MED_LINE_REGIMEN': 'SUBTYPE'}, inplace=True)
    dfm = dfm.merge(Timeline, on='PATIENT_ID', how='left')
    dfm.drop(['AGE_AT_LAST_UPDATE'], axis=1, inplace = True)
    dfm.replace({'Unknown/Not Applicable': '', 'Age 90 or older': 91}, inplace=True)
    ##Note start date cannot be blank
    dfm.insert(1, column='START_DATE', value=0)
    dfm.insert(2, column='STOP_DATE', value='')
    dfm.insert(3, column='EVENT_TYPE', value='TREATMENT')
    dfm.insert(4, column='TREATMENT_TYPE', value='Medical Therapy')
    #reorder columns
    dfm = dfm[['PATIENT_ID', 'START_DATE', 'STOP_DATE', 'EVENT_TYPE', 'TREATMENT_TYPE', 'AGENT(S)', 'SUBTYPE', 'MED_PRIMARY_DIAGNOSIS_SITE_CODE', 'MED_PRIMARY_DIAGNOSIS_SITE', 'YEAR_OF_MED_START', 'YEAR_OF_MED_STOP', 'MED_CONTINUING', 'CHANGE_OF_TREATMENT', 'SYSTEMIC_SURGERY_SEQUENCE', 'AGE_AT_MED_START', 'AGE_AT_MED_STOP', 'AGE_AT_DIAGNOSIS']]
    print(*dfm)
    #convert strings to floats for start and stop date calculations
    dfm['AGE_AT_MED_START'] = pd.to_numeric(dfm.AGE_AT_MED_START, errors='coerce')
    dfm['AGE_AT_MED_STOP'] = pd.to_numeric(dfm.AGE_AT_MED_STOP, errors='coerce')
    dfm['AGE_AT_DIAGNOSIS'] = pd.to_numeric(dfm.AGE_AT_DIAGNOSIS, errors='coerce')
    dfm['START_DATE'] = (dfm['AGE_AT_DIAGNOSIS'] - dfm['AGE_AT_MED_START'])*365.25
    dfm['deltaT'] = dfm['AGE_AT_MED_STOP'] - dfm['AGE_AT_MED_START']
    dfm['STOP_DATE'] = (dfm['AGE_AT_DIAGNOSIS'] - dfm['AGE_AT_MED_START'] +dfm['deltaT'])*365.25
    print('\n\n')
    dfm.to_csv('MedicationsPFS.txt', sep='\t', index=False)
    #print(dfm.dtypes)
    dfm.drop(['AGE_AT_MED_START', 'AGE_AT_MED_STOP', 'AGE_AT_DIAGNOSIS', 'deltaT'], axis=1, inplace = True)
    print('\n\n')
    #print(dfm.dtypes)
    #cut all times prior to creation of clinical record. 
    columns = ['START_DATE', 'STOP_DATE']
    for col in columns:
        dfm[col][dfm[col] < 0] = 0
    dfm['START_DATE'] = dfm['START_DATE'].round(decimals = 0)
    dfm['STOP_DATE'] = dfm['STOP_DATE'].round(decimals = 0)
    dfm.replace(r'^\s+$', np.nan, regex=True)
    ##Note start date cannot be blank
    dfm[['START_DATE']] = dfm[['START_DATE']].fillna('0')
    dfm[['STOP_DATE']] = dfm[['STOP_DATE']].fillna('')
    dfm['START_DATE'] = dfm['START_DATE'].astype(str).apply(lambda x: x.replace('.0',''))
    dfm['STOP_DATE'] = dfm['STOP_DATE'].astype(str).apply(lambda x: x.replace('.0',''))
    print('\n\n')
    #print(dfm.dtypes)
    #Patient Filter
    #dfm = dfm[dfm['PATIENT_ID'].isin(['03GCDX0K5N', '06B7Y294XL', '07EF826VTZ', '0C6NW2TC07', '0G5X064EB6', '0GFLCECW75', '0MJYHKZ6AK', '0RA5URUMVK', '0RFVUNS0M6', '0XBPD7NGNT', '15B16LU7OM', '1H50O85Z25', '1HJFO6EFWS', '1IQTFPMLQR', '1K0JG7QIX8', '1K0JG7QIX8', '1Q58ISCSYD', '1RQGB4W7T4', '1YUCT7U163', '274XK9OLE5', '282INULRAK', '28MTFJ7TIV', '2C9ADE4KW6', '2HNKJUZDQR', '2JRJWC9O7M', '2L2ZRV1LZY', '2N04GOV5XK', '2R23NCXQE0', '2UN4T04H08', '37SHF80OE0', '3C2HGYMBKU', '3E15DF0JVG', '3TLQLPNZP8', '41PP6CPRV7', '45PTETBCPX', '45VCVC3MKN', '4AFTCY1HP0', '4B7HJMVOOV', '4FLIDTYACP', '4FLIDTYACP', '4FTJM67GLA', '4OGM5KMQX2', '4S1A75B0CZ', '4UOOLFGTJC', '4UY5RRTZDA', '506HNR0L9Y', '56QVAODIUI', '5CE11AIS9P', '5JWOSJEGEQ', '5UUQPS8Y5W', '5VD97K8TUD', '63JFZJ9HYR', '64F1NAL85E', '66LL4MEA0D', '69NEPW50R9', '6A9ORBLX71', '6DRRKFWIL1', '6JD9VGBWZB', '6PX4GRURXV', '6T5GP2FH7P', '6U163EB0I2', '70HLF9IQQQ', '70P236XSX2', '731EJZUG4W', '73BBJCCFCG', '75HWMQ9THU', '75QKEG7QT5', '761ATF11EM', '78SU3ENGQ3', '78SU3ENGQ3', '797KYKKNTB', '7DQNA9XN5I', '7JT8CIRL35', '7O8FXBKCMP', '7P2075M6NA', '7XZV89ZKWL', '80II7FE25J', '80TN2OH9U2', '8A9JJV56BE', '8C9KM403EG', '8I52ZJY44W', '8IPGJJGS93', '8K3F2U4AEG', '8SSG3I43TG', '8SWRWW6ZR8', '8U1UG57FYC', '8XOABM89XD', '8YEVOS4E0I', '8ZVLI42YTE', '913EQAJLBI', '923EC0J63U', '928IO9LH0Z', '960N6W495Z', '988C2SHC0G', '995JAM1PGX', '9BZRRRQBLI', '9DFUAVXQTU', '9FUJS5X3V0', '9H5WMVOS0K', '9K443LZVX9', '9RQJ8JFILJ', '9SKYB6E4B4', '9ZXCC6ULKK', '9ZXCC6ULKK', 'A049JYM87K', 'A3AM9T7VPE', 'ACLOR9FFS5', 'ALL5VBRLOM', 'AP93ZGZX2W', 'AVMUPL2GS0', 'B3D07F2E2S', 'BCZK5PTSG9', 'BK6DKQTSI5', 'BTHE1B9DZ2', 'BU8GCS3MYT', 'BV2N69FYI6', 'BV2N69FYI6', 'BXXKG0HZIK', 'C5A7KUNVCS', 'CADAN6NOJ4', 'CGDPFH8IBW', 'CK9QL0A5ET', 'CM03MT997O', 'CM03MT997O', 'CM03MT997O', 'CTGHR4UZ18', 'CW0ZVIVU0P', 'CZE113UZIR', 'D0WP0B97Q7', 'D2C5F4UMXY', 'D7KAPTXFOG', 'DCPKIAYXDB', 'DQF8CHSKPG', 'DYKRF4KH3H', 'E2H89HF5KZ', 'E4H964SOU9', 'EF7AOZI79O', 'EHOVA8D07G', 'EKSOBYXBVG', 'EKXLPR0080', 'EMVBZPMWST', 'EOZ9LELKBM', 'EQNGL8N9W7', 'EU4WQTA1KL', 'F1XUR9R1RK', 'F2TQDDM3HH', 'F905RYZFIO', 'FNVA6L545Q', 'FOTCBZ3OMJ', 'FPLKDOPOH3', 'FVZV8LMNI5', 'FXGR65TB4V', 'FXSR99KNNT', 'GA7GNCV7KM', 'GF6S6LZVF5', 'GH6KKNEPWX', 'GIJZDVWC2Q', 'GNJVFE3HRM', 'GNXLDBYWN1', 'GOYPZUZNAW', 'GPLOBS1HRX', 'GU9WE974XH', 'GWNOJYMTTB', 'HG7EO5096W', 'HUCN0XE91R', 'HWECMUCW0Z', 'HWIG0FSC6Z', 'HXKS3QWP6F', 'HY2ZUQI0SD', 'I2EM49R0AO', 'I65GB8ZEJ2', 'IMH5IF30JJ', 'IW3OCCZBPJ', 'IW7K4RUJG9', 'IXEPJTR5B9', 'J0D6I815XU', 'J15Q9C8MIE', 'J3C9QS0P78', 'J60EKEW4OF', 'J81JUVHO9L', 'J9NHJFS8IN', 'JBLZPPZSFH', 'JG2Z3K6UC5', 'JHBPEVZ48F', 'JHZSJNIKIP', 'JJGV4LPENY', 'JLGA08C3K5', 'JM3MBDWG5H', 'JNXYDUFEED', 'JPSK03CCDH', 'JPSK03CCDH', 'JQHVTFHS1X', 'JRC4RTZNTK', 'JSWZX87JYV', 'JZKKZ1OTQJ', 'K1DVC6O9C7', 'K60LO6G1FG', 'K8DFO8FTT1', 'K8OFSLEHIX', 'KAN1H1VWQG', 'KC8HPZPX7T', 'KPCHILTX70', 'KSTUVZCJIF', 'KW650F0557', 'L410MYBLDN', 'L6RMSU0BXG', 'L8HCOD7F64', 'LATEJ509E1', 'LLXVLRWEAE', 'LM4ZZPQQFO', 'LS6UIWKFEJ', 'LX5K9JZKLH', 'LXELBL78SE', 'M1D0DA5RWT', 'M9G3Y6A4X7', 'MG1SUQR5K8', 'MIVJ67B9F5', 'MJYLGG6W0K', 'MRZW12GMST', 'MYG3XG6SAQ', 'N4K1N7UUMM', 'N675BI25Y6', 'N8USV2HKJ6', 'NAMDA2GZXA', 'NAR8MDP01R', 'NDV4VM27O6', 'NE17FOVUC5', 'NK65MD354P', 'NOD926GZ41', 'NPUNSPR11V', 'O8RFN0K0CC', 'OA1K7SKWYW', 'OFMS9RN5W8', 'OFOSHGCQE4', 'OJQ5IBJCDB', 'ONHYTARPID', 'OPI88VIT01', 'OQYBJZXP4D', 'ORL6LF2IMD', 'OWNR0K9J5O', 'OXKYPBVO0Q', 'OXKYPBVO0Q', 'P0YOHVOGGY', 'P6F7GJOYOS', 'PCM0C3L6G1', 'PGVXTE18CS', 'PJ316U0AO6', 'PQDMEIJU8A', 'PTN0DBFVOH', 'Q1LFWVL6GK', 'Q1TM0Z8WWH', 'Q3XX8RO9CY', 'QB511DI6QC', 'QD6U65T2FH', 'QDTQI9DG2U', 'QJ80IGAAAA', 'QO27BH1WTV', 'QUEJJTOM01', 'QVEQP2FBPS', 'QVLKW6YJ9S', 'QZ68W571DZ', 'R0JJZPXHE6', 'R9EFJQNQD7', 'RC5DO9WG9A', 'RDNMW7PS5V', 'RL11WPPBNV', 'RN9YUHZDTE', 'RRWTZZSW3F', 'S414XWRR6V', 'S5AEWMYN2T', 'S5R5RW4L05', 'S7GI3EL7N1', 'S9CLR0UMQD', 'S9GIV41A0F', 'SF953NJFYV', 'SQV8WBN78P', 'STZJJQ7LKL', 'SUYGY0MFVB', 'SW0L8B0LS4', 'T3BB8UT9KP', 'T653MJNWUY', 'T8IWENLDFE', 'TEEDQQPZOI', 'TFJ3UKACLL', 'TG4OSS40W2', 'TGCXP2RE55', 'TTMBZYSHU8', 'TTN1Q2AA35', 'TTN1Q2AA35', 'TTUKMQJUXT', 'TX7XF0XMBL', 'TXJ9VXE4QC', 'TYUGO1UOO3', 'U0JP7MO0CU', 'U1KHQ8NI5I', 'U3YN74R7B9', 'U6SEJCBUKI', 'U86PUE8URG', 'UB8PDQ9GI1', 'UD8WP9Y59H', 'UFKCGRWOE6', 'ULH5Q1UQZ9', 'UMZ44YQOYC', 'UMZ44YQOYC', 'UOHLJNY3WH', 'UOWQ0WV182', 'UYBI9N0WBE', 'V3V21CTLZR', 'VE1W0KD4IA', 'VF0HIW4X63', 'VI0WMPQNAD', 'VIMLBEI2ZQ', 'VKOQ16NG0M', 'VKOQ16NG0M', 'VMKJ7V6SRD', 'VMKJ7V6SRD', 'VQS7PT887V', 'VV8VUVF9XA', 'W3W1ROPU95', 'W6XNV1LQB1', 'W7UAIA7T0S', 'W80YCK8SOZ', 'WCYEL3APQW', 'WIYGVFIKGG', 'WK78CI5LDO', 'WMJ41DPJ3Y', 'WPYHWCFES9', 'WR6V5879WH', 'WXJA24FUA7', 'WYMINM1LOF', 'X0W1ONC8B2', 'X1MS489ILB', 'X2AGMD6P59', 'X4HHLWE7SC', 'X4JDSBA0MD', 'X4JDSBA0MD', 'X5A3YUL1OI', 'X5WMJ9MUOV', 'X7566XEPLR', 'X8NG9V0HJ6', 'XASH2D2L7T', 'XG77M6C95Z', 'XLDFSI23XN', 'XTIR5SW5E8', 'Y0S85J0AZR', 'Y1DJA21YHZ', 'Y4HVM00UR9', 'YBIJBXII7J', 'YF6JV0C9MC', 'YP49T16HM1', 'YV85IDK7KQ', 'YY88X645NE', 'Z0075688AQ', 'Z0RE2PT1R9', 'Z1WO6PGN3E', 'Z9478QITNI', 'ZA6WQHDGBM', 'ZBG7MKFOZ1', 'ZC60EAOUMH', 'ZE7RW2IJMR', 'ZHRR08KKVM', 'ZIN73MJOAB', 'ZINL8265JQ', 'ZKG95SXJUU', 'ZKZ02PIDIU', 'ZM2Y4LL5QP', 'ZQRZQA3DZT', 'ZW84KJ7LJA'])]
    dfm.to_csv('Medications.txt', sep='\t', index=False)
##Collapse medications to 2 types: ICI, Platinum, and Drug Conjugates - merge with patient file

dfm = dfm[['PATIENT_ID', 'AGENT(S)']]
dfmICI = dfm[dfm['AGENT(S)'].isin(['Atezolizumab', 'Nivolumab', 'Pembrolizumab', 'Durvalumab', 'Avelumab'])]
dfmPlat = dfm[dfm['AGENT(S)'].isin(['Cisplatin', 'Carboplatin', 'Oxaliplatin'])]
dfmConjug = dfm[dfm['AGENT(S)'].isin(['Enfortumab Vedotin', 'Sacituzumab Govitecan'])]

dfm.drop_duplicates(subset=['PATIENT_ID', 'AGENT(S)'], keep='first', inplace=True)
dfmPlat.drop_duplicates(subset=['PATIENT_ID', 'AGENT(S)'], keep='first', inplace=True)
dfmConjug.drop_duplicates(subset=['PATIENT_ID', 'AGENT(S)'], keep='first', inplace=True)

dfmA = dfmICI.loc[dfmICI['AGENT(S)'].isin(['Atezolizumab'])]
dfmA.rename(columns={'AGENT(S)': 'ATEZOLIZUMAB'}, inplace=True)
dfmAv = dfmICI.loc[dfmICI['AGENT(S)'].isin(['Avelumab'])]
dfmAv.rename(columns={'AGENT(S)': 'AVELUMAB'}, inplace=True)
dfmN = dfmICI.loc[dfmICI['AGENT(S)'].isin(['Nivolumab'])]
dfmN.rename(columns={'AGENT(S)': 'NIVOLUMAB'}, inplace=True)
dfmP = dfmICI.loc[dfmICI['AGENT(S)'].isin(['Pembrolizumab'])]
dfmP.rename(columns={'AGENT(S)': 'PEMBROLIZUMAB'}, inplace=True)

dfmCis = dfmPlat.loc[dfmPlat['AGENT(S)'].isin(['Cisplatin'])]
dfmCis.rename(columns={'AGENT(S)': 'CISPLATIN'}, inplace=True)
dfmCar = dfmPlat.loc[dfmPlat['AGENT(S)'].isin(['Cisplatin'])]
dfmCar.rename(columns={'AGENT(S)': 'CISPLATIN'}, inplace=True)
dfmCis = dfmPlat.loc[dfmPlat['AGENT(S)'].isin(['Carboplatin'])]
dfmCis.rename(columns={'AGENT(S)': 'CARBOPLATIN'}, inplace=True)
dfmOxa = dfmPlat.loc[dfmPlat['AGENT(S)'].isin(['Oxaliplatin'])]
dfmOxa.rename(columns={'AGENT(S)': 'OXALIPLATIN'}, inplace=True)

dfmEn = dfmConjug.loc[dfmConjug['AGENT(S)'].isin(['Enfortumab Vedotin'])]
dfmEn.rename(columns={'AGENT(S)': 'Enfortumab Vedotin'}, inplace=True)
dfmSa = dfmConjug.loc[dfmConjug['AGENT(S)'].isin(['Sacituzumab Govitecan'])]
dfmSa.rename(columns={'AGENT(S)': 'Sacituzumab Govitecan'}, inplace=True)

dfm1 = dfmA.merge(dfmN, on='PATIENT_ID', how='outer').merge(dfmP, on='PATIENT_ID', how='outer').merge(dfmAv, on='PATIENT_ID', how='outer')
dfm1['IMMUNOTHERAPY'] = 'Yes'

dfm2 = dfmCis.merge(dfmCar, on='PATIENT_ID', how='outer').merge(dfmOxa, on='PATIENT_ID', how='outer')
dfm2['PLATINUM TREATMENT'] = 'Yes'

dfm3 = dfmEn.merge(dfmSa, on='PATIENT_ID', how='outer')
dfm3['CONJUGATE TREATMENT'] = 'Yes'

dfm1 = dfm1.sort_values(by=['PATIENT_ID'])
dfm2 = dfm2.sort_values(by=['PATIENT_ID'])
dfm3 = dfm3.sort_values(by=['PATIENT_ID'])
#Patient Filter
#dfm1 = dfm1[dfm1['PATIENT_ID'].isin(['03GCDX0K5N', '06B7Y294XL', '07EF826VTZ', '0C6NW2TC07', '0G5X064EB6', '0GFLCECW75', '0MJYHKZ6AK', '0RA5URUMVK', '0RFVUNS0M6', '0XBPD7NGNT', '15B16LU7OM', '1H50O85Z25', '1HJFO6EFWS', '1IQTFPMLQR', '1K0JG7QIX8', '1K0JG7QIX8', '1Q58ISCSYD', '1RQGB4W7T4', '1YUCT7U163', '274XK9OLE5', '282INULRAK', '28MTFJ7TIV', '2C9ADE4KW6', '2HNKJUZDQR', '2JRJWC9O7M', '2L2ZRV1LZY', '2N04GOV5XK', '2R23NCXQE0', '2UN4T04H08', '37SHF80OE0', '3C2HGYMBKU', '3E15DF0JVG', '3TLQLPNZP8', '41PP6CPRV7', '45PTETBCPX', '45VCVC3MKN', '4AFTCY1HP0', '4B7HJMVOOV', '4FLIDTYACP', '4FLIDTYACP', '4FTJM67GLA', '4OGM5KMQX2', '4S1A75B0CZ', '4UOOLFGTJC', '4UY5RRTZDA', '506HNR0L9Y', '56QVAODIUI', '5CE11AIS9P', '5JWOSJEGEQ', '5UUQPS8Y5W', '5VD97K8TUD', '63JFZJ9HYR', '64F1NAL85E', '66LL4MEA0D', '69NEPW50R9', '6A9ORBLX71', '6DRRKFWIL1', '6JD9VGBWZB', '6PX4GRURXV', '6T5GP2FH7P', '6U163EB0I2', '70HLF9IQQQ', '70P236XSX2', '731EJZUG4W', '73BBJCCFCG', '75HWMQ9THU', '75QKEG7QT5', '761ATF11EM', '78SU3ENGQ3', '78SU3ENGQ3', '797KYKKNTB', '7DQNA9XN5I', '7JT8CIRL35', '7O8FXBKCMP', '7P2075M6NA', '7XZV89ZKWL', '80II7FE25J', '80TN2OH9U2', '8A9JJV56BE', '8C9KM403EG', '8I52ZJY44W', '8IPGJJGS93', '8K3F2U4AEG', '8SSG3I43TG', '8SWRWW6ZR8', '8U1UG57FYC', '8XOABM89XD', '8YEVOS4E0I', '8ZVLI42YTE', '913EQAJLBI', '923EC0J63U', '928IO9LH0Z', '960N6W495Z', '988C2SHC0G', '995JAM1PGX', '9BZRRRQBLI', '9DFUAVXQTU', '9FUJS5X3V0', '9H5WMVOS0K', '9K443LZVX9', '9RQJ8JFILJ', '9SKYB6E4B4', '9ZXCC6ULKK', '9ZXCC6ULKK', 'A049JYM87K', 'A3AM9T7VPE', 'ACLOR9FFS5', 'ALL5VBRLOM', 'AP93ZGZX2W', 'AVMUPL2GS0', 'B3D07F2E2S', 'BCZK5PTSG9', 'BK6DKQTSI5', 'BTHE1B9DZ2', 'BU8GCS3MYT', 'BV2N69FYI6', 'BV2N69FYI6', 'BXXKG0HZIK', 'C5A7KUNVCS', 'CADAN6NOJ4', 'CGDPFH8IBW', 'CK9QL0A5ET', 'CM03MT997O', 'CM03MT997O', 'CM03MT997O', 'CTGHR4UZ18', 'CW0ZVIVU0P', 'CZE113UZIR', 'D0WP0B97Q7', 'D2C5F4UMXY', 'D7KAPTXFOG', 'DCPKIAYXDB', 'DQF8CHSKPG', 'DYKRF4KH3H', 'E2H89HF5KZ', 'E4H964SOU9', 'EF7AOZI79O', 'EHOVA8D07G', 'EKSOBYXBVG', 'EKXLPR0080', 'EMVBZPMWST', 'EOZ9LELKBM', 'EQNGL8N9W7', 'EU4WQTA1KL', 'F1XUR9R1RK', 'F2TQDDM3HH', 'F905RYZFIO', 'FNVA6L545Q', 'FOTCBZ3OMJ', 'FPLKDOPOH3', 'FVZV8LMNI5', 'FXGR65TB4V', 'FXSR99KNNT', 'GA7GNCV7KM', 'GF6S6LZVF5', 'GH6KKNEPWX', 'GIJZDVWC2Q', 'GNJVFE3HRM', 'GNXLDBYWN1', 'GOYPZUZNAW', 'GPLOBS1HRX', 'GU9WE974XH', 'GWNOJYMTTB', 'HG7EO5096W', 'HUCN0XE91R', 'HWECMUCW0Z', 'HWIG0FSC6Z', 'HXKS3QWP6F', 'HY2ZUQI0SD', 'I2EM49R0AO', 'I65GB8ZEJ2', 'IMH5IF30JJ', 'IW3OCCZBPJ', 'IW7K4RUJG9', 'IXEPJTR5B9', 'J0D6I815XU', 'J15Q9C8MIE', 'J3C9QS0P78', 'J60EKEW4OF', 'J81JUVHO9L', 'J9NHJFS8IN', 'JBLZPPZSFH', 'JG2Z3K6UC5', 'JHBPEVZ48F', 'JHZSJNIKIP', 'JJGV4LPENY', 'JLGA08C3K5', 'JM3MBDWG5H', 'JNXYDUFEED', 'JPSK03CCDH', 'JPSK03CCDH', 'JQHVTFHS1X', 'JRC4RTZNTK', 'JSWZX87JYV', 'JZKKZ1OTQJ', 'K1DVC6O9C7', 'K60LO6G1FG', 'K8DFO8FTT1', 'K8OFSLEHIX', 'KAN1H1VWQG', 'KC8HPZPX7T', 'KPCHILTX70', 'KSTUVZCJIF', 'KW650F0557', 'L410MYBLDN', 'L6RMSU0BXG', 'L8HCOD7F64', 'LATEJ509E1', 'LLXVLRWEAE', 'LM4ZZPQQFO', 'LS6UIWKFEJ', 'LX5K9JZKLH', 'LXELBL78SE', 'M1D0DA5RWT', 'M9G3Y6A4X7', 'MG1SUQR5K8', 'MIVJ67B9F5', 'MJYLGG6W0K', 'MRZW12GMST', 'MYG3XG6SAQ', 'N4K1N7UUMM', 'N675BI25Y6', 'N8USV2HKJ6', 'NAMDA2GZXA', 'NAR8MDP01R', 'NDV4VM27O6', 'NE17FOVUC5', 'NK65MD354P', 'NOD926GZ41', 'NPUNSPR11V', 'O8RFN0K0CC', 'OA1K7SKWYW', 'OFMS9RN5W8', 'OFOSHGCQE4', 'OJQ5IBJCDB', 'ONHYTARPID', 'OPI88VIT01', 'OQYBJZXP4D', 'ORL6LF2IMD', 'OWNR0K9J5O', 'OXKYPBVO0Q', 'OXKYPBVO0Q', 'P0YOHVOGGY', 'P6F7GJOYOS', 'PCM0C3L6G1', 'PGVXTE18CS', 'PJ316U0AO6', 'PQDMEIJU8A', 'PTN0DBFVOH', 'Q1LFWVL6GK', 'Q1TM0Z8WWH', 'Q3XX8RO9CY', 'QB511DI6QC', 'QD6U65T2FH', 'QDTQI9DG2U', 'QJ80IGAAAA', 'QO27BH1WTV', 'QUEJJTOM01', 'QVEQP2FBPS', 'QVLKW6YJ9S', 'QZ68W571DZ', 'R0JJZPXHE6', 'R9EFJQNQD7', 'RC5DO9WG9A', 'RDNMW7PS5V', 'RL11WPPBNV', 'RN9YUHZDTE', 'RRWTZZSW3F', 'S414XWRR6V', 'S5AEWMYN2T', 'S5R5RW4L05', 'S7GI3EL7N1', 'S9CLR0UMQD', 'S9GIV41A0F', 'SF953NJFYV', 'SQV8WBN78P', 'STZJJQ7LKL', 'SUYGY0MFVB', 'SW0L8B0LS4', 'T3BB8UT9KP', 'T653MJNWUY', 'T8IWENLDFE', 'TEEDQQPZOI', 'TFJ3UKACLL', 'TG4OSS40W2', 'TGCXP2RE55', 'TTMBZYSHU8', 'TTN1Q2AA35', 'TTN1Q2AA35', 'TTUKMQJUXT', 'TX7XF0XMBL', 'TXJ9VXE4QC', 'TYUGO1UOO3', 'U0JP7MO0CU', 'U1KHQ8NI5I', 'U3YN74R7B9', 'U6SEJCBUKI', 'U86PUE8URG', 'UB8PDQ9GI1', 'UD8WP9Y59H', 'UFKCGRWOE6', 'ULH5Q1UQZ9', 'UMZ44YQOYC', 'UMZ44YQOYC', 'UOHLJNY3WH', 'UOWQ0WV182', 'UYBI9N0WBE', 'V3V21CTLZR', 'VE1W0KD4IA', 'VF0HIW4X63', 'VI0WMPQNAD', 'VIMLBEI2ZQ', 'VKOQ16NG0M', 'VKOQ16NG0M', 'VMKJ7V6SRD', 'VMKJ7V6SRD', 'VQS7PT887V', 'VV8VUVF9XA', 'W3W1ROPU95', 'W6XNV1LQB1', 'W7UAIA7T0S', 'W80YCK8SOZ', 'WCYEL3APQW', 'WIYGVFIKGG', 'WK78CI5LDO', 'WMJ41DPJ3Y', 'WPYHWCFES9', 'WR6V5879WH', 'WXJA24FUA7', 'WYMINM1LOF', 'X0W1ONC8B2', 'X1MS489ILB', 'X2AGMD6P59', 'X4HHLWE7SC', 'X4JDSBA0MD', 'X4JDSBA0MD', 'X5A3YUL1OI', 'X5WMJ9MUOV', 'X7566XEPLR', 'X8NG9V0HJ6', 'XASH2D2L7T', 'XG77M6C95Z', 'XLDFSI23XN', 'XTIR5SW5E8', 'Y0S85J0AZR', 'Y1DJA21YHZ', 'Y4HVM00UR9', 'YBIJBXII7J', 'YF6JV0C9MC', 'YP49T16HM1', 'YV85IDK7KQ', 'YY88X645NE', 'Z0075688AQ', 'Z0RE2PT1R9', 'Z1WO6PGN3E', 'Z9478QITNI', 'ZA6WQHDGBM', 'ZBG7MKFOZ1', 'ZC60EAOUMH', 'ZE7RW2IJMR', 'ZHRR08KKVM', 'ZIN73MJOAB', 'ZINL8265JQ', 'ZKG95SXJUU', 'ZKZ02PIDIU', 'ZM2Y4LL5QP', 'ZQRZQA3DZT', 'ZW84KJ7LJA'])]
PatientCDE = PatientCDE.merge(dfm1, how='left', on='PATIENT_ID').merge(dfm2, how='left', on='PATIENT_ID').merge(dfm3, how='left', on='PATIENT_ID')
dfm1.to_csv('Medications2.txt', sep='\t', index=False)
dfm2.to_csv('Medications_PLATINUM.txt', sep='\t', index=False)
dfm3.to_csv('Medications_CONJUGATES.txt', sep='\t', index=False)
PatientCDE.to_csv('PatientCDE.txt', sep='\t', index=False)

#ADD Patient Weight, BMI
tBMI = []
FileSample1 = glob.glob('*PhysicalAssessment_V4.csv')
for f in FileSample1: 
    dfstart2 = pd.read_csv(f, index_col=None, header=0)
    dfstart2.rename(columns={'AvatarKey': 'PATIENT_ID'}, inplace=True)
    tBMI.append(dfstart2)
    dfBMI = pd.concat(tBMI, axis=0, ignore_index=True)
    ##Strip all Unknown/Not Applicable
    dfBMI = dfBMI.replace(regex={r'Unknown/Not Applicable': ''})
    dfBMI = dfBMI.replace(regex={r'Age 90 or older': '91'})
    #Convert headers to snake_case and upper case
    dfBMI.columns = [camel_to_snake(col) for col in dfBMI.columns]
    dfBMI.rename(columns = str.upper , inplace = True)
    dfBMI.rename(columns={'P_A_T_I_E_N_T__I_D':'PATIENT_ID', 'B_M_I': 'BMI'}, inplace=True)
    #dfBMI['AverageBMI'] = dfBMI['BMI']
    dfBMI.to_csv('DiffBMI.txt', sep='\t', index=False)
    dfBMI['BMI'] = pd.to_numeric(dfBMI['BMI'], errors='coerce')
    dfBMI['BODY_WEIGHT'] = pd.to_numeric(dfBMI['BODY_WEIGHT'], errors='coerce')
    AverageBMI = dfBMI
    DiffBMI = dfBMI
    AverageBMI['BMI'] = AverageBMI.groupby('PATIENT_ID')['BMI'].transform('mean')
    AverageBMI['BODY_WEIGHT'] = AverageBMI.groupby('PATIENT_ID')['BODY_WEIGHT'].transform('mean')
    AverageBMI.drop_duplicates(subset=['PATIENT_ID'], keep='first', inplace=True)
    AverageBMI = AverageBMI[['PATIENT_ID', 'BODY_WEIGHT', 'BMI']]
    AverageBMI.rename(columns={'BODY_WEIGHT': 'AVG_BODY_WEIGHT', 'BMI': 'AVG_BMI'}, inplace=True)
    AverageBMI.to_csv('AverageBMI.txt', sep='\t', index=False)
    #DiffBMI['BMI'] = DiffBMI.groupby(['PATIENT_ID'])['BMI'].transform(lambda x: x.diff())
    #DiffBMI['DELTA_BMI'] = DiffBMI.groupby('PATIENT_ID').apply(lambda x: x['BMI'].diff()).values
    DiagnosisDate = pd.read_csv('SampleCDEdata.txt', sep='\t', header=0)
    DiagnosisDate =DiagnosisDate[['PATIENT_ID', 'AGE_AT_DIAGNOSIS']]
    DiffBMI = DiffBMI.merge(DiagnosisDate, how='left', on='PATIENT_ID')
    DiffBMI.to_csv('DiffBMI.txt', sep='\t', index=False)
    PatientCDE = PatientCDE.merge(dfBMI, how='left', on='PATIENT_ID').merge(AverageBMI, how='left', on='PATIENT_ID')
    PatientCDE.rename(columns={'BMI':'BMI_AT_DIAGNOSIS', 'BODY_WEIGHT': 'BODY_WEIGHT_AT_DIAGNOSIS'}, inplace=True)
    PatientCDE.to_csv('PatientCDE.txt', sep='\t', index=False)
print('\n\n\n######################################################## COMPLETE PATIENT HEADERS  ###################################################\n\n')

print('\n\n\n######################################################## STARTING PFS MEDS EXTRACTION  ###################################################\n\n')
##Use change in medication to add PFS data to outcomes
t9 = []
FilePatient6 = glob.glob('*_Medications_V4.csv')
for f in FilePatient6:
    dfstart9 = pd.read_csv(f, index_col=None, header=0)
    t9.append(dfstart9)
    pfsMed2 = pd.concat(t9, axis=0, ignore_index=True)
#FileSample2 = glob.glob('*_Medications_V4.csv')
#for f in FileSample2:
#    pfsMed2 = pd.read_csv(f)    
    pfsMed2 = pfsMed2[['AvatarKey', 'MedicationInd', 'MedReasonNoneGiven', 'MedPrimaryDiagnosisSiteCode', 'MedPrimaryDiagnosisSite', 'Medication', 'MedLineRegimen', 'AgeAtMedStart', 'AgeAtMedStartFlag', 'YearOfMedStart', 'MedContinuing', 'AgeAtMedStop', 'AgeAtMedStopFlag', 'YearOfMedStop', 'ChangeOfTreatment']]
    print('\n\n************************************************************ Medication Change from Progression ********************************\n\n')
    pfsMed2.columns = [camel_to_snake(col) for col in pfsMed2.columns]
    pfsMed2.rename(columns = str.upper , inplace = True)
    pfsMed2.rename(columns={'AVATAR_KEY': 'PATIENT_ID'}, inplace=True)
    #Patient Filter
    #pfsMed2 = pfsMed2[pfsMed2['PATIENT_ID'].isin(['03GCDX0K5N', '06B7Y294XL', '07EF826VTZ', '0C6NW2TC07', '0G5X064EB6', '0GFLCECW75', '0MJYHKZ6AK', '0RA5URUMVK', '0RFVUNS0M6', '0XBPD7NGNT', '15B16LU7OM', '1H50O85Z25', '1HJFO6EFWS', '1IQTFPMLQR', '1K0JG7QIX8', '1K0JG7QIX8', '1Q58ISCSYD', '1RQGB4W7T4', '1YUCT7U163', '274XK9OLE5', '282INULRAK', '28MTFJ7TIV', '2C9ADE4KW6', '2HNKJUZDQR', '2JRJWC9O7M', '2L2ZRV1LZY', '2N04GOV5XK', '2R23NCXQE0', '2UN4T04H08', '37SHF80OE0', '3C2HGYMBKU', '3E15DF0JVG', '3TLQLPNZP8', '41PP6CPRV7', '45PTETBCPX', '45VCVC3MKN', '4AFTCY1HP0', '4B7HJMVOOV', '4FLIDTYACP', '4FLIDTYACP', '4FTJM67GLA', '4OGM5KMQX2', '4S1A75B0CZ', '4UOOLFGTJC', '4UY5RRTZDA', '506HNR0L9Y', '56QVAODIUI', '5CE11AIS9P', '5JWOSJEGEQ', '5UUQPS8Y5W', '5VD97K8TUD', '63JFZJ9HYR', '64F1NAL85E', '66LL4MEA0D', '69NEPW50R9', '6A9ORBLX71', '6DRRKFWIL1', '6JD9VGBWZB', '6PX4GRURXV', '6T5GP2FH7P', '6U163EB0I2', '70HLF9IQQQ', '70P236XSX2', '731EJZUG4W', '73BBJCCFCG', '75HWMQ9THU', '75QKEG7QT5', '761ATF11EM', '78SU3ENGQ3', '78SU3ENGQ3', '797KYKKNTB', '7DQNA9XN5I', '7JT8CIRL35', '7O8FXBKCMP', '7P2075M6NA', '7XZV89ZKWL', '80II7FE25J', '80TN2OH9U2', '8A9JJV56BE', '8C9KM403EG', '8I52ZJY44W', '8IPGJJGS93', '8K3F2U4AEG', '8SSG3I43TG', '8SWRWW6ZR8', '8U1UG57FYC', '8XOABM89XD', '8YEVOS4E0I', '8ZVLI42YTE', '913EQAJLBI', '923EC0J63U', '928IO9LH0Z', '960N6W495Z', '988C2SHC0G', '995JAM1PGX', '9BZRRRQBLI', '9DFUAVXQTU', '9FUJS5X3V0', '9H5WMVOS0K', '9K443LZVX9', '9RQJ8JFILJ', '9SKYB6E4B4', '9ZXCC6ULKK', '9ZXCC6ULKK', 'A049JYM87K', 'A3AM9T7VPE', 'ACLOR9FFS5', 'ALL5VBRLOM', 'AP93ZGZX2W', 'AVMUPL2GS0', 'B3D07F2E2S', 'BCZK5PTSG9', 'BK6DKQTSI5', 'BTHE1B9DZ2', 'BU8GCS3MYT', 'BV2N69FYI6', 'BV2N69FYI6', 'BXXKG0HZIK', 'C5A7KUNVCS', 'CADAN6NOJ4', 'CGDPFH8IBW', 'CK9QL0A5ET', 'CM03MT997O', 'CM03MT997O', 'CM03MT997O', 'CTGHR4UZ18', 'CW0ZVIVU0P', 'CZE113UZIR', 'D0WP0B97Q7', 'D2C5F4UMXY', 'D7KAPTXFOG', 'DCPKIAYXDB', 'DQF8CHSKPG', 'DYKRF4KH3H', 'E2H89HF5KZ', 'E4H964SOU9', 'EF7AOZI79O', 'EHOVA8D07G', 'EKSOBYXBVG', 'EKXLPR0080', 'EMVBZPMWST', 'EOZ9LELKBM', 'EQNGL8N9W7', 'EU4WQTA1KL', 'F1XUR9R1RK', 'F2TQDDM3HH', 'F905RYZFIO', 'FNVA6L545Q', 'FOTCBZ3OMJ', 'FPLKDOPOH3', 'FVZV8LMNI5', 'FXGR65TB4V', 'FXSR99KNNT', 'GA7GNCV7KM', 'GF6S6LZVF5', 'GH6KKNEPWX', 'GIJZDVWC2Q', 'GNJVFE3HRM', 'GNXLDBYWN1', 'GOYPZUZNAW', 'GPLOBS1HRX', 'GU9WE974XH', 'GWNOJYMTTB', 'HG7EO5096W', 'HUCN0XE91R', 'HWECMUCW0Z', 'HWIG0FSC6Z', 'HXKS3QWP6F', 'HY2ZUQI0SD', 'I2EM49R0AO', 'I65GB8ZEJ2', 'IMH5IF30JJ', 'IW3OCCZBPJ', 'IW7K4RUJG9', 'IXEPJTR5B9', 'J0D6I815XU', 'J15Q9C8MIE', 'J3C9QS0P78', 'J60EKEW4OF', 'J81JUVHO9L', 'J9NHJFS8IN', 'JBLZPPZSFH', 'JG2Z3K6UC5', 'JHBPEVZ48F', 'JHZSJNIKIP', 'JJGV4LPENY', 'JLGA08C3K5', 'JM3MBDWG5H', 'JNXYDUFEED', 'JPSK03CCDH', 'JPSK03CCDH', 'JQHVTFHS1X', 'JRC4RTZNTK', 'JSWZX87JYV', 'JZKKZ1OTQJ', 'K1DVC6O9C7', 'K60LO6G1FG', 'K8DFO8FTT1', 'K8OFSLEHIX', 'KAN1H1VWQG', 'KC8HPZPX7T', 'KPCHILTX70', 'KSTUVZCJIF', 'KW650F0557', 'L410MYBLDN', 'L6RMSU0BXG', 'L8HCOD7F64', 'LATEJ509E1', 'LLXVLRWEAE', 'LM4ZZPQQFO', 'LS6UIWKFEJ', 'LX5K9JZKLH', 'LXELBL78SE', 'M1D0DA5RWT', 'M9G3Y6A4X7', 'MG1SUQR5K8', 'MIVJ67B9F5', 'MJYLGG6W0K', 'MRZW12GMST', 'MYG3XG6SAQ', 'N4K1N7UUMM', 'N675BI25Y6', 'N8USV2HKJ6', 'NAMDA2GZXA', 'NAR8MDP01R', 'NDV4VM27O6', 'NE17FOVUC5', 'NK65MD354P', 'NOD926GZ41', 'NPUNSPR11V', 'O8RFN0K0CC', 'OA1K7SKWYW', 'OFMS9RN5W8', 'OFOSHGCQE4', 'OJQ5IBJCDB', 'ONHYTARPID', 'OPI88VIT01', 'OQYBJZXP4D', 'ORL6LF2IMD', 'OWNR0K9J5O', 'OXKYPBVO0Q', 'OXKYPBVO0Q', 'P0YOHVOGGY', 'P6F7GJOYOS', 'PCM0C3L6G1', 'PGVXTE18CS', 'PJ316U0AO6', 'PQDMEIJU8A', 'PTN0DBFVOH', 'Q1LFWVL6GK', 'Q1TM0Z8WWH', 'Q3XX8RO9CY', 'QB511DI6QC', 'QD6U65T2FH', 'QDTQI9DG2U', 'QJ80IGAAAA', 'QO27BH1WTV', 'QUEJJTOM01', 'QVEQP2FBPS', 'QVLKW6YJ9S', 'QZ68W571DZ', 'R0JJZPXHE6', 'R9EFJQNQD7', 'RC5DO9WG9A', 'RDNMW7PS5V', 'RL11WPPBNV', 'RN9YUHZDTE', 'RRWTZZSW3F', 'S414XWRR6V', 'S5AEWMYN2T', 'S5R5RW4L05', 'S7GI3EL7N1', 'S9CLR0UMQD', 'S9GIV41A0F', 'SF953NJFYV', 'SQV8WBN78P', 'STZJJQ7LKL', 'SUYGY0MFVB', 'SW0L8B0LS4', 'T3BB8UT9KP', 'T653MJNWUY', 'T8IWENLDFE', 'TEEDQQPZOI', 'TFJ3UKACLL', 'TG4OSS40W2', 'TGCXP2RE55', 'TTMBZYSHU8', 'TTN1Q2AA35', 'TTN1Q2AA35', 'TTUKMQJUXT', 'TX7XF0XMBL', 'TXJ9VXE4QC', 'TYUGO1UOO3', 'U0JP7MO0CU', 'U1KHQ8NI5I', 'U3YN74R7B9', 'U6SEJCBUKI', 'U86PUE8URG', 'UB8PDQ9GI1', 'UD8WP9Y59H', 'UFKCGRWOE6', 'ULH5Q1UQZ9', 'UMZ44YQOYC', 'UMZ44YQOYC', 'UOHLJNY3WH', 'UOWQ0WV182', 'UYBI9N0WBE', 'V3V21CTLZR', 'VE1W0KD4IA', 'VF0HIW4X63', 'VI0WMPQNAD', 'VIMLBEI2ZQ', 'VKOQ16NG0M', 'VKOQ16NG0M', 'VMKJ7V6SRD', 'VMKJ7V6SRD', 'VQS7PT887V', 'VV8VUVF9XA', 'W3W1ROPU95', 'W6XNV1LQB1', 'W7UAIA7T0S', 'W80YCK8SOZ', 'WCYEL3APQW', 'WIYGVFIKGG', 'WK78CI5LDO', 'WMJ41DPJ3Y', 'WPYHWCFES9', 'WR6V5879WH', 'WXJA24FUA7', 'WYMINM1LOF', 'X0W1ONC8B2', 'X1MS489ILB', 'X2AGMD6P59', 'X4HHLWE7SC', 'X4JDSBA0MD', 'X4JDSBA0MD', 'X5A3YUL1OI', 'X5WMJ9MUOV', 'X7566XEPLR', 'X8NG9V0HJ6', 'XASH2D2L7T', 'XG77M6C95Z', 'XLDFSI23XN', 'XTIR5SW5E8', 'Y0S85J0AZR', 'Y1DJA21YHZ', 'Y4HVM00UR9', 'YBIJBXII7J', 'YF6JV0C9MC', 'YP49T16HM1', 'YV85IDK7KQ', 'YY88X645NE', 'Z0075688AQ', 'Z0RE2PT1R9', 'Z1WO6PGN3E', 'Z9478QITNI', 'ZA6WQHDGBM', 'ZBG7MKFOZ1', 'ZC60EAOUMH', 'ZE7RW2IJMR', 'ZHRR08KKVM', 'ZIN73MJOAB', 'ZINL8265JQ', 'ZKG95SXJUU', 'ZKZ02PIDIU', 'ZM2Y4LL5QP', 'ZQRZQA3DZT', 'ZW84KJ7LJA'])]

    pfsMed2 = pfsMed2[['PATIENT_ID', 'MEDICATION', 'MED_LINE_REGIMEN', 'AGE_AT_MED_START', 'AGE_AT_MED_STOP', 'CHANGE_OF_TREATMENT']]
   ####################################### 
    #BUILD 3 MED CLASS DATAFRAMES FOR PFS
    pfsMed_plat = pfsMed2[pfsMed2['MEDICATION'].isin(['Cisplatin', 'Carboplatin', 'Oxaliplatin'])]
    pfsMed_conjug = pfsMed2[pfsMed2['MEDICATION'].isin(['Enfortumab Vedotin', 'Sacituzumab Govitecan'])]
    pfsMed_ICI = pfsMed2[pfsMed2['MEDICATION'].isin(['Atezolizumab', 'Avelumab', 'Nivolumab', 'Pembrolizumab', 'Durvalumab'])]
    
    #keep first ICI medication
    pfsMed_ICI.sort_values(by=['AGE_AT_MED_START'], inplace = True)
    pfsMed_ICI.drop_duplicates(subset=['PATIENT_ID'], keep='first',  inplace=True)
    #Rename columns for PFS status
    pfsMed_ICI.rename(columns={'CHANGE_OF_TREATMENT': 'PFS_STATUS'}, inplace=True)
    pfsMed_ICI['AGE_AT_MED_START'] = pd.to_numeric(pfsMed_ICI.AGE_AT_MED_START, errors='coerce')
    pfsMed_ICI['AGE_AT_MED_STOP'] = pd.to_numeric(pfsMed_ICI.AGE_AT_MED_STOP, errors='coerce')
    pfsMed_ICI['PFS_MONTHS'] = (pfsMed_ICI['AGE_AT_MED_STOP'] - pfsMed_ICI['AGE_AT_MED_START'])*12
    pfsMed_ICI['MEDICATION_INDICATION'] = 'Medication(s) Change'
    pfsMed_ICI.rename(columns={'AGE_AT_MED_START': 'AGE_AT_MED_START_Med', 'AGE_AT_MED_STOP': 'STOP_DATE_Med', 'MEDICATION': 'AGENT(S)'}, inplace=True)
    
    #Save medical file for other PFS intersection with therapy
    pfsMed_ICI.to_csv('MedicationsPFStempICI.txt', sep='\t', index=False)
    #Build Progression file from change in treatment
    pfsMed_ICI = pfsMed_ICI[pfsMed_ICI['PFS_STATUS'].isin(['Yes, due to progression'])]
    pfsMed_ICI['PFS_STATUS'] = pfsMed_ICI['PFS_STATUS'].str.replace('Yes, due to progression', '1:Progressed')
    #pfsMed_ICI.drop(['MED_LINE_REGIMEN'], axis=1, inplace = True)
    #add source column
    pfsMed_ICI['SOURCE'] = 'MEDICATIONS-CHANGE_OF_TREATMENT'
    pfsMed_ICI = pfsMed_ICI[['PATIENT_ID', 'AGENT(S)', 'PFS_STATUS', 'PFS_MONTHS', 'AGE_AT_MED_START_Med', 'STOP_DATE_Med', 'SOURCE', 'MED_LINE_REGIMEN']]
    pfsMed_ICI.to_csv('PFSmedICI.txt', sep='\t', index=False)
    

print('\n\n\n######################################################## STARTING PFS OUTCOMES EXTRACTION  ###################################################\n\n')

t8 = []
FilePatient5 = glob.glob('*_Outcomes_V4.csv')
for f in FilePatient5:
    dfstart8 = pd.read_csv(f, index_col=None, header=0)
    t8.append(dfstart8)
    pfs = pd.concat(t8, axis=0, ignore_index=True)
#FileSample2 = glob.glob('*_Outcomes_V4.csv')
#for f in FileSample2:
#    pfs = pd.read_csv(f)
    #cleanup unused columns
    pfs = pfs[['AvatarKey', 'ProgRecurInd', 'AgeAtProgRecur', 'AgeAtCurrentDiseaseStatus', 'CurrentDiseaseStatus']]
    #Convert headers to snake_case and upper case
    pfs.columns = [camel_to_snake(col) for col in pfs.columns]
    pfs.rename(columns = str.upper , inplace = True)
    pfs.rename(columns={'AVATAR_KEY': 'PATIENT_ID'}, inplace=True)
    #Patient Filter
    #pfs = pfs[pfs['PATIENT_ID'].isin(['03GCDX0K5N', '06B7Y294XL', '07EF826VTZ', '0C6NW2TC07', '0G5X064EB6', '0GFLCECW75', '0MJYHKZ6AK', '0RA5URUMVK', '0RFVUNS0M6', '0XBPD7NGNT', '15B16LU7OM', '1H50O85Z25', '1HJFO6EFWS', '1IQTFPMLQR', '1K0JG7QIX8', '1K0JG7QIX8', '1Q58ISCSYD', '1RQGB4W7T4', '1YUCT7U163', '274XK9OLE5', '282INULRAK', '28MTFJ7TIV', '2C9ADE4KW6', '2HNKJUZDQR', '2JRJWC9O7M', '2L2ZRV1LZY', '2N04GOV5XK', '2R23NCXQE0', '2UN4T04H08', '37SHF80OE0', '3C2HGYMBKU', '3E15DF0JVG', '3TLQLPNZP8', '41PP6CPRV7', '45PTETBCPX', '45VCVC3MKN', '4AFTCY1HP0', '4B7HJMVOOV', '4FLIDTYACP', '4FLIDTYACP', '4FTJM67GLA', '4OGM5KMQX2', '4S1A75B0CZ', '4UOOLFGTJC', '4UY5RRTZDA', '506HNR0L9Y', '56QVAODIUI', '5CE11AIS9P', '5JWOSJEGEQ', '5UUQPS8Y5W', '5VD97K8TUD', '63JFZJ9HYR', '64F1NAL85E', '66LL4MEA0D', '69NEPW50R9', '6A9ORBLX71', '6DRRKFWIL1', '6JD9VGBWZB', '6PX4GRURXV', '6T5GP2FH7P', '6U163EB0I2', '70HLF9IQQQ', '70P236XSX2', '731EJZUG4W', '73BBJCCFCG', '75HWMQ9THU', '75QKEG7QT5', '761ATF11EM', '78SU3ENGQ3', '78SU3ENGQ3', '797KYKKNTB', '7DQNA9XN5I', '7JT8CIRL35', '7O8FXBKCMP', '7P2075M6NA', '7XZV89ZKWL', '80II7FE25J', '80TN2OH9U2', '8A9JJV56BE', '8C9KM403EG', '8I52ZJY44W', '8IPGJJGS93', '8K3F2U4AEG', '8SSG3I43TG', '8SWRWW6ZR8', '8U1UG57FYC', '8XOABM89XD', '8YEVOS4E0I', '8ZVLI42YTE', '913EQAJLBI', '923EC0J63U', '928IO9LH0Z', '960N6W495Z', '988C2SHC0G', '995JAM1PGX', '9BZRRRQBLI', '9DFUAVXQTU', '9FUJS5X3V0', '9H5WMVOS0K', '9K443LZVX9', '9RQJ8JFILJ', '9SKYB6E4B4', '9ZXCC6ULKK', '9ZXCC6ULKK', 'A049JYM87K', 'A3AM9T7VPE', 'ACLOR9FFS5', 'ALL5VBRLOM', 'AP93ZGZX2W', 'AVMUPL2GS0', 'B3D07F2E2S', 'BCZK5PTSG9', 'BK6DKQTSI5', 'BTHE1B9DZ2', 'BU8GCS3MYT', 'BV2N69FYI6', 'BV2N69FYI6', 'BXXKG0HZIK', 'C5A7KUNVCS', 'CADAN6NOJ4', 'CGDPFH8IBW', 'CK9QL0A5ET', 'CM03MT997O', 'CM03MT997O', 'CM03MT997O', 'CTGHR4UZ18', 'CW0ZVIVU0P', 'CZE113UZIR', 'D0WP0B97Q7', 'D2C5F4UMXY', 'D7KAPTXFOG', 'DCPKIAYXDB', 'DQF8CHSKPG', 'DYKRF4KH3H', 'E2H89HF5KZ', 'E4H964SOU9', 'EF7AOZI79O', 'EHOVA8D07G', 'EKSOBYXBVG', 'EKXLPR0080', 'EMVBZPMWST', 'EOZ9LELKBM', 'EQNGL8N9W7', 'EU4WQTA1KL', 'F1XUR9R1RK', 'F2TQDDM3HH', 'F905RYZFIO', 'FNVA6L545Q', 'FOTCBZ3OMJ', 'FPLKDOPOH3', 'FVZV8LMNI5', 'FXGR65TB4V', 'FXSR99KNNT', 'GA7GNCV7KM', 'GF6S6LZVF5', 'GH6KKNEPWX', 'GIJZDVWC2Q', 'GNJVFE3HRM', 'GNXLDBYWN1', 'GOYPZUZNAW', 'GPLOBS1HRX', 'GU9WE974XH', 'GWNOJYMTTB', 'HG7EO5096W', 'HUCN0XE91R', 'HWECMUCW0Z', 'HWIG0FSC6Z', 'HXKS3QWP6F', 'HY2ZUQI0SD', 'I2EM49R0AO', 'I65GB8ZEJ2', 'IMH5IF30JJ', 'IW3OCCZBPJ', 'IW7K4RUJG9', 'IXEPJTR5B9', 'J0D6I815XU', 'J15Q9C8MIE', 'J3C9QS0P78', 'J60EKEW4OF', 'J81JUVHO9L', 'J9NHJFS8IN', 'JBLZPPZSFH', 'JG2Z3K6UC5', 'JHBPEVZ48F', 'JHZSJNIKIP', 'JJGV4LPENY', 'JLGA08C3K5', 'JM3MBDWG5H', 'JNXYDUFEED', 'JPSK03CCDH', 'JPSK03CCDH', 'JQHVTFHS1X', 'JRC4RTZNTK', 'JSWZX87JYV', 'JZKKZ1OTQJ', 'K1DVC6O9C7', 'K60LO6G1FG', 'K8DFO8FTT1', 'K8OFSLEHIX', 'KAN1H1VWQG', 'KC8HPZPX7T', 'KPCHILTX70', 'KSTUVZCJIF', 'KW650F0557', 'L410MYBLDN', 'L6RMSU0BXG', 'L8HCOD7F64', 'LATEJ509E1', 'LLXVLRWEAE', 'LM4ZZPQQFO', 'LS6UIWKFEJ', 'LX5K9JZKLH', 'LXELBL78SE', 'M1D0DA5RWT', 'M9G3Y6A4X7', 'MG1SUQR5K8', 'MIVJ67B9F5', 'MJYLGG6W0K', 'MRZW12GMST', 'MYG3XG6SAQ', 'N4K1N7UUMM', 'N675BI25Y6', 'N8USV2HKJ6', 'NAMDA2GZXA', 'NAR8MDP01R', 'NDV4VM27O6', 'NE17FOVUC5', 'NK65MD354P', 'NOD926GZ41', 'NPUNSPR11V', 'O8RFN0K0CC', 'OA1K7SKWYW', 'OFMS9RN5W8', 'OFOSHGCQE4', 'OJQ5IBJCDB', 'ONHYTARPID', 'OPI88VIT01', 'OQYBJZXP4D', 'ORL6LF2IMD', 'OWNR0K9J5O', 'OXKYPBVO0Q', 'OXKYPBVO0Q', 'P0YOHVOGGY', 'P6F7GJOYOS', 'PCM0C3L6G1', 'PGVXTE18CS', 'PJ316U0AO6', 'PQDMEIJU8A', 'PTN0DBFVOH', 'Q1LFWVL6GK', 'Q1TM0Z8WWH', 'Q3XX8RO9CY', 'QB511DI6QC', 'QD6U65T2FH', 'QDTQI9DG2U', 'QJ80IGAAAA', 'QO27BH1WTV', 'QUEJJTOM01', 'QVEQP2FBPS', 'QVLKW6YJ9S', 'QZ68W571DZ', 'R0JJZPXHE6', 'R9EFJQNQD7', 'RC5DO9WG9A', 'RDNMW7PS5V', 'RL11WPPBNV', 'RN9YUHZDTE', 'RRWTZZSW3F', 'S414XWRR6V', 'S5AEWMYN2T', 'S5R5RW4L05', 'S7GI3EL7N1', 'S9CLR0UMQD', 'S9GIV41A0F', 'SF953NJFYV', 'SQV8WBN78P', 'STZJJQ7LKL', 'SUYGY0MFVB', 'SW0L8B0LS4', 'T3BB8UT9KP', 'T653MJNWUY', 'T8IWENLDFE', 'TEEDQQPZOI', 'TFJ3UKACLL', 'TG4OSS40W2', 'TGCXP2RE55', 'TTMBZYSHU8', 'TTN1Q2AA35', 'TTN1Q2AA35', 'TTUKMQJUXT', 'TX7XF0XMBL', 'TXJ9VXE4QC', 'TYUGO1UOO3', 'U0JP7MO0CU', 'U1KHQ8NI5I', 'U3YN74R7B9', 'U6SEJCBUKI', 'U86PUE8URG', 'UB8PDQ9GI1', 'UD8WP9Y59H', 'UFKCGRWOE6', 'ULH5Q1UQZ9', 'UMZ44YQOYC', 'UMZ44YQOYC', 'UOHLJNY3WH', 'UOWQ0WV182', 'UYBI9N0WBE', 'V3V21CTLZR', 'VE1W0KD4IA', 'VF0HIW4X63', 'VI0WMPQNAD', 'VIMLBEI2ZQ', 'VKOQ16NG0M', 'VKOQ16NG0M', 'VMKJ7V6SRD', 'VMKJ7V6SRD', 'VQS7PT887V', 'VV8VUVF9XA', 'W3W1ROPU95', 'W6XNV1LQB1', 'W7UAIA7T0S', 'W80YCK8SOZ', 'WCYEL3APQW', 'WIYGVFIKGG', 'WK78CI5LDO', 'WMJ41DPJ3Y', 'WPYHWCFES9', 'WR6V5879WH', 'WXJA24FUA7', 'WYMINM1LOF', 'X0W1ONC8B2', 'X1MS489ILB', 'X2AGMD6P59', 'X4HHLWE7SC', 'X4JDSBA0MD', 'X4JDSBA0MD', 'X5A3YUL1OI', 'X5WMJ9MUOV', 'X7566XEPLR', 'X8NG9V0HJ6', 'XASH2D2L7T', 'XG77M6C95Z', 'XLDFSI23XN', 'XTIR5SW5E8', 'Y0S85J0AZR', 'Y1DJA21YHZ', 'Y4HVM00UR9', 'YBIJBXII7J', 'YF6JV0C9MC', 'YP49T16HM1', 'YV85IDK7KQ', 'YY88X645NE', 'Z0075688AQ', 'Z0RE2PT1R9', 'Z1WO6PGN3E', 'Z9478QITNI', 'ZA6WQHDGBM', 'ZBG7MKFOZ1', 'ZC60EAOUMH', 'ZE7RW2IJMR', 'ZHRR08KKVM', 'ZIN73MJOAB', 'ZINL8265JQ', 'ZKG95SXJUU', 'ZKZ02PIDIU', 'ZM2Y4LL5QP', 'ZQRZQA3DZT', 'ZW84KJ7LJA'])]

    print(*pfs)
    #############################################
    #Fill in progression recurrence with info from current disease state if unknown
    pfs.loc[(pfs['PROG_RECUR_IND'].str.contains("Unknown/Not Applicable"))&(pfs['CURRENT_DISEASE_STATUS'].notnull()),'PROG_RECUR_IND'] = pfs['CURRENT_DISEASE_STATUS']
    #Grab current age if not listed
    pfs.loc[(pfs['AGE_AT_PROG_RECUR'].str.contains("Unknown/Not Applicable"))&(pfs['AGE_AT_CURRENT_DISEASE_STATUS'].notnull()),'AGE_AT_PROG_RECUR'] = pfs['AGE_AT_CURRENT_DISEASE_STATUS']
    #Replace PFS status strings to Progressed or disease free: Note has to be in the order with the occurence of "no" frequently. 
    pfs['AGE_AT_PROG_RECUR'] = pfs['AGE_AT_PROG_RECUR'].str.replace('Non-Evaluable', '')
    pfs['SOURCE'] = 'OUTCOMES-PROG_RECUR_IND-' + pfs['PROG_RECUR_IND']
    pfs['PROG_RECUR_IND'] = pfs['PROG_RECUR_IND'].str.replace('Non-Evaluable', '')
    pfs['PROG_RECUR_IND'] = pfs['PROG_RECUR_IND'].str.replace('Unknown/Not Applicable', '')
    #change active disease to account for filter below where active disease change within one month of med stop is progression. 
    pfs['PROG_RECUR_IND'] = pfs['PROG_RECUR_IND'].str.replace('Active Disease', '0:DiseaseFree')
    pfs['PROG_RECUR_IND'] = pfs['PROG_RECUR_IND'].str.replace('Progression', '1:Progressed')
    pfs['PROG_RECUR_IND'] = pfs['PROG_RECUR_IND'].str.replace('Recurrence', '1:Progressed')

    #Clean up no evidence of disease
    pfs['PROG_RECUR_IND'] = pfs['PROG_RECUR_IND'].str.replace('No', 'Replace')
    pfs['PROG_RECUR_IND'] = pfs['PROG_RECUR_IND'].str.replace('Replace Evidence of Disease', '0:DiseaseFree')
    pfs['PROG_RECUR_IND'] = pfs['PROG_RECUR_IND'].str.replace('Replace', '0:DiseaseFree')
    
    #Clenaup string ages
    pfs['AGE_AT_PROG_RECUR'] = pfs['AGE_AT_PROG_RECUR'].str.replace('Unknown/Not Applicable', '')
    pfs['AGE_AT_PROG_RECUR'] = pfs['AGE_AT_PROG_RECUR'].str.replace('Age 90 or older', '91')
    
    pfs['PROG_RECUR_IND'] = pfs['PROG_RECUR_IND'].str.replace('Non-Evaluable', '')
    pfs['PROG_RECUR_IND'] = pfs['PROG_RECUR_IND'].str.replace('Unknown/Not Applicable', '')
    #Renmane progression column
    pfs.rename(columns={'PROG_RECUR_IND': 'PFS_STATUS'}, inplace=True)
    
    
    #******************************************************
    #Drop duplicate progressions etc for respective start dates later on - duplicate below thru **************
    pfs.sort_values(by=['AGE_AT_PROG_RECUR'], ascending=False).groupby('PATIENT_ID')
    pfs_plat = pfs
    pfs_plat.to_csv('pfs_plat.txt', sep='\t', index=False)
    pfs_conjug = pfs
    pfs_conjug.to_csv('pfs_conjug.txt', sep='\t', index=False)
    
    #Merge in ICI start date to calculate PFS months with AGE_AT_PROG_RECUR' from above
    pfsMed_ICI = pd.read_csv('MedicationsPFStempICI.txt', sep='\t')
    pfsMed_ICI.rename(columns={'AGENT(S)': 'IMMUNOTHERAPY_AGENT'}, inplace=True)
    pfsMed_ICI.drop(['PFS_STATUS', 'PFS_MONTHS', 'MED_LINE_REGIMEN'], axis=1, inplace = True)
    #Merge ICI df with pfs dataframe. Medications_V4.csv
    pfs_ICI = pfs.merge(pfsMed_ICI, how='left', on='PATIENT_ID')
    pfs_ICI.to_csv('ICI_Lookheretemp.txt', sep='\t', index=False)
    #Convert objects to numeric
    pfs_ICI['AGE_AT_PROG_RECUR'] = pd.to_numeric(pfs_ICI.AGE_AT_PROG_RECUR, errors='coerce')
    pfs_ICI['AGE_AT_MED_START_Med'] = pd.to_numeric(pfs_ICI.AGE_AT_MED_START_Med, errors='coerce')
    #calculate pfs_months
    pfs_ICI['PFS_MONTHS'] = (pfs_ICI['AGE_AT_PROG_RECUR'] - pfs_ICI['AGE_AT_MED_START_Med'])*12
    #pfs.drop(['AGE_AT_CURRENT_DISEASE_STATUS', 'CHANGE_OF_TREATMENT', 'CURRENT_DISEASE_STATUS'], axis=1, inplace = True)
    pfs_ICI.drop(['CURRENT_DISEASE_STATUS'], axis=1, inplace = True)
    #note 'AGE_AT_PROG_RECUR' = pfs 'STOP_DATE'}, inplace=True
    pfs_ICI.rename(columns={'AGE_AT_PROG_RECUR': 'STOP_DATE'}, inplace=True)
    
    print('******************************************************** KEEP COLUMNS *******************************\n\n')
    print(*pfs)
    print('******************************************************** KEEP COLUMNS *******************************\n\n')
    
    pfs_ICI = pfs_ICI[['PATIENT_ID', 'IMMUNOTHERAPY_AGENT', 'PFS_STATUS', 'PFS_MONTHS', 'AGE_AT_MED_START_Med', 'STOP_DATE_Med', 'STOP_DATE', 'SOURCE']]
    pfs_ICI.rename(columns={'AGE_AT_MED_START_Med': 'AGE_AT_MED_START'}, inplace=True)
    print(*pfs_ICI)
    #Filter out patients without status
    pfs_ICI = pfs_ICI[pfs_ICI['PFS_STATUS'].isin(['0:DiseaseFree', '1:Progressed'])]
    #remove rows where progression etc happened before med start date.
    pfs_ICI = pfs_ICI[(pfs_ICI['PFS_MONTHS'] > 0)]
    pfs_ICI.sort_values(by=['STOP_DATE'], inplace = True, ascending=False)
    #Add calculation to select active disease progression within a month of med stop:
    pfs_ICI['PROGMEDSTOP'] = (pfs_ICI['STOP_DATE'] - pfs_ICI['STOP_DATE_Med'])*12
    pfs_ICI.loc[(pfs_ICI['SOURCE'].str.contains("OUTCOMES-PROG_RECUR_IND-Active Disease"))&(pfs_ICI['PROGMEDSTOP'].between(-1, 1)),'PFS_STATUS'] = '1:Progressed'
    pfs_ICI.drop(['PROGMEDSTOP'], axis=1, inplace = True)
    #pfs.groupby('PATIENT_ID', 'PFS_STATUS').apply(lambda x: x.sort_values('STOP_DATE'), inplace = True, ascending=False)
    print('\n\n************************************************************ Progression FILTER ********************************\n\n')
    #Select progressed lines for filtering
    df_PFS_outcomesProg = pfs_ICI.loc[pfs['PFS_STATUS'].isin(['1:Progressed'])]
    #df_PFS_outcomesProg.query('STOP_DATE' > 'AGE_AT_MED_STOP')
    #df_PFS_outcomesProg[df_PFS_outcomesProg.STOP_DATE.gt(df_PFS_outcomesProg.AGE_AT_MED_STOP)]
    df_PFS_outcomesProg = df_PFS_outcomesProg.sort_values(by=['PATIENT_ID', 'PFS_MONTHS'], ascending=True)
    ####keep earliest progression
    #df_PFS_outcomesProg.drop_duplicates(subset=['PATIENT_ID'], keep='first', inplace=True)
    df_PFS_outcomesProg.to_csv('PFScombinedProg.txt', sep='\t', index=False)
    print('\n\n************************************************************ Disease Free FILTER ********************************\n\n')

    df_PFS_outcomesDisFree = pfs_ICI.loc[pfs['PFS_STATUS'].isin(['0:DiseaseFree'])]
    #Eliminate outcomes where the stop date is 
    #df_PFS_outcomesDisFree[df_PFS_outcomesDisFree.STOP_DATE.gt(df_PFS_outcomesDisFree.AGE_AT_MED_STOP)]
    #Sort disease free with longest months on top 
    df_PFS_outcomesDisFree.sort_values(by=['PATIENT_ID', 'PFS_MONTHS'], ascending=False).groupby('PATIENT_ID')
    ####Keep longest disease free PFS
    #df_PFS_outcomesDisFree.drop_duplicates(subset=['PATIENT_ID'], keep='first', inplace=True)
    df_PFS_outcomesDisFree.to_csv('PFScombinedDFree.txt', sep='\t', index=False)

    Merge = [df_PFS_outcomesDisFree, df_PFS_outcomesProg]
    dfOutcomes = pd.concat(Merge)
    dfOutcomes.sort_values(by=['PATIENT_ID'], inplace=True)
    dfOutcomes = dfOutcomes.sort_values(by=['PATIENT_ID', 'PFS_MONTHS'], ascending=False)
    #keep the overal longest PFS shows response or recurrence
    dfOutcomes.drop_duplicates(subset=['PATIENT_ID'], keep='first', inplace=True)
    ######SAVE separate temp files for pfsMed_plat pfsMed_conjug pfsMed_ICI = pfsM
    dfOutcomes.rename(columns={'AGENT(S)': 'IMMUNOTHERAPY_AGENT'}, inplace=True)
    dfOutcomes.to_csv('PFS_OutcomesICI.txt', sep='\t', index=False)
#********************************* Outcomes with med start merge duplicate dataframes. 

tX = []
FilePatient7 = glob.glob('*_Imaging_V4.csv')
for f in FilePatient7:
    dfstartX = pd.read_csv(f, index_col=None, header=0)
    tX.append(dfstartX)
    pfsImag = pd.concat(tX, axis=0, ignore_index=True)    
#FileSample3 = glob.glob('*_Imaging_V4.csv')
#for f in FileSample3:
#    pfsImag = pd.read_csv(f)
    pfsImag = pfsImag[['AvatarKey', 'AgeAtImageScan', 'ImageLesionGrowth', 'ImagingModalityMethod']]
    print('\n\n************************************************************ Image for Progression ********************************\n\n')
    pfsImag.columns = [camel_to_snake(col) for col in pfsImag.columns]
    pfsImag.rename(columns = str.upper , inplace = True)
    #Rename imaging to account for datain source column
    pfsImag.rename(columns={'AVATAR_KEY': 'PATIENT_ID', 'IMAGING_MODALITY_METHOD': 'SOURCE'}, inplace=True)
    #Patient Filter
    #pfsImag = pfsImag[pfsImag['PATIENT_ID'].isin(['03GCDX0K5N', '06B7Y294XL', '07EF826VTZ', '0C6NW2TC07', '0G5X064EB6', '0GFLCECW75', '0MJYHKZ6AK', '0RA5URUMVK', '0RFVUNS0M6', '0XBPD7NGNT', '15B16LU7OM', '1H50O85Z25', '1HJFO6EFWS', '1IQTFPMLQR', '1K0JG7QIX8', '1K0JG7QIX8', '1Q58ISCSYD', '1RQGB4W7T4', '1YUCT7U163', '274XK9OLE5', '282INULRAK', '28MTFJ7TIV', '2C9ADE4KW6', '2HNKJUZDQR', '2JRJWC9O7M', '2L2ZRV1LZY', '2N04GOV5XK', '2R23NCXQE0', '2UN4T04H08', '37SHF80OE0', '3C2HGYMBKU', '3E15DF0JVG', '3TLQLPNZP8', '41PP6CPRV7', '45PTETBCPX', '45VCVC3MKN', '4AFTCY1HP0', '4B7HJMVOOV', '4FLIDTYACP', '4FLIDTYACP', '4FTJM67GLA', '4OGM5KMQX2', '4S1A75B0CZ', '4UOOLFGTJC', '4UY5RRTZDA', '506HNR0L9Y', '56QVAODIUI', '5CE11AIS9P', '5JWOSJEGEQ', '5UUQPS8Y5W', '5VD97K8TUD', '63JFZJ9HYR', '64F1NAL85E', '66LL4MEA0D', '69NEPW50R9', '6A9ORBLX71', '6DRRKFWIL1', '6JD9VGBWZB', '6PX4GRURXV', '6T5GP2FH7P', '6U163EB0I2', '70HLF9IQQQ', '70P236XSX2', '731EJZUG4W', '73BBJCCFCG', '75HWMQ9THU', '75QKEG7QT5', '761ATF11EM', '78SU3ENGQ3', '78SU3ENGQ3', '797KYKKNTB', '7DQNA9XN5I', '7JT8CIRL35', '7O8FXBKCMP', '7P2075M6NA', '7XZV89ZKWL', '80II7FE25J', '80TN2OH9U2', '8A9JJV56BE', '8C9KM403EG', '8I52ZJY44W', '8IPGJJGS93', '8K3F2U4AEG', '8SSG3I43TG', '8SWRWW6ZR8', '8U1UG57FYC', '8XOABM89XD', '8YEVOS4E0I', '8ZVLI42YTE', '913EQAJLBI', '923EC0J63U', '928IO9LH0Z', '960N6W495Z', '988C2SHC0G', '995JAM1PGX', '9BZRRRQBLI', '9DFUAVXQTU', '9FUJS5X3V0', '9H5WMVOS0K', '9K443LZVX9', '9RQJ8JFILJ', '9SKYB6E4B4', '9ZXCC6ULKK', '9ZXCC6ULKK', 'A049JYM87K', 'A3AM9T7VPE', 'ACLOR9FFS5', 'ALL5VBRLOM', 'AP93ZGZX2W', 'AVMUPL2GS0', 'B3D07F2E2S', 'BCZK5PTSG9', 'BK6DKQTSI5', 'BTHE1B9DZ2', 'BU8GCS3MYT', 'BV2N69FYI6', 'BV2N69FYI6', 'BXXKG0HZIK', 'C5A7KUNVCS', 'CADAN6NOJ4', 'CGDPFH8IBW', 'CK9QL0A5ET', 'CM03MT997O', 'CM03MT997O', 'CM03MT997O', 'CTGHR4UZ18', 'CW0ZVIVU0P', 'CZE113UZIR', 'D0WP0B97Q7', 'D2C5F4UMXY', 'D7KAPTXFOG', 'DCPKIAYXDB', 'DQF8CHSKPG', 'DYKRF4KH3H', 'E2H89HF5KZ', 'E4H964SOU9', 'EF7AOZI79O', 'EHOVA8D07G', 'EKSOBYXBVG', 'EKXLPR0080', 'EMVBZPMWST', 'EOZ9LELKBM', 'EQNGL8N9W7', 'EU4WQTA1KL', 'F1XUR9R1RK', 'F2TQDDM3HH', 'F905RYZFIO', 'FNVA6L545Q', 'FOTCBZ3OMJ', 'FPLKDOPOH3', 'FVZV8LMNI5', 'FXGR65TB4V', 'FXSR99KNNT', 'GA7GNCV7KM', 'GF6S6LZVF5', 'GH6KKNEPWX', 'GIJZDVWC2Q', 'GNJVFE3HRM', 'GNXLDBYWN1', 'GOYPZUZNAW', 'GPLOBS1HRX', 'GU9WE974XH', 'GWNOJYMTTB', 'HG7EO5096W', 'HUCN0XE91R', 'HWECMUCW0Z', 'HWIG0FSC6Z', 'HXKS3QWP6F', 'HY2ZUQI0SD', 'I2EM49R0AO', 'I65GB8ZEJ2', 'IMH5IF30JJ', 'IW3OCCZBPJ', 'IW7K4RUJG9', 'IXEPJTR5B9', 'J0D6I815XU', 'J15Q9C8MIE', 'J3C9QS0P78', 'J60EKEW4OF', 'J81JUVHO9L', 'J9NHJFS8IN', 'JBLZPPZSFH', 'JG2Z3K6UC5', 'JHBPEVZ48F', 'JHZSJNIKIP', 'JJGV4LPENY', 'JLGA08C3K5', 'JM3MBDWG5H', 'JNXYDUFEED', 'JPSK03CCDH', 'JPSK03CCDH', 'JQHVTFHS1X', 'JRC4RTZNTK', 'JSWZX87JYV', 'JZKKZ1OTQJ', 'K1DVC6O9C7', 'K60LO6G1FG', 'K8DFO8FTT1', 'K8OFSLEHIX', 'KAN1H1VWQG', 'KC8HPZPX7T', 'KPCHILTX70', 'KSTUVZCJIF', 'KW650F0557', 'L410MYBLDN', 'L6RMSU0BXG', 'L8HCOD7F64', 'LATEJ509E1', 'LLXVLRWEAE', 'LM4ZZPQQFO', 'LS6UIWKFEJ', 'LX5K9JZKLH', 'LXELBL78SE', 'M1D0DA5RWT', 'M9G3Y6A4X7', 'MG1SUQR5K8', 'MIVJ67B9F5', 'MJYLGG6W0K', 'MRZW12GMST', 'MYG3XG6SAQ', 'N4K1N7UUMM', 'N675BI25Y6', 'N8USV2HKJ6', 'NAMDA2GZXA', 'NAR8MDP01R', 'NDV4VM27O6', 'NE17FOVUC5', 'NK65MD354P', 'NOD926GZ41', 'NPUNSPR11V', 'O8RFN0K0CC', 'OA1K7SKWYW', 'OFMS9RN5W8', 'OFOSHGCQE4', 'OJQ5IBJCDB', 'ONHYTARPID', 'OPI88VIT01', 'OQYBJZXP4D', 'ORL6LF2IMD', 'OWNR0K9J5O', 'OXKYPBVO0Q', 'OXKYPBVO0Q', 'P0YOHVOGGY', 'P6F7GJOYOS', 'PCM0C3L6G1', 'PGVXTE18CS', 'PJ316U0AO6', 'PQDMEIJU8A', 'PTN0DBFVOH', 'Q1LFWVL6GK', 'Q1TM0Z8WWH', 'Q3XX8RO9CY', 'QB511DI6QC', 'QD6U65T2FH', 'QDTQI9DG2U', 'QJ80IGAAAA', 'QO27BH1WTV', 'QUEJJTOM01', 'QVEQP2FBPS', 'QVLKW6YJ9S', 'QZ68W571DZ', 'R0JJZPXHE6', 'R9EFJQNQD7', 'RC5DO9WG9A', 'RDNMW7PS5V', 'RL11WPPBNV', 'RN9YUHZDTE', 'RRWTZZSW3F', 'S414XWRR6V', 'S5AEWMYN2T', 'S5R5RW4L05', 'S7GI3EL7N1', 'S9CLR0UMQD', 'S9GIV41A0F', 'SF953NJFYV', 'SQV8WBN78P', 'STZJJQ7LKL', 'SUYGY0MFVB', 'SW0L8B0LS4', 'T3BB8UT9KP', 'T653MJNWUY', 'T8IWENLDFE', 'TEEDQQPZOI', 'TFJ3UKACLL', 'TG4OSS40W2', 'TGCXP2RE55', 'TTMBZYSHU8', 'TTN1Q2AA35', 'TTN1Q2AA35', 'TTUKMQJUXT', 'TX7XF0XMBL', 'TXJ9VXE4QC', 'TYUGO1UOO3', 'U0JP7MO0CU', 'U1KHQ8NI5I', 'U3YN74R7B9', 'U6SEJCBUKI', 'U86PUE8URG', 'UB8PDQ9GI1', 'UD8WP9Y59H', 'UFKCGRWOE6', 'ULH5Q1UQZ9', 'UMZ44YQOYC', 'UMZ44YQOYC', 'UOHLJNY3WH', 'UOWQ0WV182', 'UYBI9N0WBE', 'V3V21CTLZR', 'VE1W0KD4IA', 'VF0HIW4X63', 'VI0WMPQNAD', 'VIMLBEI2ZQ', 'VKOQ16NG0M', 'VKOQ16NG0M', 'VMKJ7V6SRD', 'VMKJ7V6SRD', 'VQS7PT887V', 'VV8VUVF9XA', 'W3W1ROPU95', 'W6XNV1LQB1', 'W7UAIA7T0S', 'W80YCK8SOZ', 'WCYEL3APQW', 'WIYGVFIKGG', 'WK78CI5LDO', 'WMJ41DPJ3Y', 'WPYHWCFES9', 'WR6V5879WH', 'WXJA24FUA7', 'WYMINM1LOF', 'X0W1ONC8B2', 'X1MS489ILB', 'X2AGMD6P59', 'X4HHLWE7SC', 'X4JDSBA0MD', 'X4JDSBA0MD', 'X5A3YUL1OI', 'X5WMJ9MUOV', 'X7566XEPLR', 'X8NG9V0HJ6', 'XASH2D2L7T', 'XG77M6C95Z', 'XLDFSI23XN', 'XTIR5SW5E8', 'Y0S85J0AZR', 'Y1DJA21YHZ', 'Y4HVM00UR9', 'YBIJBXII7J', 'YF6JV0C9MC', 'YP49T16HM1', 'YV85IDK7KQ', 'YY88X645NE', 'Z0075688AQ', 'Z0RE2PT1R9', 'Z1WO6PGN3E', 'Z9478QITNI', 'ZA6WQHDGBM', 'ZBG7MKFOZ1', 'ZC60EAOUMH', 'ZE7RW2IJMR', 'ZHRR08KKVM', 'ZIN73MJOAB', 'ZINL8265JQ', 'ZKG95SXJUU', 'ZKZ02PIDIU', 'ZM2Y4LL5QP', 'ZQRZQA3DZT', 'ZW84KJ7LJA'])]

    pfsImag = pfsImag[pfsImag['IMAGE_LESION_GROWTH'].isin(['Yes'])]
    pfsImag['SOURCE'] = 'IMAGING-' +  pfsImag['SOURCE']  
    #Build/change PFS status
    pfsImag.rename(columns={'SOURCE': 'IMAGE_INDICATION', 'IMAGE_LESION_GROWTH': 'PFS_STATUS'}, inplace=True)
    pfsImag['PFS_STATUS'] = pfsImag['PFS_STATUS'].str.replace('Yes', '1:Progressed')
    #Drop duplicate progressions etc for respective start dates later on - duplicate below thru **************
    pfs_platImag = pfsImag
    pfs_platImag.to_csv('pfs_platImag.txt', sep='\t', index=False)
    pfs_conjugImag = pfsImag
    pfs_conjugImag.to_csv('pfs_conjugImag.txt', sep='\t', index=False)
    
#**************************************************   DUPlicate DF for other drugs    
    #***** ICI - specific start dates
    pfsMed3 = pd.read_csv('MedicationsPFStempICI.txt', sep='\t')
    pfsMed3.drop(['PFS_STATUS', 'PFS_MONTHS'], axis=1, inplace = True)
    pfsImagICI = pfsImag.merge(pfsMed3, how='left', on='PATIENT_ID')
    print(*pfsImagICI)
    
    #pfsImagICI.drop(['PFS_STATUS', 'PFS_MONTHS', 'STOP_DATE_Med'], axis=1, inplace = True)
    #Drop patients with no ICI 
    pfsImagICI = pfsImagICI[pfsImagICI['AGENT(S)'].isin(['Atezolizumab', 'Avelumab', 'Nivolumab', 'Pembrolizumab', 'Durvalumab'])]
    pfsImagICI['AGE_AT_IMAGE_SCAN'] = pfsImag['AGE_AT_IMAGE_SCAN'].str.replace('Age 90 or older', '91')
    pfsImagICI['AGE_AT_IMAGE_SCAN'] = pd.to_numeric(pfsImag.AGE_AT_IMAGE_SCAN, errors='coerce')
    pfsImagICI = pfsImagICI.sort_values(by=['PATIENT_ID', 'AGE_AT_IMAGE_SCAN'], ascending=True)
    #remove rows where image-based progression etc happened before med start date.
    pfsImagICI['TempFilter'] = (pfsImagICI['AGE_AT_IMAGE_SCAN'] - pfsImagICI['AGE_AT_MED_START_Med'])*12
    pfsImagICI['TempFilter'] = pd.to_numeric(pfsImagICI.TempFilter, errors='coerce')
    pfsImagICI['TempFilter'] = pfsImagICI['TempFilter'].fillna(0)
    pfsImagICI.to_csv('pfsImag.txt', sep='\t', index=False)
    #keep first PFS event and those from after med start
    pfsImagICI = pfsImagICI[(pfsImagICI['TempFilter'] > 0)]
    pfsImagICI = pfsImagICI.sort_values(by=['PATIENT_ID', 'TempFilter'], ascending=True)
    #pfsImagICI.drop_duplicates(subset=['PATIENT_ID'], keep='first', inplace=True)
    pfsImagICI.drop(['TempFilter'], axis=1, inplace = True)
    
    print('\n\n********************************************************\n\n')
    print(*pfsImagICI)
    print('\n\n********************************************************\n\n')
    
    #pfsImagICI.to_csv('PFS_IMAGsort.txt', sep='\t', index=False)
    pfsImagICI.drop_duplicates(subset=['PATIENT_ID', 'AGE_AT_MED_START_Med'], keep='first', inplace=True)
    #Calculate PFS
    pfsImagICI['AGE_AT_IMAGE_SCAN'] = pd.to_numeric(pfsImagICI.AGE_AT_IMAGE_SCAN, errors='coerce')
    pfsImagICI['AGE_AT_MED_START_Med'] = pd.to_numeric(pfsImagICI.AGE_AT_MED_START_Med, errors='coerce')
    pfsImagICI['PFS_MONTHS'] = (pfsImagICI['AGE_AT_IMAGE_SCAN'] - pfsImagICI['AGE_AT_MED_START_Med'])*12
    pfsImagICI.to_csv('PFS_IMAGsort.txt', sep='\t', index=False)

    #Cleanup
    pfsImagICI.rename(columns={'AGE_AT_IMAGE_SCAN': 'STOP_DATE'}, inplace=True)
    #pfsImagICI['STOP_DATE_OTHER'] =  pfsImagICI['STOP_DATE']
    pfsImagICI.rename(columns={'AGE_AT_MED_START': 'AGE_AT_MED_START_Imag'}, inplace=True)
    pfsImagICI = pfsImagICI[['PATIENT_ID', 'AGENT(S)', 'PFS_STATUS', 'PFS_MONTHS', 'AGE_AT_MED_START_Med', 'STOP_DATE', 'IMAGE_INDICATION']]
    #cleanup Imaging
    pfsImagICI.to_csv('PFS_IMAG_ICI.txt', sep='\t', index=False)
#****************************************    
tXI = []
FilePatient8 = glob.glob('*_MetastaticDisease_V4.csv')
for f in FilePatient8:
    dfstartXI = pd.read_csv(f, index_col=None, header=0)
    tXI.append(dfstartXI)
    pfsMet = pd.concat(tXI, axis=0, ignore_index=True)    
#FileSample4 = glob.glob('*_MetastaticDisease_V4.csv')
#for f in FileSample4:
#    pfsMet = pd.read_csv(f)
    pfsMet = pfsMet[['AvatarKey', 'MetastaticDiseaseInd', 'AgeAtMetastaticSite']]
    print('\n\n************************************************************ Met for Progression ********************************\n\n')
    pfsMet.columns = [camel_to_snake(col) for col in pfsMet.columns]
    pfsMet.rename(columns = str.upper , inplace = True)
    pfsMet.rename(columns={'AVATAR_KEY': 'PATIENT_ID'}, inplace=True)
    #Patient Filter
    #pfsMet = pfsMet[pfsMet['PATIENT_ID'].isin(['03GCDX0K5N', '06B7Y294XL', '07EF826VTZ', '0C6NW2TC07', '0G5X064EB6', '0GFLCECW75', '0MJYHKZ6AK', '0RA5URUMVK', '0RFVUNS0M6', '0XBPD7NGNT', '15B16LU7OM', '1H50O85Z25', '1HJFO6EFWS', '1IQTFPMLQR', '1K0JG7QIX8', '1K0JG7QIX8', '1Q58ISCSYD', '1RQGB4W7T4', '1YUCT7U163', '274XK9OLE5', '282INULRAK', '28MTFJ7TIV', '2C9ADE4KW6', '2HNKJUZDQR', '2JRJWC9O7M', '2L2ZRV1LZY', '2N04GOV5XK', '2R23NCXQE0', '2UN4T04H08', '37SHF80OE0', '3C2HGYMBKU', '3E15DF0JVG', '3TLQLPNZP8', '41PP6CPRV7', '45PTETBCPX', '45VCVC3MKN', '4AFTCY1HP0', '4B7HJMVOOV', '4FLIDTYACP', '4FLIDTYACP', '4FTJM67GLA', '4OGM5KMQX2', '4S1A75B0CZ', '4UOOLFGTJC', '4UY5RRTZDA', '506HNR0L9Y', '56QVAODIUI', '5CE11AIS9P', '5JWOSJEGEQ', '5UUQPS8Y5W', '5VD97K8TUD', '63JFZJ9HYR', '64F1NAL85E', '66LL4MEA0D', '69NEPW50R9', '6A9ORBLX71', '6DRRKFWIL1', '6JD9VGBWZB', '6PX4GRURXV', '6T5GP2FH7P', '6U163EB0I2', '70HLF9IQQQ', '70P236XSX2', '731EJZUG4W', '73BBJCCFCG', '75HWMQ9THU', '75QKEG7QT5', '761ATF11EM', '78SU3ENGQ3', '78SU3ENGQ3', '797KYKKNTB', '7DQNA9XN5I', '7JT8CIRL35', '7O8FXBKCMP', '7P2075M6NA', '7XZV89ZKWL', '80II7FE25J', '80TN2OH9U2', '8A9JJV56BE', '8C9KM403EG', '8I52ZJY44W', '8IPGJJGS93', '8K3F2U4AEG', '8SSG3I43TG', '8SWRWW6ZR8', '8U1UG57FYC', '8XOABM89XD', '8YEVOS4E0I', '8ZVLI42YTE', '913EQAJLBI', '923EC0J63U', '928IO9LH0Z', '960N6W495Z', '988C2SHC0G', '995JAM1PGX', '9BZRRRQBLI', '9DFUAVXQTU', '9FUJS5X3V0', '9H5WMVOS0K', '9K443LZVX9', '9RQJ8JFILJ', '9SKYB6E4B4', '9ZXCC6ULKK', '9ZXCC6ULKK', 'A049JYM87K', 'A3AM9T7VPE', 'ACLOR9FFS5', 'ALL5VBRLOM', 'AP93ZGZX2W', 'AVMUPL2GS0', 'B3D07F2E2S', 'BCZK5PTSG9', 'BK6DKQTSI5', 'BTHE1B9DZ2', 'BU8GCS3MYT', 'BV2N69FYI6', 'BV2N69FYI6', 'BXXKG0HZIK', 'C5A7KUNVCS', 'CADAN6NOJ4', 'CGDPFH8IBW', 'CK9QL0A5ET', 'CM03MT997O', 'CM03MT997O', 'CM03MT997O', 'CTGHR4UZ18', 'CW0ZVIVU0P', 'CZE113UZIR', 'D0WP0B97Q7', 'D2C5F4UMXY', 'D7KAPTXFOG', 'DCPKIAYXDB', 'DQF8CHSKPG', 'DYKRF4KH3H', 'E2H89HF5KZ', 'E4H964SOU9', 'EF7AOZI79O', 'EHOVA8D07G', 'EKSOBYXBVG', 'EKXLPR0080', 'EMVBZPMWST', 'EOZ9LELKBM', 'EQNGL8N9W7', 'EU4WQTA1KL', 'F1XUR9R1RK', 'F2TQDDM3HH', 'F905RYZFIO', 'FNVA6L545Q', 'FOTCBZ3OMJ', 'FPLKDOPOH3', 'FVZV8LMNI5', 'FXGR65TB4V', 'FXSR99KNNT', 'GA7GNCV7KM', 'GF6S6LZVF5', 'GH6KKNEPWX', 'GIJZDVWC2Q', 'GNJVFE3HRM', 'GNXLDBYWN1', 'GOYPZUZNAW', 'GPLOBS1HRX', 'GU9WE974XH', 'GWNOJYMTTB', 'HG7EO5096W', 'HUCN0XE91R', 'HWECMUCW0Z', 'HWIG0FSC6Z', 'HXKS3QWP6F', 'HY2ZUQI0SD', 'I2EM49R0AO', 'I65GB8ZEJ2', 'IMH5IF30JJ', 'IW3OCCZBPJ', 'IW7K4RUJG9', 'IXEPJTR5B9', 'J0D6I815XU', 'J15Q9C8MIE', 'J3C9QS0P78', 'J60EKEW4OF', 'J81JUVHO9L', 'J9NHJFS8IN', 'JBLZPPZSFH', 'JG2Z3K6UC5', 'JHBPEVZ48F', 'JHZSJNIKIP', 'JJGV4LPENY', 'JLGA08C3K5', 'JM3MBDWG5H', 'JNXYDUFEED', 'JPSK03CCDH', 'JPSK03CCDH', 'JQHVTFHS1X', 'JRC4RTZNTK', 'JSWZX87JYV', 'JZKKZ1OTQJ', 'K1DVC6O9C7', 'K60LO6G1FG', 'K8DFO8FTT1', 'K8OFSLEHIX', 'KAN1H1VWQG', 'KC8HPZPX7T', 'KPCHILTX70', 'KSTUVZCJIF', 'KW650F0557', 'L410MYBLDN', 'L6RMSU0BXG', 'L8HCOD7F64', 'LATEJ509E1', 'LLXVLRWEAE', 'LM4ZZPQQFO', 'LS6UIWKFEJ', 'LX5K9JZKLH', 'LXELBL78SE', 'M1D0DA5RWT', 'M9G3Y6A4X7', 'MG1SUQR5K8', 'MIVJ67B9F5', 'MJYLGG6W0K', 'MRZW12GMST', 'MYG3XG6SAQ', 'N4K1N7UUMM', 'N675BI25Y6', 'N8USV2HKJ6', 'NAMDA2GZXA', 'NAR8MDP01R', 'NDV4VM27O6', 'NE17FOVUC5', 'NK65MD354P', 'NOD926GZ41', 'NPUNSPR11V', 'O8RFN0K0CC', 'OA1K7SKWYW', 'OFMS9RN5W8', 'OFOSHGCQE4', 'OJQ5IBJCDB', 'ONHYTARPID', 'OPI88VIT01', 'OQYBJZXP4D', 'ORL6LF2IMD', 'OWNR0K9J5O', 'OXKYPBVO0Q', 'OXKYPBVO0Q', 'P0YOHVOGGY', 'P6F7GJOYOS', 'PCM0C3L6G1', 'PGVXTE18CS', 'PJ316U0AO6', 'PQDMEIJU8A', 'PTN0DBFVOH', 'Q1LFWVL6GK', 'Q1TM0Z8WWH', 'Q3XX8RO9CY', 'QB511DI6QC', 'QD6U65T2FH', 'QDTQI9DG2U', 'QJ80IGAAAA', 'QO27BH1WTV', 'QUEJJTOM01', 'QVEQP2FBPS', 'QVLKW6YJ9S', 'QZ68W571DZ', 'R0JJZPXHE6', 'R9EFJQNQD7', 'RC5DO9WG9A', 'RDNMW7PS5V', 'RL11WPPBNV', 'RN9YUHZDTE', 'RRWTZZSW3F', 'S414XWRR6V', 'S5AEWMYN2T', 'S5R5RW4L05', 'S7GI3EL7N1', 'S9CLR0UMQD', 'S9GIV41A0F', 'SF953NJFYV', 'SQV8WBN78P', 'STZJJQ7LKL', 'SUYGY0MFVB', 'SW0L8B0LS4', 'T3BB8UT9KP', 'T653MJNWUY', 'T8IWENLDFE', 'TEEDQQPZOI', 'TFJ3UKACLL', 'TG4OSS40W2', 'TGCXP2RE55', 'TTMBZYSHU8', 'TTN1Q2AA35', 'TTN1Q2AA35', 'TTUKMQJUXT', 'TX7XF0XMBL', 'TXJ9VXE4QC', 'TYUGO1UOO3', 'U0JP7MO0CU', 'U1KHQ8NI5I', 'U3YN74R7B9', 'U6SEJCBUKI', 'U86PUE8URG', 'UB8PDQ9GI1', 'UD8WP9Y59H', 'UFKCGRWOE6', 'ULH5Q1UQZ9', 'UMZ44YQOYC', 'UMZ44YQOYC', 'UOHLJNY3WH', 'UOWQ0WV182', 'UYBI9N0WBE', 'V3V21CTLZR', 'VE1W0KD4IA', 'VF0HIW4X63', 'VI0WMPQNAD', 'VIMLBEI2ZQ', 'VKOQ16NG0M', 'VKOQ16NG0M', 'VMKJ7V6SRD', 'VMKJ7V6SRD', 'VQS7PT887V', 'VV8VUVF9XA', 'W3W1ROPU95', 'W6XNV1LQB1', 'W7UAIA7T0S', 'W80YCK8SOZ', 'WCYEL3APQW', 'WIYGVFIKGG', 'WK78CI5LDO', 'WMJ41DPJ3Y', 'WPYHWCFES9', 'WR6V5879WH', 'WXJA24FUA7', 'WYMINM1LOF', 'X0W1ONC8B2', 'X1MS489ILB', 'X2AGMD6P59', 'X4HHLWE7SC', 'X4JDSBA0MD', 'X4JDSBA0MD', 'X5A3YUL1OI', 'X5WMJ9MUOV', 'X7566XEPLR', 'X8NG9V0HJ6', 'XASH2D2L7T', 'XG77M6C95Z', 'XLDFSI23XN', 'XTIR5SW5E8', 'Y0S85J0AZR', 'Y1DJA21YHZ', 'Y4HVM00UR9', 'YBIJBXII7J', 'YF6JV0C9MC', 'YP49T16HM1', 'YV85IDK7KQ', 'YY88X645NE', 'Z0075688AQ', 'Z0RE2PT1R9', 'Z1WO6PGN3E', 'Z9478QITNI', 'ZA6WQHDGBM', 'ZBG7MKFOZ1', 'ZC60EAOUMH', 'ZE7RW2IJMR', 'ZHRR08KKVM', 'ZIN73MJOAB', 'ZINL8265JQ', 'ZKG95SXJUU', 'ZKZ02PIDIU', 'ZM2Y4LL5QP', 'ZQRZQA3DZT', 'ZW84KJ7LJA'])]

    print(*pfsMet)
    #Take all available mets by filtering those without Mets - used double negative as comments were variable for "yes/positive"
    pfsMet = pfsMet[~pfsMet['METASTATIC_DISEASE_IND'].isin(['No'])]
    #drop dulicate mets
    pfsMet.drop_duplicates(subset=['PATIENT_ID', 'AGE_AT_METASTATIC_SITE'], keep='first', inplace=True)
    pfsMet.sort_values(by=['AGE_AT_METASTATIC_SITE'], ascending=False).groupby('PATIENT_ID')
    #Drop samples with no dates
    pfsMet = pfsMet[~pfsMet['AGE_AT_METASTATIC_SITE'].isin(['Unknown/Not Applicable'])]
    #BUILD duplicated DFs for other drug calsses. 
    pfsMet_plat = pfsMet
    pfsMet_plat.to_csv('pfsMet_plat.txt', sep='\t', index=False)
    pfsMet_conjug = pfsMet
    pfsMet_conjug.to_csv('pfsMet_conjug.txt', sep='\t', index=False)

#MET DUPLICATE dataframe ************************************               
    #Merge in medication/ICI start dates
    pfsMed4 = pd.read_csv('MedicationsPFStempICI.txt', sep='\t')
    pfsMetICI = pfsMet.merge(pfsMed3, how='left', on='PATIENT_ID')
    pfsMetICI.drop_duplicates(subset=['PATIENT_ID', 'AGE_AT_METASTATIC_SITE'], inplace=True)
    pfsMetICI.to_csv('pfsMetICITemperr.txt', sep='\t', index=False)
    #Calculate MFS
    pfsMetICI['AGE_AT_METASTATIC_SITE'] = pd.to_numeric(pfsMetICI.AGE_AT_METASTATIC_SITE, errors='coerce')
    pfsMetICI['AGE_AT_MED_START_Med'] = pd.to_numeric(pfsMetICI.AGE_AT_MED_START_Med, errors='coerce')
    pfsMetICI['PFS_MONTHS'] = (pfsMetICI['AGE_AT_METASTATIC_SITE'] - pfsMetICI['AGE_AT_MED_START_Med'])*12
    #clear all met events before med start
    pfsMetICI = pfsMetICI[(pfsMetICI['PFS_MONTHS'] > 0)]
    #Drop patients with no ICI 
    pfsMetICI = pfsMetICI[pfsMetICI['AGENT(S)'].isin(['Atezolizumab', 'Avelumab', 'Nivolumab', 'Pembrolizumab', 'Durvalumab'])]
    #Cleanup
    pfsMetICI['PFS_STATUS'] = '1:Progressed'
    pfsMetICI.rename(columns={'METASTATIC_DISEASE_IND': 'SOURCE', 'AGE_AT_METASTATIC_SITE': 'STOP_DATE'}, inplace=True)
    pfsMetICI['STOP_DATE_OTHER'] =  pfsMetICI['STOP_DATE']
    pfsMetICI['SOURCE'] = pfsMetICI['SOURCE'].str.replace('Yes', 'METASTATIC_DISEASE-METASTATIC_DISEASE_IND')
    pfsMetICI.to_csv('pfsMetICI.txt', sep='\t', index=False)
    pfsMetICI = pfsMetICI[['PATIENT_ID', 'AGENT(S)', 'PFS_STATUS', 'PFS_MONTHS', 'AGE_AT_MED_START_Med', 'STOP_DATE_Med', 'STOP_DATE', 'SOURCE', 'STOP_DATE_OTHER']]
    pfsMetICI.rename(columns={'AGE_AT_MED_START': 'AGE_AT_MED_START_Met', 'AGE_AT_MED_STOP': 'AGE_AT_MED_STOP_Met'}, inplace=True)
    pfsMetICI = pfsMetICI[(pfsMetICI['PFS_MONTHS'] > 0)]
    pfsMetICI.rename(columns={'SOURCE': 'MET_INDICATION'}, inplace=True)
    pfsMetICI.to_csv('pfsMetICI.txt', sep='\t', index=False)
#MET DUPLICATE dataframe ************************************

dfOutcomes.rename(columns={'SOURCE': 'OUTCOMES_INDICATION'}, inplace=True)
#cleanup meds
pfsMed_ICI = pd.read_csv('MedicationsPFStempICI.txt', sep='\t')
pfsMed_ICI['PFS_MONTHS'] = pfsMed_ICI[pfsMed_ICI['PFS_MONTHS'] > 0]['PFS_MONTHS']
pfsMed_ICI['PFS_MONTHS'].replace('', np.nan, inplace=True)
pfsMed_ICI = pfsMed_ICI[pfsMed_ICI['PFS_MONTHS'].notna()]
pfsMed_ICI.rename(columns={'SOURCE': 'MEDICATION_INDICATION'}, inplace=True)

#Merge pfs, pfsMed, pfsImag, pfsMet
df_PFSTransposed = dfOutcomes.merge(pfsMed_ICI, how='outer', on='PATIENT_ID').merge(pfsImagICI, how='outer', on='PATIENT_ID').merge(pfsMetICI, how='outer', on='PATIENT_ID')
df_PFSTransposed.to_csv('PFScombinedTemp.txt', sep='\t', index=False)
#Sort dataframe to filter simple rows without multiple PFS indications
df_PFSTransposed.sort_values(by=['MEDICATION_INDICATION', 'IMAGE_INDICATION', 'MET_INDICATION'], inplace=True)
df_PFSTransposed.to_csv('PFScombined.txt', sep='\t', index=False)

####MAKE PROGRESSION FILE and merge all progression data
print('\n\n************************************************************  Progression ICI Merge ********************************\n\n')
pfsMed_ICI.rename(columns={'AGENT(S)': 'IMMUNOTHERAPY_AGENT', 'AGE_AT_MED_START_Med': 'AGE_AT_MED_START', 'STOP_DATE_Med': 'STOP_DATE', 'MEDICATION_INDICATION': 'INDICATION'}, inplace=True)
print(*pfsMed_ICI)
pfsMed_ICI.to_csv('pfsMed_ICI.FInal.txt', sep='\t', index=False)
pfsMetICI.rename(columns={'AGENT(S)': 'IMMUNOTHERAPY_AGENT', 'AGE_AT_MED_START_Med': 'AGE_AT_MED_START', 'MET_INDICATION': 'INDICATION'}, inplace=True)
#pfsMet.drop(['AGE_AT_METASTATIC_SITE'], axis=1, inplace = True)
print(*pfsMetICI)
pfsMetICI.to_csv('pfsMetICI.FInal.txt', sep='\t', index=False)
pfsImagICI.rename(columns={'AGENT(S)': 'IMMUNOTHERAPY_AGENT', 'AGE_AT_MED_START_Med': 'AGE_AT_MED_START', 'STOP_DATE_Med': 'STOP_DATE', 'IMAGE_INDICATION': 'INDICATION'}, inplace=True)
#pfsImag.drop(['AGE_AT_IMAGE_SCAN'], axis=1, inplace = True)
print(*pfsImagICI)
pfsImagICI.to_csv('pfsImagICI.FInal.txt', sep='\t', index=False)
df_PFS_outcomesProg.drop(['STOP_DATE_Med'], axis=1, inplace = True)
df_PFS_outcomesProg.rename(columns={'SOURCE': 'INDICATION'}, inplace=True)
print(*df_PFS_outcomesProg)
#concatenate and take earliest. 
Progression = [df_PFS_outcomesProg, pfsMed_ICI, pfsImagICI, pfsMetICI]
dfProgFinal = pd.concat(Progression)
dfProgFinal = dfProgFinal.sort_values(by=['PATIENT_ID', 'PFS_MONTHS'], ascending=False)
#dfProgFinal['IMMUNOTHERAPY_AGENT'] = dfProgFinal['IMMUNOTHERAPY_AGENT'].replace('', pd.NA).fillna(dfProgFinal['AGENT(S)'])
dfProgFinal.to_csv('dfProgFinalICI.txt', sep='\t', index=False)
####Keep longest disease free PFS
#dfProgFinal.drop_duplicates(subset=['PATIENT_ID'], keep='first', inplace=True)
print(*dfProgFinal)
#Bring in disease free
df_PFS_outcomesDisFree.drop(['STOP_DATE_Med'], axis=1, inplace = True)
df_PFS_outcomesDisFree.rename(columns={'SOURCE': 'INDICATION'}, inplace=True)
print(*df_PFS_outcomesDisFree)
Outcomes = [df_PFS_outcomesDisFree, dfProgFinal]
pfsFinal = pd.concat(Outcomes)
#Sort disease free with longest months on top 
pfsFinal = pfsFinal.sort_values(by=['PATIENT_ID', 'PFS_MONTHS'], ascending=False)
#pfsFinal = pfsFinal.sort_values(by=['PFS_MONTHS'], ascending=False)
####Keep longest disease free PFS
#pfsFinal.drop_duplicates(subset=['PATIENT_ID'], keep='first', inplace=True)

#remove patients with age estimations
pfsFinal = pfsFinal[~pfsFinal['STOP_DATE'].isin([91])]
pfsFinal.rename(columns={'INDICATION': 'PFS_INDICATION_ICI'}, inplace=True)


pfsFinal.to_csv('PFS_Final.txt', sep='\t')
#Split PFS into separate docs for filtering logic grabbing unique patientIDs for temp file
for patient in pfsFinal['PATIENT_ID'].unique():
    sub_dff = pfsFinal[pfsFinal['PATIENT_ID'] == patient]
    if sub_dff['PFS_STATUS'].isin(['1:Progressed']).all():
    #if sub_dff['PFS_STATUS'].eq(sub_dff['PFS_STATUS'].iloc[0]).all():
    #if len(np.unique(sub_dff.PFS_STATUS))=='1:Progressed':
        sub_dff = sub_dff.drop_duplicates(subset=['PATIENT_ID'], keep='last')
    if sub_dff['PFS_STATUS'].isin(['0:DiseaseFree']).all():
        sub_dff = sub_dff.drop_duplicates(subset=['PATIENT_ID'], keep='first')
    #if sub_dff[sub_dff.PFS_STATUS == (['1:Progressed'])].iloc[0]:
        #sub_dff = sub_dff.drop_duplicates(subset=['PATIENT_ID'], keep='first')
    sub_dff.loc[:(sub_dff == '0:DiseaseFree').any(1).idxmax()]
sub_dff.to_csv(f'{patient}_ICI_data.txt', sep='\t', index=False)
    
filenames = glob.glob('*ICI_data.txt')

dff2 = (pd.read_csv(f, index_col=None, header=0, nrows=1, sep='\t') for f in filenames)
concatenated_dfICI = pd.concat(dff2, ignore_index=True)
concatenated_dfICI.to_csv('PFS_FinalTrimICI.txt', sep='\t', index=False)
concatenated_dfICI.rename(columns={'PFS_STATUS': 'PFS_STATUS_ICI', 'PFS_MONTHS': 'PFS_MONTHS_ICI', 'AGE_AT_MED_START': 'AGE_AT_MED_START_ICI', 'STOP_DATE_Med': 'STOP_DATE_Med_ICI', 'MED_LINE_REGIMEN': 'MED_LINE_REGIMEN_ICI'}, inplace=True)
concatenated_dfICI.drop(['STOP_DATE_OTHER'], axis=1, inplace = True)


print('\n\n**********************************   MERGE ALL PFS data      ****************************\n\n\n\n')

#Merge to PatientCDE
PatientCDE = PatientCDE.merge(concatenated_dfICI, how='left', on='PATIENT_ID')
PatientCDE.drop(['AGE_AT_DIAGNOSIS'], axis=1, inplace = True)






#
#
#
###################PLATINUM
#keep first Platinum medication
pfsMed_plat.sort_values(by=['AGE_AT_MED_START'], inplace = True)
pfsMed_plat.drop_duplicates(subset=['PATIENT_ID'], keep='first',  inplace=True)
#Rename columns for PFS status
pfsMed_plat.rename(columns={'CHANGE_OF_TREATMENT': 'PFS_STATUS'}, inplace=True)
pfsMed_plat['AGE_AT_MED_START'] = pd.to_numeric(pfsMed_plat.AGE_AT_MED_START, errors='coerce')
pfsMed_plat['AGE_AT_MED_STOP'] = pd.to_numeric(pfsMed_plat.AGE_AT_MED_STOP, errors='coerce')
pfsMed_plat['PFS_MONTHS'] = (pfsMed_plat['AGE_AT_MED_STOP'] - pfsMed_plat['AGE_AT_MED_START'])*12
pfsMed_plat['MEDICATION_INDICATION'] = 'Medication(s) Change'
pfsMed_plat.rename(columns={'AGE_AT_MED_START': 'AGE_AT_MED_START_Med', 'AGE_AT_MED_STOP': 'STOP_DATE_Med', 'MEDICATION': 'AGENT(S)'}, inplace=True)
    
#Save medical file for other PFS intersection with therapy    
pfsMed_plat.to_csv('MedicationsPFStempPLAT.txt', sep='\t', index=False)
#Build Progression file from change in treatment
pfsMed_plat = pfsMed_plat[pfsMed_plat['PFS_STATUS'].isin(['Yes, due to progression'])]
pfsMed_plat['PFS_STATUS'] = pfsMed_plat['PFS_STATUS'].str.replace('Yes, due to progression', '1:Progressed')
pfsMed_plat.drop(['MED_LINE_REGIMEN'], axis=1, inplace = True)
#add source column
pfsMed_plat['SOURCE'] = 'MEDICATIONS-CHANGE_OF_TREATMENT'
pfsMed_plat = pfsMed_plat[['PATIENT_ID', 'AGENT(S)', 'PFS_STATUS', 'PFS_MONTHS', 'AGE_AT_MED_START_Med', 'STOP_DATE_Med', 'SOURCE']]
pfsMed_plat.to_csv('PFSmed.txt', sep='\t', index=False)
    

print('\n\n\n######################################################## STARTING PFS OUTCOMES EXTRACTION  ###################################################\n\n')

t8 = []
FilePatient5 = glob.glob('*_Outcomes_V4.csv')
for f in FilePatient5:
    dfstart8 = pd.read_csv(f, index_col=None, header=0)
    t8.append(dfstart8)
    pfs = pd.concat(t8, axis=0, ignore_index=True)
#FileSample2 = glob.glob('*_Outcomes_V4.csv')
#for f in FileSample2:
#    pfs = pd.read_csv(f)
    #cleanup unused columns
    pfs = pfs[['AvatarKey', 'ProgRecurInd', 'AgeAtProgRecur', 'AgeAtCurrentDiseaseStatus', 'CurrentDiseaseStatus']]
    #Convert headers to snake_case and upper case
    pfs.columns = [camel_to_snake(col) for col in pfs.columns]
    pfs.rename(columns = str.upper , inplace = True)
    pfs.rename(columns={'AVATAR_KEY': 'PATIENT_ID'}, inplace=True)
    #Patient Filter
    #pfs = pfs[pfs['PATIENT_ID'].isin(['03GCDX0K5N', '06B7Y294XL', '07EF826VTZ', '0C6NW2TC07', '0G5X064EB6', '0GFLCECW75', '0MJYHKZ6AK', '0RA5URUMVK', '0RFVUNS0M6', '0XBPD7NGNT', '15B16LU7OM', '1H50O85Z25', '1HJFO6EFWS', '1IQTFPMLQR', '1K0JG7QIX8', '1K0JG7QIX8', '1Q58ISCSYD', '1RQGB4W7T4', '1YUCT7U163', '274XK9OLE5', '282INULRAK', '28MTFJ7TIV', '2C9ADE4KW6', '2HNKJUZDQR', '2JRJWC9O7M', '2L2ZRV1LZY', '2N04GOV5XK', '2R23NCXQE0', '2UN4T04H08', '37SHF80OE0', '3C2HGYMBKU', '3E15DF0JVG', '3TLQLPNZP8', '41PP6CPRV7', '45PTETBCPX', '45VCVC3MKN', '4AFTCY1HP0', '4B7HJMVOOV', '4FLIDTYACP', '4FLIDTYACP', '4FTJM67GLA', '4OGM5KMQX2', '4S1A75B0CZ', '4UOOLFGTJC', '4UY5RRTZDA', '506HNR0L9Y', '56QVAODIUI', '5CE11AIS9P', '5JWOSJEGEQ', '5UUQPS8Y5W', '5VD97K8TUD', '63JFZJ9HYR', '64F1NAL85E', '66LL4MEA0D', '69NEPW50R9', '6A9ORBLX71', '6DRRKFWIL1', '6JD9VGBWZB', '6PX4GRURXV', '6T5GP2FH7P', '6U163EB0I2', '70HLF9IQQQ', '70P236XSX2', '731EJZUG4W', '73BBJCCFCG', '75HWMQ9THU', '75QKEG7QT5', '761ATF11EM', '78SU3ENGQ3', '78SU3ENGQ3', '797KYKKNTB', '7DQNA9XN5I', '7JT8CIRL35', '7O8FXBKCMP', '7P2075M6NA', '7XZV89ZKWL', '80II7FE25J', '80TN2OH9U2', '8A9JJV56BE', '8C9KM403EG', '8I52ZJY44W', '8IPGJJGS93', '8K3F2U4AEG', '8SSG3I43TG', '8SWRWW6ZR8', '8U1UG57FYC', '8XOABM89XD', '8YEVOS4E0I', '8ZVLI42YTE', '913EQAJLBI', '923EC0J63U', '928IO9LH0Z', '960N6W495Z', '988C2SHC0G', '995JAM1PGX', '9BZRRRQBLI', '9DFUAVXQTU', '9FUJS5X3V0', '9H5WMVOS0K', '9K443LZVX9', '9RQJ8JFILJ', '9SKYB6E4B4', '9ZXCC6ULKK', '9ZXCC6ULKK', 'A049JYM87K', 'A3AM9T7VPE', 'ACLOR9FFS5', 'ALL5VBRLOM', 'AP93ZGZX2W', 'AVMUPL2GS0', 'B3D07F2E2S', 'BCZK5PTSG9', 'BK6DKQTSI5', 'BTHE1B9DZ2', 'BU8GCS3MYT', 'BV2N69FYI6', 'BV2N69FYI6', 'BXXKG0HZIK', 'C5A7KUNVCS', 'CADAN6NOJ4', 'CGDPFH8IBW', 'CK9QL0A5ET', 'CM03MT997O', 'CM03MT997O', 'CM03MT997O', 'CTGHR4UZ18', 'CW0ZVIVU0P', 'CZE113UZIR', 'D0WP0B97Q7', 'D2C5F4UMXY', 'D7KAPTXFOG', 'DCPKIAYXDB', 'DQF8CHSKPG', 'DYKRF4KH3H', 'E2H89HF5KZ', 'E4H964SOU9', 'EF7AOZI79O', 'EHOVA8D07G', 'EKSOBYXBVG', 'EKXLPR0080', 'EMVBZPMWST', 'EOZ9LELKBM', 'EQNGL8N9W7', 'EU4WQTA1KL', 'F1XUR9R1RK', 'F2TQDDM3HH', 'F905RYZFIO', 'FNVA6L545Q', 'FOTCBZ3OMJ', 'FPLKDOPOH3', 'FVZV8LMNI5', 'FXGR65TB4V', 'FXSR99KNNT', 'GA7GNCV7KM', 'GF6S6LZVF5', 'GH6KKNEPWX', 'GIJZDVWC2Q', 'GNJVFE3HRM', 'GNXLDBYWN1', 'GOYPZUZNAW', 'GPLOBS1HRX', 'GU9WE974XH', 'GWNOJYMTTB', 'HG7EO5096W', 'HUCN0XE91R', 'HWECMUCW0Z', 'HWIG0FSC6Z', 'HXKS3QWP6F', 'HY2ZUQI0SD', 'I2EM49R0AO', 'I65GB8ZEJ2', 'IMH5IF30JJ', 'IW3OCCZBPJ', 'IW7K4RUJG9', 'IXEPJTR5B9', 'J0D6I815XU', 'J15Q9C8MIE', 'J3C9QS0P78', 'J60EKEW4OF', 'J81JUVHO9L', 'J9NHJFS8IN', 'JBLZPPZSFH', 'JG2Z3K6UC5', 'JHBPEVZ48F', 'JHZSJNIKIP', 'JJGV4LPENY', 'JLGA08C3K5', 'JM3MBDWG5H', 'JNXYDUFEED', 'JPSK03CCDH', 'JPSK03CCDH', 'JQHVTFHS1X', 'JRC4RTZNTK', 'JSWZX87JYV', 'JZKKZ1OTQJ', 'K1DVC6O9C7', 'K60LO6G1FG', 'K8DFO8FTT1', 'K8OFSLEHIX', 'KAN1H1VWQG', 'KC8HPZPX7T', 'KPCHILTX70', 'KSTUVZCJIF', 'KW650F0557', 'L410MYBLDN', 'L6RMSU0BXG', 'L8HCOD7F64', 'LATEJ509E1', 'LLXVLRWEAE', 'LM4ZZPQQFO', 'LS6UIWKFEJ', 'LX5K9JZKLH', 'LXELBL78SE', 'M1D0DA5RWT', 'M9G3Y6A4X7', 'MG1SUQR5K8', 'MIVJ67B9F5', 'MJYLGG6W0K', 'MRZW12GMST', 'MYG3XG6SAQ', 'N4K1N7UUMM', 'N675BI25Y6', 'N8USV2HKJ6', 'NAMDA2GZXA', 'NAR8MDP01R', 'NDV4VM27O6', 'NE17FOVUC5', 'NK65MD354P', 'NOD926GZ41', 'NPUNSPR11V', 'O8RFN0K0CC', 'OA1K7SKWYW', 'OFMS9RN5W8', 'OFOSHGCQE4', 'OJQ5IBJCDB', 'ONHYTARPID', 'OPI88VIT01', 'OQYBJZXP4D', 'ORL6LF2IMD', 'OWNR0K9J5O', 'OXKYPBVO0Q', 'OXKYPBVO0Q', 'P0YOHVOGGY', 'P6F7GJOYOS', 'PCM0C3L6G1', 'PGVXTE18CS', 'PJ316U0AO6', 'PQDMEIJU8A', 'PTN0DBFVOH', 'Q1LFWVL6GK', 'Q1TM0Z8WWH', 'Q3XX8RO9CY', 'QB511DI6QC', 'QD6U65T2FH', 'QDTQI9DG2U', 'QJ80IGAAAA', 'QO27BH1WTV', 'QUEJJTOM01', 'QVEQP2FBPS', 'QVLKW6YJ9S', 'QZ68W571DZ', 'R0JJZPXHE6', 'R9EFJQNQD7', 'RC5DO9WG9A', 'RDNMW7PS5V', 'RL11WPPBNV', 'RN9YUHZDTE', 'RRWTZZSW3F', 'S414XWRR6V', 'S5AEWMYN2T', 'S5R5RW4L05', 'S7GI3EL7N1', 'S9CLR0UMQD', 'S9GIV41A0F', 'SF953NJFYV', 'SQV8WBN78P', 'STZJJQ7LKL', 'SUYGY0MFVB', 'SW0L8B0LS4', 'T3BB8UT9KP', 'T653MJNWUY', 'T8IWENLDFE', 'TEEDQQPZOI', 'TFJ3UKACLL', 'TG4OSS40W2', 'TGCXP2RE55', 'TTMBZYSHU8', 'TTN1Q2AA35', 'TTN1Q2AA35', 'TTUKMQJUXT', 'TX7XF0XMBL', 'TXJ9VXE4QC', 'TYUGO1UOO3', 'U0JP7MO0CU', 'U1KHQ8NI5I', 'U3YN74R7B9', 'U6SEJCBUKI', 'U86PUE8URG', 'UB8PDQ9GI1', 'UD8WP9Y59H', 'UFKCGRWOE6', 'ULH5Q1UQZ9', 'UMZ44YQOYC', 'UMZ44YQOYC', 'UOHLJNY3WH', 'UOWQ0WV182', 'UYBI9N0WBE', 'V3V21CTLZR', 'VE1W0KD4IA', 'VF0HIW4X63', 'VI0WMPQNAD', 'VIMLBEI2ZQ', 'VKOQ16NG0M', 'VKOQ16NG0M', 'VMKJ7V6SRD', 'VMKJ7V6SRD', 'VQS7PT887V', 'VV8VUVF9XA', 'W3W1ROPU95', 'W6XNV1LQB1', 'W7UAIA7T0S', 'W80YCK8SOZ', 'WCYEL3APQW', 'WIYGVFIKGG', 'WK78CI5LDO', 'WMJ41DPJ3Y', 'WPYHWCFES9', 'WR6V5879WH', 'WXJA24FUA7', 'WYMINM1LOF', 'X0W1ONC8B2', 'X1MS489ILB', 'X2AGMD6P59', 'X4HHLWE7SC', 'X4JDSBA0MD', 'X4JDSBA0MD', 'X5A3YUL1OI', 'X5WMJ9MUOV', 'X7566XEPLR', 'X8NG9V0HJ6', 'XASH2D2L7T', 'XG77M6C95Z', 'XLDFSI23XN', 'XTIR5SW5E8', 'Y0S85J0AZR', 'Y1DJA21YHZ', 'Y4HVM00UR9', 'YBIJBXII7J', 'YF6JV0C9MC', 'YP49T16HM1', 'YV85IDK7KQ', 'YY88X645NE', 'Z0075688AQ', 'Z0RE2PT1R9', 'Z1WO6PGN3E', 'Z9478QITNI', 'ZA6WQHDGBM', 'ZBG7MKFOZ1', 'ZC60EAOUMH', 'ZE7RW2IJMR', 'ZHRR08KKVM', 'ZIN73MJOAB', 'ZINL8265JQ', 'ZKG95SXJUU', 'ZKZ02PIDIU', 'ZM2Y4LL5QP', 'ZQRZQA3DZT', 'ZW84KJ7LJA'])]

    print(*pfs)
    #############################################
    #Fill in progression recurrence with info from current disease state if unknown
    pfs.loc[(pfs['PROG_RECUR_IND'].str.contains("Unknown/Not Applicable"))&(pfs['CURRENT_DISEASE_STATUS'].notnull()),'PROG_RECUR_IND'] = pfs['CURRENT_DISEASE_STATUS']
    #Grab current age if not listed
    pfs.loc[(pfs['AGE_AT_PROG_RECUR'].str.contains("Unknown/Not Applicable"))&(pfs['AGE_AT_CURRENT_DISEASE_STATUS'].notnull()),'AGE_AT_PROG_RECUR'] = pfs['AGE_AT_CURRENT_DISEASE_STATUS']
    #Replace PFS status strings to Progressed or disease free: Note has to be in the order with the occurence of "no" frequently. 
    pfs['AGE_AT_PROG_RECUR'] = pfs['AGE_AT_PROG_RECUR'].str.replace('Non-Evaluable', '')
    pfs['SOURCE'] = 'OUTCOMES-PROG_RECUR_IND-' + pfs['PROG_RECUR_IND']
    pfs['PROG_RECUR_IND'] = pfs['PROG_RECUR_IND'].str.replace('Non-Evaluable', '')
    pfs['PROG_RECUR_IND'] = pfs['PROG_RECUR_IND'].str.replace('Unknown/Not Applicable', '')
    #change active disease to account for filter below where active disease change within one month of med stop is progression. 
    pfs['PROG_RECUR_IND'] = pfs['PROG_RECUR_IND'].str.replace('Active Disease', '0:DiseaseFree')
    pfs['PROG_RECUR_IND'] = pfs['PROG_RECUR_IND'].str.replace('Progression', '1:Progressed')
    pfs['PROG_RECUR_IND'] = pfs['PROG_RECUR_IND'].str.replace('Recurrence', '1:Progressed')

    #Clean up no evidence of disease
    pfs['PROG_RECUR_IND'] = pfs['PROG_RECUR_IND'].str.replace('No', 'Replace')
    pfs['PROG_RECUR_IND'] = pfs['PROG_RECUR_IND'].str.replace('Replace Evidence of Disease', '0:DiseaseFree')
    pfs['PROG_RECUR_IND'] = pfs['PROG_RECUR_IND'].str.replace('Replace', '0:DiseaseFree')
    
    #Clenaup string ages
    pfs['AGE_AT_PROG_RECUR'] = pfs['AGE_AT_PROG_RECUR'].str.replace('Unknown/Not Applicable', '')
    pfs['AGE_AT_PROG_RECUR'] = pfs['AGE_AT_PROG_RECUR'].str.replace('Age 90 or older', '91')
    
    pfs['PROG_RECUR_IND'] = pfs['PROG_RECUR_IND'].str.replace('Non-Evaluable', '')
    pfs['PROG_RECUR_IND'] = pfs['PROG_RECUR_IND'].str.replace('Unknown/Not Applicable', '')
    #Renmane progression column
    pfs.rename(columns={'PROG_RECUR_IND': 'PFS_STATUS'}, inplace=True)
    
    
    #******************************************************
    #Drop duplicate progressions etc for respective start dates later on - duplicate below thru **************
    pfs.sort_values(by=['AGE_AT_PROG_RECUR'], ascending=False).groupby('PATIENT_ID')
    
    #Merge in PLAT start date to calculate PFS months with AGE_AT_PROG_RECUR' from above
    pfsMed_plat = pd.read_csv('MedicationsPFStempPLAT.txt', sep='\t')
    pfsMed_plat.rename(columns={'AGENT(S)': 'PLATINUM_AGENT'}, inplace=True)
    pfsMed_plat.drop(['PFS_STATUS', 'PFS_MONTHS', 'MED_LINE_REGIMEN'], axis=1, inplace = True)
    #Merge PLAT df with pfs dataframe. Medications_V4.csv
    pfs_PLAT = pfs.merge(pfsMed_plat, how='left', on='PATIENT_ID')
    pfs_PLAT.to_csv('PLAT_Lookheretemp.txt', sep='\t', index=False)
    #Convert objects to numeric
    pfs_PLAT['AGE_AT_PROG_RECUR'] = pd.to_numeric(pfs_PLAT.AGE_AT_PROG_RECUR, errors='coerce')
    pfs_PLAT['AGE_AT_MED_START_Med'] = pd.to_numeric(pfs_PLAT.AGE_AT_MED_START_Med, errors='coerce')
    #calculate pfs_months
    pfs_PLAT['PFS_MONTHS'] = (pfs_PLAT['AGE_AT_PROG_RECUR'] - pfs_PLAT['AGE_AT_MED_START_Med'])*12
    #pfs.drop(['AGE_AT_CURRENT_DISEASE_STATUS', 'CHANGE_OF_TREATMENT', 'CURRENT_DISEASE_STATUS'], axis=1, inplace = True)
    pfs_PLAT.drop(['CURRENT_DISEASE_STATUS'], axis=1, inplace = True)
    #note 'AGE_AT_PROG_RECUR' = pfs 'STOP_DATE'}, inplace=True
    pfs_PLAT.rename(columns={'AGE_AT_PROG_RECUR': 'STOP_DATE'}, inplace=True)
    
    print('******************************************************** KEEP COLUMNS *******************************\n\n')
    print(*pfs)
    print('******************************************************** KEEP COLUMNS *******************************\n\n')
    
    pfs_PLAT = pfs_PLAT[['PATIENT_ID', 'PLATINUM_AGENT', 'PFS_STATUS', 'PFS_MONTHS', 'AGE_AT_MED_START_Med', 'STOP_DATE_Med', 'STOP_DATE', 'SOURCE']]
    pfs_PLAT.rename(columns={'AGE_AT_MED_START_Med': 'AGE_AT_MED_START'}, inplace=True)
    print(*pfs_PLAT)
    #Filter out patients without status
    pfs_PLAT = pfs_PLAT[pfs_PLAT['PFS_STATUS'].isin(['0:DiseaseFree', '1:Progressed'])]
    #remove rows where progression etc happened before med start date.
    pfs_PLAT = pfs_PLAT[(pfs_PLAT['PFS_MONTHS'] > 0)]
    pfs_PLAT.sort_values(by=['STOP_DATE'], inplace = True, ascending=False)
    #Add calculation to select active disease progression within a month of med stop:
    pfs_PLAT['PROGMEDSTOP'] = (pfs_PLAT['STOP_DATE'] - pfs_PLAT['STOP_DATE_Med'])*12
    pfs_PLAT.loc[(pfs_PLAT['SOURCE'].str.contains("OUTCOMES-PROG_RECUR_IND-Active Disease"))&(pfs_PLAT['PROGMEDSTOP'].between(-1, 1)),'PFS_STATUS'] = '1:Progressed'
    pfs_PLAT.drop(['PROGMEDSTOP'], axis=1, inplace = True)
    #pfs.groupby('PATIENT_ID', 'PFS_STATUS').apply(lambda x: x.sort_values('STOP_DATE'), inplace = True, ascending=False)
    print('\n\n************************************************************ Progression FILTER ********************************\n\n')
    #Select progressed lines for filtering
    df_PFS_outcomesProg = pfs_PLAT.loc[pfs['PFS_STATUS'].isin(['1:Progressed'])]
    #df_PFS_outcomesProg.query('STOP_DATE' > 'AGE_AT_MED_STOP')
    #df_PFS_outcomesProg[df_PFS_outcomesProg.STOP_DATE.gt(df_PFS_outcomesProg.AGE_AT_MED_STOP)]
    df_PFS_outcomesProg = df_PFS_outcomesProg.sort_values(by=['PATIENT_ID', 'PFS_MONTHS'], ascending=True)
    ####keep earliest progression
    #df_PFS_outcomesProg.drop_duplicates(subset=['PATIENT_ID'], keep='first', inplace=True)
    df_PFS_outcomesProg.to_csv('PFScombinedProg.txt', sep='\t', index=False)
    print('\n\n************************************************************ Disease Free FILTER ********************************\n\n')

    df_PFS_outcomesDisFree = pfs_PLAT.loc[pfs['PFS_STATUS'].isin(['0:DiseaseFree'])]
    #Eliminate outcomes where the stop date is 
    #df_PFS_outcomesDisFree[df_PFS_outcomesDisFree.STOP_DATE.gt(df_PFS_outcomesDisFree.AGE_AT_MED_STOP)]
    #Sort disease free with longest months on top 
    df_PFS_outcomesDisFree.sort_values(by=['PATIENT_ID', 'PFS_MONTHS'], ascending=False).groupby('PATIENT_ID')
    ####Keep longest disease free PFS
    #df_PFS_outcomesDisFree.drop_duplicates(subset=['PATIENT_ID'], keep='first', inplace=True)
    df_PFS_outcomesDisFree.to_csv('PFScombinedDFree.txt', sep='\t', index=False)

    Merge = [df_PFS_outcomesDisFree, df_PFS_outcomesProg]
    dfOutcomes = pd.concat(Merge)
    dfOutcomes.sort_values(by=['PATIENT_ID'], inplace=True)
    dfOutcomes = dfOutcomes.sort_values(by=['PATIENT_ID', 'PFS_MONTHS'], ascending=False)
    #keep the overal longest PFS shows response or recurrence
    dfOutcomes.drop_duplicates(subset=['PATIENT_ID'], keep='first', inplace=True)
    ######SAVE separate temp files for pfsMed_plat pfsMed_conjug pfsMed_plat = pfsM
    dfOutcomes.rename(columns={'AGENT(S)': 'PLATINUM_AGENT'}, inplace=True)
    dfOutcomes.to_csv('PFS_OutcomesPLAT.txt', sep='\t', index=False)
#********************************* Outcomes with med start merge duplicate dataframes. 

tX = []
FilePatient7 = glob.glob('*_Imaging_V4.csv')
for f in FilePatient7:
    dfstartX = pd.read_csv(f, index_col=None, header=0)
    tX.append(dfstartX)
    pfsImag = pd.concat(tX, axis=0, ignore_index=True)    
#FileSample3 = glob.glob('*_Imaging_V4.csv')
#for f in FileSample3:
#    pfsImag = pd.read_csv(f)
    pfsImag = pfsImag[['AvatarKey', 'AgeAtImageScan', 'ImageLesionGrowth', 'ImagingModalityMethod']]
    print('\n\n************************************************************ Image for Progression ********************************\n\n')
    pfsImag.columns = [camel_to_snake(col) for col in pfsImag.columns]
    pfsImag.rename(columns = str.upper , inplace = True)
    #Rename imaging to account for datain source column
    pfsImag.rename(columns={'AVATAR_KEY': 'PATIENT_ID', 'IMAGING_MODALITY_METHOD': 'SOURCE'}, inplace=True)
    #Patient Filter
    #pfsImag = pfsImag[pfsImag['PATIENT_ID'].isin(['03GCDX0K5N', '06B7Y294XL', '07EF826VTZ', '0C6NW2TC07', '0G5X064EB6', '0GFLCECW75', '0MJYHKZ6AK', '0RA5URUMVK', '0RFVUNS0M6', '0XBPD7NGNT', '15B16LU7OM', '1H50O85Z25', '1HJFO6EFWS', '1IQTFPMLQR', '1K0JG7QIX8', '1K0JG7QIX8', '1Q58ISCSYD', '1RQGB4W7T4', '1YUCT7U163', '274XK9OLE5', '282INULRAK', '28MTFJ7TIV', '2C9ADE4KW6', '2HNKJUZDQR', '2JRJWC9O7M', '2L2ZRV1LZY', '2N04GOV5XK', '2R23NCXQE0', '2UN4T04H08', '37SHF80OE0', '3C2HGYMBKU', '3E15DF0JVG', '3TLQLPNZP8', '41PP6CPRV7', '45PTETBCPX', '45VCVC3MKN', '4AFTCY1HP0', '4B7HJMVOOV', '4FLIDTYACP', '4FLIDTYACP', '4FTJM67GLA', '4OGM5KMQX2', '4S1A75B0CZ', '4UOOLFGTJC', '4UY5RRTZDA', '506HNR0L9Y', '56QVAODIUI', '5CE11AIS9P', '5JWOSJEGEQ', '5UUQPS8Y5W', '5VD97K8TUD', '63JFZJ9HYR', '64F1NAL85E', '66LL4MEA0D', '69NEPW50R9', '6A9ORBLX71', '6DRRKFWIL1', '6JD9VGBWZB', '6PX4GRURXV', '6T5GP2FH7P', '6U163EB0I2', '70HLF9IQQQ', '70P236XSX2', '731EJZUG4W', '73BBJCCFCG', '75HWMQ9THU', '75QKEG7QT5', '761ATF11EM', '78SU3ENGQ3', '78SU3ENGQ3', '797KYKKNTB', '7DQNA9XN5I', '7JT8CIRL35', '7O8FXBKCMP', '7P2075M6NA', '7XZV89ZKWL', '80II7FE25J', '80TN2OH9U2', '8A9JJV56BE', '8C9KM403EG', '8I52ZJY44W', '8IPGJJGS93', '8K3F2U4AEG', '8SSG3I43TG', '8SWRWW6ZR8', '8U1UG57FYC', '8XOABM89XD', '8YEVOS4E0I', '8ZVLI42YTE', '913EQAJLBI', '923EC0J63U', '928IO9LH0Z', '960N6W495Z', '988C2SHC0G', '995JAM1PGX', '9BZRRRQBLI', '9DFUAVXQTU', '9FUJS5X3V0', '9H5WMVOS0K', '9K443LZVX9', '9RQJ8JFILJ', '9SKYB6E4B4', '9ZXCC6ULKK', '9ZXCC6ULKK', 'A049JYM87K', 'A3AM9T7VPE', 'ACLOR9FFS5', 'ALL5VBRLOM', 'AP93ZGZX2W', 'AVMUPL2GS0', 'B3D07F2E2S', 'BCZK5PTSG9', 'BK6DKQTSI5', 'BTHE1B9DZ2', 'BU8GCS3MYT', 'BV2N69FYI6', 'BV2N69FYI6', 'BXXKG0HZIK', 'C5A7KUNVCS', 'CADAN6NOJ4', 'CGDPFH8IBW', 'CK9QL0A5ET', 'CM03MT997O', 'CM03MT997O', 'CM03MT997O', 'CTGHR4UZ18', 'CW0ZVIVU0P', 'CZE113UZIR', 'D0WP0B97Q7', 'D2C5F4UMXY', 'D7KAPTXFOG', 'DCPKIAYXDB', 'DQF8CHSKPG', 'DYKRF4KH3H', 'E2H89HF5KZ', 'E4H964SOU9', 'EF7AOZI79O', 'EHOVA8D07G', 'EKSOBYXBVG', 'EKXLPR0080', 'EMVBZPMWST', 'EOZ9LELKBM', 'EQNGL8N9W7', 'EU4WQTA1KL', 'F1XUR9R1RK', 'F2TQDDM3HH', 'F905RYZFIO', 'FNVA6L545Q', 'FOTCBZ3OMJ', 'FPLKDOPOH3', 'FVZV8LMNI5', 'FXGR65TB4V', 'FXSR99KNNT', 'GA7GNCV7KM', 'GF6S6LZVF5', 'GH6KKNEPWX', 'GIJZDVWC2Q', 'GNJVFE3HRM', 'GNXLDBYWN1', 'GOYPZUZNAW', 'GPLOBS1HRX', 'GU9WE974XH', 'GWNOJYMTTB', 'HG7EO5096W', 'HUCN0XE91R', 'HWECMUCW0Z', 'HWIG0FSC6Z', 'HXKS3QWP6F', 'HY2ZUQI0SD', 'I2EM49R0AO', 'I65GB8ZEJ2', 'IMH5IF30JJ', 'IW3OCCZBPJ', 'IW7K4RUJG9', 'IXEPJTR5B9', 'J0D6I815XU', 'J15Q9C8MIE', 'J3C9QS0P78', 'J60EKEW4OF', 'J81JUVHO9L', 'J9NHJFS8IN', 'JBLZPPZSFH', 'JG2Z3K6UC5', 'JHBPEVZ48F', 'JHZSJNIKIP', 'JJGV4LPENY', 'JLGA08C3K5', 'JM3MBDWG5H', 'JNXYDUFEED', 'JPSK03CCDH', 'JPSK03CCDH', 'JQHVTFHS1X', 'JRC4RTZNTK', 'JSWZX87JYV', 'JZKKZ1OTQJ', 'K1DVC6O9C7', 'K60LO6G1FG', 'K8DFO8FTT1', 'K8OFSLEHIX', 'KAN1H1VWQG', 'KC8HPZPX7T', 'KPCHILTX70', 'KSTUVZCJIF', 'KW650F0557', 'L410MYBLDN', 'L6RMSU0BXG', 'L8HCOD7F64', 'LATEJ509E1', 'LLXVLRWEAE', 'LM4ZZPQQFO', 'LS6UIWKFEJ', 'LX5K9JZKLH', 'LXELBL78SE', 'M1D0DA5RWT', 'M9G3Y6A4X7', 'MG1SUQR5K8', 'MIVJ67B9F5', 'MJYLGG6W0K', 'MRZW12GMST', 'MYG3XG6SAQ', 'N4K1N7UUMM', 'N675BI25Y6', 'N8USV2HKJ6', 'NAMDA2GZXA', 'NAR8MDP01R', 'NDV4VM27O6', 'NE17FOVUC5', 'NK65MD354P', 'NOD926GZ41', 'NPUNSPR11V', 'O8RFN0K0CC', 'OA1K7SKWYW', 'OFMS9RN5W8', 'OFOSHGCQE4', 'OJQ5IBJCDB', 'ONHYTARPID', 'OPI88VIT01', 'OQYBJZXP4D', 'ORL6LF2IMD', 'OWNR0K9J5O', 'OXKYPBVO0Q', 'OXKYPBVO0Q', 'P0YOHVOGGY', 'P6F7GJOYOS', 'PCM0C3L6G1', 'PGVXTE18CS', 'PJ316U0AO6', 'PQDMEIJU8A', 'PTN0DBFVOH', 'Q1LFWVL6GK', 'Q1TM0Z8WWH', 'Q3XX8RO9CY', 'QB511DI6QC', 'QD6U65T2FH', 'QDTQI9DG2U', 'QJ80IGAAAA', 'QO27BH1WTV', 'QUEJJTOM01', 'QVEQP2FBPS', 'QVLKW6YJ9S', 'QZ68W571DZ', 'R0JJZPXHE6', 'R9EFJQNQD7', 'RC5DO9WG9A', 'RDNMW7PS5V', 'RL11WPPBNV', 'RN9YUHZDTE', 'RRWTZZSW3F', 'S414XWRR6V', 'S5AEWMYN2T', 'S5R5RW4L05', 'S7GI3EL7N1', 'S9CLR0UMQD', 'S9GIV41A0F', 'SF953NJFYV', 'SQV8WBN78P', 'STZJJQ7LKL', 'SUYGY0MFVB', 'SW0L8B0LS4', 'T3BB8UT9KP', 'T653MJNWUY', 'T8IWENLDFE', 'TEEDQQPZOI', 'TFJ3UKACLL', 'TG4OSS40W2', 'TGCXP2RE55', 'TTMBZYSHU8', 'TTN1Q2AA35', 'TTN1Q2AA35', 'TTUKMQJUXT', 'TX7XF0XMBL', 'TXJ9VXE4QC', 'TYUGO1UOO3', 'U0JP7MO0CU', 'U1KHQ8NI5I', 'U3YN74R7B9', 'U6SEJCBUKI', 'U86PUE8URG', 'UB8PDQ9GI1', 'UD8WP9Y59H', 'UFKCGRWOE6', 'ULH5Q1UQZ9', 'UMZ44YQOYC', 'UMZ44YQOYC', 'UOHLJNY3WH', 'UOWQ0WV182', 'UYBI9N0WBE', 'V3V21CTLZR', 'VE1W0KD4IA', 'VF0HIW4X63', 'VI0WMPQNAD', 'VIMLBEI2ZQ', 'VKOQ16NG0M', 'VKOQ16NG0M', 'VMKJ7V6SRD', 'VMKJ7V6SRD', 'VQS7PT887V', 'VV8VUVF9XA', 'W3W1ROPU95', 'W6XNV1LQB1', 'W7UAIA7T0S', 'W80YCK8SOZ', 'WCYEL3APQW', 'WIYGVFIKGG', 'WK78CI5LDO', 'WMJ41DPJ3Y', 'WPYHWCFES9', 'WR6V5879WH', 'WXJA24FUA7', 'WYMINM1LOF', 'X0W1ONC8B2', 'X1MS489ILB', 'X2AGMD6P59', 'X4HHLWE7SC', 'X4JDSBA0MD', 'X4JDSBA0MD', 'X5A3YUL1OI', 'X5WMJ9MUOV', 'X7566XEPLR', 'X8NG9V0HJ6', 'XASH2D2L7T', 'XG77M6C95Z', 'XLDFSI23XN', 'XTIR5SW5E8', 'Y0S85J0AZR', 'Y1DJA21YHZ', 'Y4HVM00UR9', 'YBIJBXII7J', 'YF6JV0C9MC', 'YP49T16HM1', 'YV85IDK7KQ', 'YY88X645NE', 'Z0075688AQ', 'Z0RE2PT1R9', 'Z1WO6PGN3E', 'Z9478QITNI', 'ZA6WQHDGBM', 'ZBG7MKFOZ1', 'ZC60EAOUMH', 'ZE7RW2IJMR', 'ZHRR08KKVM', 'ZIN73MJOAB', 'ZINL8265JQ', 'ZKG95SXJUU', 'ZKZ02PIDIU', 'ZM2Y4LL5QP', 'ZQRZQA3DZT', 'ZW84KJ7LJA'])]

    pfsImag = pfsImag[pfsImag['IMAGE_LESION_GROWTH'].isin(['Yes'])]
    pfsImag['SOURCE'] = 'IMAGING-' +  pfsImag['SOURCE']  
    #Build/change PFS status
    pfsImag.rename(columns={'SOURCE': 'IMAGE_INDICATION', 'IMAGE_LESION_GROWTH': 'PFS_STATUS'}, inplace=True)
    pfsImag['PFS_STATUS'] = pfsImag['PFS_STATUS'].str.replace('Yes', '1:Progressed')
    #Drop duplicate progressions etc for respective start dates later on - duplicate below thru **************
    pfs_platImag = pfsImag
    pfs_platImag.to_csv('pfs_platImag.txt', sep='\t', index=False)
    pfs_conjugImag = pfsImag
    pfs_conjugImag.to_csv('pfs_conjugImag.txt', sep='\t', index=False)



#**************************************************   DUPlicate DF for other drugs    
    #***** PLAT - specific start dates
    pfsMed3 = pd.read_csv('MedicationsPFStempPLAT.txt', sep='\t')
    pfsMed3.drop(['PFS_STATUS', 'PFS_MONTHS'], axis=1, inplace = True)
    pfsImagPLAT = pfsImag.merge(pfsMed3, how='left', on='PATIENT_ID')
    print(*pfsImagPLAT)
    
    #pfsImagPLAT.drop(['PFS_STATUS', 'PFS_MONTHS', 'STOP_DATE_Med'], axis=1, inplace = True)
    #Drop patients with no PLAT 
    pfsImagPLAT = pfsImagPLAT[pfsImagPLAT['AGENT(S)'].isin(['Cisplatin', 'Carboplatin', 'Oxaliplatin'])]
    pfsImagPLAT['AGE_AT_IMAGE_SCAN'] = pfsImag['AGE_AT_IMAGE_SCAN'].str.replace('Age 90 or older', '91')
    pfsImagPLAT['AGE_AT_IMAGE_SCAN'] = pd.to_numeric(pfsImag.AGE_AT_IMAGE_SCAN, errors='coerce')
    pfsImagPLAT = pfsImagPLAT.sort_values(by=['PATIENT_ID', 'AGE_AT_IMAGE_SCAN'], ascending=True)
    #remove rows where image-based progression etc happened before med start date.
    pfsImagPLAT['TempFilter'] = (pfsImagPLAT['AGE_AT_IMAGE_SCAN'] - pfsImagPLAT['AGE_AT_MED_START_Med'])*12
    pfsImagPLAT['TempFilter'] = pd.to_numeric(pfsImagPLAT.TempFilter, errors='coerce')
    pfsImagPLAT['TempFilter'] = pfsImagPLAT['TempFilter'].fillna(0)
    pfsImagPLAT.to_csv('pfsImag.txt', sep='\t', index=False)
    #keep first PFS event and those from after med start
    pfsImagPLAT = pfsImagPLAT[(pfsImagPLAT['TempFilter'] > 0)]
    pfsImagPLAT = pfsImagPLAT.sort_values(by=['PATIENT_ID', 'TempFilter'], ascending=True)
    #pfsImagPLAT.drop_duplicates(subset=['PATIENT_ID'], keep='first', inplace=True)
    pfsImagPLAT.drop(['TempFilter'], axis=1, inplace = True)
    
    print('\n\n********************************************************\n\n')
    print(*pfsImagPLAT)
    print('\n\n********************************************************\n\n')
    
    #pfsImagPLAT.to_csv('PFS_IMAGsort.txt', sep='\t', index=False)
    pfsImagPLAT.drop_duplicates(subset=['PATIENT_ID', 'AGE_AT_MED_START_Med'], keep='first', inplace=True)
    #Calculate PFS
    pfsImagPLAT['AGE_AT_IMAGE_SCAN'] = pd.to_numeric(pfsImagPLAT.AGE_AT_IMAGE_SCAN, errors='coerce')
    pfsImagPLAT['AGE_AT_MED_START_Med'] = pd.to_numeric(pfsImagPLAT.AGE_AT_MED_START_Med, errors='coerce')
    pfsImagPLAT['PFS_MONTHS'] = (pfsImagPLAT['AGE_AT_IMAGE_SCAN'] - pfsImagPLAT['AGE_AT_MED_START_Med'])*12
    pfsImagPLAT.to_csv('PFS_IMAGsort.txt', sep='\t', index=False)

    #Cleanup
    pfsImagPLAT.rename(columns={'AGE_AT_IMAGE_SCAN': 'STOP_DATE'}, inplace=True)
    #pfsImagPLAT['STOP_DATE_OTHER'] =  pfsImagPLAT['STOP_DATE']
    pfsImagPLAT.rename(columns={'AGE_AT_MED_START': 'AGE_AT_MED_START_Imag'}, inplace=True)
    pfsImagPLAT = pfsImagPLAT[['PATIENT_ID', 'AGENT(S)', 'PFS_STATUS', 'PFS_MONTHS', 'AGE_AT_MED_START_Med', 'STOP_DATE', 'IMAGE_INDICATION']]
    #cleanup Imaging
    pfsImagPLAT.to_csv('PFS_IMAG_PLAT.txt', sep='\t', index=False)
#****************************************    
tXI = []
FilePatient8 = glob.glob('*_MetastaticDisease_V4.csv')
for f in FilePatient8:
    dfstartXI = pd.read_csv(f, index_col=None, header=0)
    tXI.append(dfstartXI)
    pfsMet = pd.concat(tXI, axis=0, ignore_index=True)    
#FileSample4 = glob.glob('*_MetastaticDisease_V4.csv')
#for f in FileSample4:
#    pfsMet = pd.read_csv(f)
    pfsMet = pfsMet[['AvatarKey', 'MetastaticDiseaseInd', 'AgeAtMetastaticSite']]
    print('\n\n************************************************************ Met for Progression ********************************\n\n')
    pfsMet.columns = [camel_to_snake(col) for col in pfsMet.columns]
    pfsMet.rename(columns = str.upper , inplace = True)
    pfsMet.rename(columns={'AVATAR_KEY': 'PATIENT_ID'}, inplace=True)
    #Patient Filter
    #pfsMet = pfsMet[pfsMet['PATIENT_ID'].isin(['03GCDX0K5N', '06B7Y294XL', '07EF826VTZ', '0C6NW2TC07', '0G5X064EB6', '0GFLCECW75', '0MJYHKZ6AK', '0RA5URUMVK', '0RFVUNS0M6', '0XBPD7NGNT', '15B16LU7OM', '1H50O85Z25', '1HJFO6EFWS', '1IQTFPMLQR', '1K0JG7QIX8', '1K0JG7QIX8', '1Q58ISCSYD', '1RQGB4W7T4', '1YUCT7U163', '274XK9OLE5', '282INULRAK', '28MTFJ7TIV', '2C9ADE4KW6', '2HNKJUZDQR', '2JRJWC9O7M', '2L2ZRV1LZY', '2N04GOV5XK', '2R23NCXQE0', '2UN4T04H08', '37SHF80OE0', '3C2HGYMBKU', '3E15DF0JVG', '3TLQLPNZP8', '41PP6CPRV7', '45PTETBCPX', '45VCVC3MKN', '4AFTCY1HP0', '4B7HJMVOOV', '4FLIDTYACP', '4FLIDTYACP', '4FTJM67GLA', '4OGM5KMQX2', '4S1A75B0CZ', '4UOOLFGTJC', '4UY5RRTZDA', '506HNR0L9Y', '56QVAODIUI', '5CE11AIS9P', '5JWOSJEGEQ', '5UUQPS8Y5W', '5VD97K8TUD', '63JFZJ9HYR', '64F1NAL85E', '66LL4MEA0D', '69NEPW50R9', '6A9ORBLX71', '6DRRKFWIL1', '6JD9VGBWZB', '6PX4GRURXV', '6T5GP2FH7P', '6U163EB0I2', '70HLF9IQQQ', '70P236XSX2', '731EJZUG4W', '73BBJCCFCG', '75HWMQ9THU', '75QKEG7QT5', '761ATF11EM', '78SU3ENGQ3', '78SU3ENGQ3', '797KYKKNTB', '7DQNA9XN5I', '7JT8CIRL35', '7O8FXBKCMP', '7P2075M6NA', '7XZV89ZKWL', '80II7FE25J', '80TN2OH9U2', '8A9JJV56BE', '8C9KM403EG', '8I52ZJY44W', '8IPGJJGS93', '8K3F2U4AEG', '8SSG3I43TG', '8SWRWW6ZR8', '8U1UG57FYC', '8XOABM89XD', '8YEVOS4E0I', '8ZVLI42YTE', '913EQAJLBI', '923EC0J63U', '928IO9LH0Z', '960N6W495Z', '988C2SHC0G', '995JAM1PGX', '9BZRRRQBLI', '9DFUAVXQTU', '9FUJS5X3V0', '9H5WMVOS0K', '9K443LZVX9', '9RQJ8JFILJ', '9SKYB6E4B4', '9ZXCC6ULKK', '9ZXCC6ULKK', 'A049JYM87K', 'A3AM9T7VPE', 'ACLOR9FFS5', 'ALL5VBRLOM', 'AP93ZGZX2W', 'AVMUPL2GS0', 'B3D07F2E2S', 'BCZK5PTSG9', 'BK6DKQTSI5', 'BTHE1B9DZ2', 'BU8GCS3MYT', 'BV2N69FYI6', 'BV2N69FYI6', 'BXXKG0HZIK', 'C5A7KUNVCS', 'CADAN6NOJ4', 'CGDPFH8IBW', 'CK9QL0A5ET', 'CM03MT997O', 'CM03MT997O', 'CM03MT997O', 'CTGHR4UZ18', 'CW0ZVIVU0P', 'CZE113UZIR', 'D0WP0B97Q7', 'D2C5F4UMXY', 'D7KAPTXFOG', 'DCPKIAYXDB', 'DQF8CHSKPG', 'DYKRF4KH3H', 'E2H89HF5KZ', 'E4H964SOU9', 'EF7AOZI79O', 'EHOVA8D07G', 'EKSOBYXBVG', 'EKXLPR0080', 'EMVBZPMWST', 'EOZ9LELKBM', 'EQNGL8N9W7', 'EU4WQTA1KL', 'F1XUR9R1RK', 'F2TQDDM3HH', 'F905RYZFIO', 'FNVA6L545Q', 'FOTCBZ3OMJ', 'FPLKDOPOH3', 'FVZV8LMNI5', 'FXGR65TB4V', 'FXSR99KNNT', 'GA7GNCV7KM', 'GF6S6LZVF5', 'GH6KKNEPWX', 'GIJZDVWC2Q', 'GNJVFE3HRM', 'GNXLDBYWN1', 'GOYPZUZNAW', 'GPLOBS1HRX', 'GU9WE974XH', 'GWNOJYMTTB', 'HG7EO5096W', 'HUCN0XE91R', 'HWECMUCW0Z', 'HWIG0FSC6Z', 'HXKS3QWP6F', 'HY2ZUQI0SD', 'I2EM49R0AO', 'I65GB8ZEJ2', 'IMH5IF30JJ', 'IW3OCCZBPJ', 'IW7K4RUJG9', 'IXEPJTR5B9', 'J0D6I815XU', 'J15Q9C8MIE', 'J3C9QS0P78', 'J60EKEW4OF', 'J81JUVHO9L', 'J9NHJFS8IN', 'JBLZPPZSFH', 'JG2Z3K6UC5', 'JHBPEVZ48F', 'JHZSJNIKIP', 'JJGV4LPENY', 'JLGA08C3K5', 'JM3MBDWG5H', 'JNXYDUFEED', 'JPSK03CCDH', 'JPSK03CCDH', 'JQHVTFHS1X', 'JRC4RTZNTK', 'JSWZX87JYV', 'JZKKZ1OTQJ', 'K1DVC6O9C7', 'K60LO6G1FG', 'K8DFO8FTT1', 'K8OFSLEHIX', 'KAN1H1VWQG', 'KC8HPZPX7T', 'KPCHILTX70', 'KSTUVZCJIF', 'KW650F0557', 'L410MYBLDN', 'L6RMSU0BXG', 'L8HCOD7F64', 'LATEJ509E1', 'LLXVLRWEAE', 'LM4ZZPQQFO', 'LS6UIWKFEJ', 'LX5K9JZKLH', 'LXELBL78SE', 'M1D0DA5RWT', 'M9G3Y6A4X7', 'MG1SUQR5K8', 'MIVJ67B9F5', 'MJYLGG6W0K', 'MRZW12GMST', 'MYG3XG6SAQ', 'N4K1N7UUMM', 'N675BI25Y6', 'N8USV2HKJ6', 'NAMDA2GZXA', 'NAR8MDP01R', 'NDV4VM27O6', 'NE17FOVUC5', 'NK65MD354P', 'NOD926GZ41', 'NPUNSPR11V', 'O8RFN0K0CC', 'OA1K7SKWYW', 'OFMS9RN5W8', 'OFOSHGCQE4', 'OJQ5IBJCDB', 'ONHYTARPID', 'OPI88VIT01', 'OQYBJZXP4D', 'ORL6LF2IMD', 'OWNR0K9J5O', 'OXKYPBVO0Q', 'OXKYPBVO0Q', 'P0YOHVOGGY', 'P6F7GJOYOS', 'PCM0C3L6G1', 'PGVXTE18CS', 'PJ316U0AO6', 'PQDMEIJU8A', 'PTN0DBFVOH', 'Q1LFWVL6GK', 'Q1TM0Z8WWH', 'Q3XX8RO9CY', 'QB511DI6QC', 'QD6U65T2FH', 'QDTQI9DG2U', 'QJ80IGAAAA', 'QO27BH1WTV', 'QUEJJTOM01', 'QVEQP2FBPS', 'QVLKW6YJ9S', 'QZ68W571DZ', 'R0JJZPXHE6', 'R9EFJQNQD7', 'RC5DO9WG9A', 'RDNMW7PS5V', 'RL11WPPBNV', 'RN9YUHZDTE', 'RRWTZZSW3F', 'S414XWRR6V', 'S5AEWMYN2T', 'S5R5RW4L05', 'S7GI3EL7N1', 'S9CLR0UMQD', 'S9GIV41A0F', 'SF953NJFYV', 'SQV8WBN78P', 'STZJJQ7LKL', 'SUYGY0MFVB', 'SW0L8B0LS4', 'T3BB8UT9KP', 'T653MJNWUY', 'T8IWENLDFE', 'TEEDQQPZOI', 'TFJ3UKACLL', 'TG4OSS40W2', 'TGCXP2RE55', 'TTMBZYSHU8', 'TTN1Q2AA35', 'TTN1Q2AA35', 'TTUKMQJUXT', 'TX7XF0XMBL', 'TXJ9VXE4QC', 'TYUGO1UOO3', 'U0JP7MO0CU', 'U1KHQ8NI5I', 'U3YN74R7B9', 'U6SEJCBUKI', 'U86PUE8URG', 'UB8PDQ9GI1', 'UD8WP9Y59H', 'UFKCGRWOE6', 'ULH5Q1UQZ9', 'UMZ44YQOYC', 'UMZ44YQOYC', 'UOHLJNY3WH', 'UOWQ0WV182', 'UYBI9N0WBE', 'V3V21CTLZR', 'VE1W0KD4IA', 'VF0HIW4X63', 'VI0WMPQNAD', 'VIMLBEI2ZQ', 'VKOQ16NG0M', 'VKOQ16NG0M', 'VMKJ7V6SRD', 'VMKJ7V6SRD', 'VQS7PT887V', 'VV8VUVF9XA', 'W3W1ROPU95', 'W6XNV1LQB1', 'W7UAIA7T0S', 'W80YCK8SOZ', 'WCYEL3APQW', 'WIYGVFIKGG', 'WK78CI5LDO', 'WMJ41DPJ3Y', 'WPYHWCFES9', 'WR6V5879WH', 'WXJA24FUA7', 'WYMINM1LOF', 'X0W1ONC8B2', 'X1MS489ILB', 'X2AGMD6P59', 'X4HHLWE7SC', 'X4JDSBA0MD', 'X4JDSBA0MD', 'X5A3YUL1OI', 'X5WMJ9MUOV', 'X7566XEPLR', 'X8NG9V0HJ6', 'XASH2D2L7T', 'XG77M6C95Z', 'XLDFSI23XN', 'XTIR5SW5E8', 'Y0S85J0AZR', 'Y1DJA21YHZ', 'Y4HVM00UR9', 'YBIJBXII7J', 'YF6JV0C9MC', 'YP49T16HM1', 'YV85IDK7KQ', 'YY88X645NE', 'Z0075688AQ', 'Z0RE2PT1R9', 'Z1WO6PGN3E', 'Z9478QITNI', 'ZA6WQHDGBM', 'ZBG7MKFOZ1', 'ZC60EAOUMH', 'ZE7RW2IJMR', 'ZHRR08KKVM', 'ZIN73MJOAB', 'ZINL8265JQ', 'ZKG95SXJUU', 'ZKZ02PIDIU', 'ZM2Y4LL5QP', 'ZQRZQA3DZT', 'ZW84KJ7LJA'])]

    print(*pfsMet)
    #Take all available mets by filtering those without Mets - used double negative as comments were variable for "yes/positive"
    pfsMet = pfsMet[~pfsMet['METASTATIC_DISEASE_IND'].isin(['No'])]
    #drop dulicate mets
    pfsMet.drop_duplicates(subset=['PATIENT_ID', 'AGE_AT_METASTATIC_SITE'], keep='first', inplace=True)
    pfsMet.sort_values(by=['AGE_AT_METASTATIC_SITE'], ascending=False).groupby('PATIENT_ID')
    #Drop samples with no dates
    pfsMet = pfsMet[~pfsMet['AGE_AT_METASTATIC_SITE'].isin(['Unknown/Not Applicable'])]
    #BUILD duplicated DFs for other drug calsses. 
    pfsMet_plat = pfsMet
    pfsMet_plat.to_csv('pfsMet_plat.txt', sep='\t', index=False)
    pfsMet_conjug = pfsMet
    pfsMet_conjug.to_csv('pfsMet_conjug.txt', sep='\t', index=False)

#MET DUPLICATE dataframe ************************************               
    #Merge in medication/PLAT start dates
    pfsMed4 = pd.read_csv('MedicationsPFStempPLAT.txt', sep='\t')
    pfsMetPLAT = pfsMet.merge(pfsMed3, how='left', on='PATIENT_ID')
    pfsMetPLAT.drop_duplicates(subset=['PATIENT_ID', 'AGE_AT_METASTATIC_SITE'], inplace=True)
    pfsMetPLAT.to_csv('pfsMetPLATTemperr.txt', sep='\t', index=False)
    #Calculate MFS
    pfsMetPLAT['AGE_AT_METASTATIC_SITE'] = pd.to_numeric(pfsMetPLAT.AGE_AT_METASTATIC_SITE, errors='coerce')
    pfsMetPLAT['AGE_AT_MED_START_Med'] = pd.to_numeric(pfsMetPLAT.AGE_AT_MED_START_Med, errors='coerce')
    pfsMetPLAT['PFS_MONTHS'] = (pfsMetPLAT['AGE_AT_METASTATIC_SITE'] - pfsMetPLAT['AGE_AT_MED_START_Med'])*12
    #clear all met events before med start
    pfsMetPLAT = pfsMetPLAT[(pfsMetPLAT['PFS_MONTHS'] > 0)]
    #Drop patients with no PLAT 
    pfsMetPLAT = pfsMetPLAT[pfsMetPLAT['AGENT(S)'].isin(['Cisplatin', 'Carboplatin', 'Oxaliplatin'])]
    #Cleanup
    pfsMetPLAT['PFS_STATUS'] = '1:Progressed'
    pfsMetPLAT.rename(columns={'METASTATIC_DISEASE_IND': 'SOURCE', 'AGE_AT_METASTATIC_SITE': 'STOP_DATE'}, inplace=True)
    pfsMetPLAT['STOP_DATE_OTHER'] =  pfsMetPLAT['STOP_DATE']
    pfsMetPLAT['SOURCE'] = pfsMetPLAT['SOURCE'].str.replace('Yes', 'METASTATIC_DISEASE-METASTATIC_DISEASE_IND')
    pfsMetPLAT.to_csv('pfsMetPLAT.txt', sep='\t', index=False)
    pfsMetPLAT = pfsMetPLAT[['PATIENT_ID', 'AGENT(S)', 'PFS_STATUS', 'PFS_MONTHS', 'AGE_AT_MED_START_Med', 'STOP_DATE_Med', 'STOP_DATE', 'SOURCE', 'STOP_DATE_OTHER']]
    pfsMetPLAT.rename(columns={'AGE_AT_MED_START': 'AGE_AT_MED_START_Met', 'AGE_AT_MED_STOP': 'AGE_AT_MED_STOP_Met'}, inplace=True)
    pfsMetPLAT = pfsMetPLAT[(pfsMetPLAT['PFS_MONTHS'] > 0)]
    pfsMetPLAT.rename(columns={'SOURCE': 'MET_INDICATION'}, inplace=True)
    pfsMetPLAT.to_csv('pfsMetPLAT.txt', sep='\t', index=False)
#MET DUPLICATE dataframe ************************************

dfOutcomes.rename(columns={'SOURCE': 'OUTCOMES_INDICATION'}, inplace=True)
#cleanup meds
pfsMed_plat = pd.read_csv('MedicationsPFStempPLAT.txt', sep='\t')
pfsMed_plat['PFS_MONTHS'] = pfsMed_plat[pfsMed_plat['PFS_MONTHS'] > 0]['PFS_MONTHS']
pfsMed_plat['PFS_MONTHS'].replace('', np.nan, inplace=True)
pfsMed_plat = pfsMed_plat[pfsMed_plat['PFS_MONTHS'].notna()]
pfsMed_plat.rename(columns={'SOURCE': 'MEDICATION_INDICATION'}, inplace=True)

#Merge pfs, pfsMed, pfsImag, pfsMet
df_PFSTransposed = dfOutcomes.merge(pfsMed_plat, how='outer', on='PATIENT_ID').merge(pfsImagPLAT, how='outer', on='PATIENT_ID').merge(pfsMetPLAT, how='outer', on='PATIENT_ID')
df_PFSTransposed.to_csv('PFScombinedTemp.txt', sep='\t', index=False)
#Sort dataframe to filter simple rows without multiple PFS indications
df_PFSTransposed.sort_values(by=['MEDICATION_INDICATION', 'IMAGE_INDICATION', 'MET_INDICATION'], inplace=True)
df_PFSTransposed.to_csv('PFScombined.txt', sep='\t', index=False)

####MAKE PROGRESSION FILE and merge all progression data
print('\n\n************************************************************  Progression PLAT Merge ********************************\n\n')
pfsMed_plat.rename(columns={'AGENT(S)': 'PLATINUM_AGENT', 'AGE_AT_MED_START_Med': 'AGE_AT_MED_START', 'STOP_DATE_Med': 'STOP_DATE', 'MEDICATION_INDICATION': 'INDICATION'}, inplace=True)
print(*pfsMed_plat)
pfsMed_plat.to_csv('pfsMed_plat.FInal.txt', sep='\t', index=False)
pfsMetPLAT.rename(columns={'AGENT(S)': 'PLATINUM_AGENT', 'AGE_AT_MED_START_Med': 'AGE_AT_MED_START', 'MET_INDICATION': 'INDICATION'}, inplace=True)
#pfsMet.drop(['AGE_AT_METASTATIC_SITE'], axis=1, inplace = True)
print(*pfsMetPLAT)
pfsMetPLAT.to_csv('pfsMetPLAT.FInal.txt', sep='\t', index=False)
pfsImagPLAT.rename(columns={'AGENT(S)': 'PLATINUM_AGENT', 'AGE_AT_MED_START_Med': 'AGE_AT_MED_START', 'STOP_DATE_Med': 'STOP_DATE', 'IMAGE_INDICATION': 'INDICATION'}, inplace=True)
#pfsImag.drop(['AGE_AT_IMAGE_SCAN'], axis=1, inplace = True)
print(*pfsImagPLAT)
pfsImagPLAT.to_csv('pfsImagPLAT.FInal.txt', sep='\t', index=False)
df_PFS_outcomesProg.drop(['STOP_DATE_Med'], axis=1, inplace = True)
df_PFS_outcomesProg.rename(columns={'SOURCE': 'INDICATION'}, inplace=True)
print(*df_PFS_outcomesProg)
#concatenate and take earliest. 
Progression = [df_PFS_outcomesProg, pfsMed_plat, pfsImagPLAT, pfsMetPLAT]
dfProgFinal = pd.concat(Progression)
dfProgFinal = dfProgFinal.sort_values(by=['PATIENT_ID', 'PFS_MONTHS'], ascending=False)
#dfProgFinal['PLATINUM_AGENT'] = dfProgFinal['PLATINUM_AGENT'].replace('', pd.NA).fillna(dfProgFinal['AGENT(S)'])
dfProgFinal.to_csv('dfProgFinal.txt', sep='\t', index=False)
####Keep longest disease free PFS
#dfProgFinal.drop_duplicates(subset=['PATIENT_ID'], keep='first', inplace=True)
print(*dfProgFinal)
#Bring in disease free
df_PFS_outcomesDisFree.drop(['STOP_DATE_Med'], axis=1, inplace = True)
df_PFS_outcomesDisFree.rename(columns={'SOURCE': 'INDICATION'}, inplace=True)
print(*df_PFS_outcomesDisFree)
Outcomes = [df_PFS_outcomesDisFree, dfProgFinal]
pfsFinal = pd.concat(Outcomes)
#Sort disease free with longest months on top 
pfsFinal = pfsFinal.sort_values(by=['PATIENT_ID', 'PFS_MONTHS'], ascending=False)
#pfsFinal = pfsFinal.sort_values(by=['PFS_MONTHS'], ascending=False)
####Keep longest disease free PFS
#pfsFinal.drop_duplicates(subset=['PATIENT_ID'], keep='first', inplace=True)

#remove patients with age estimations
pfsFinal = pfsFinal[~pfsFinal['STOP_DATE'].isin([91])]
pfsFinal.rename(columns={'INDICATION': 'PFS_INDICATION_PLAT'}, inplace=True)


pfsFinal.to_csv('PFS_Final.txt', sep='\t')
#Split PFS into separate docs for filtering logic grabbing unique patientIDs for temp file
for patient in pfsFinal['PATIENT_ID'].unique():
    sub_dff = pfsFinal[pfsFinal['PATIENT_ID'] == patient]
    if sub_dff['PFS_STATUS'].isin(['1:Progressed']).all():
    #if sub_dff['PFS_STATUS'].eq(sub_dff['PFS_STATUS'].iloc[0]).all():
    #if len(np.unique(sub_dff.PFS_STATUS))=='1:Progressed':
        sub_dff = sub_dff.drop_duplicates(subset=['PATIENT_ID'], keep='last')
    if sub_dff['PFS_STATUS'].isin(['0:DiseaseFree']).all():
        sub_dff = sub_dff.drop_duplicates(subset=['PATIENT_ID'], keep='first')
    #if sub_dff[sub_dff.PFS_STATUS == (['1:Progressed'])].iloc[0]:
        #sub_dff = sub_dff.drop_duplicates(subset=['PATIENT_ID'], keep='first')
    sub_dff.to_csv(f'{patient}_PLAT_data.txt', sep='\t', index=False)
    
filenames = glob.glob('*PLAT_data.txt')

dff2 = (pd.read_csv(f, index_col=None, header=0, nrows=1, sep='\t') for f in filenames)
concatenated_dfPLAT = pd.concat(dff2, ignore_index=True)
concatenated_dfPLAT.to_csv('PFS_FinalTrimPLAT.txt', sep='\t', index=False)
concatenated_dfPLAT.rename(columns={'PFS_STATUS': 'PFS_STATUS_PLAT', 'PFS_MONTHS': 'PFS_MONTHS_PLAT', 'AGE_AT_MED_START': 'AGE_AT_MED_START_PLAT', 'STOP_DATE_Med': 'STOP_DATE_Med_PLAT', 'MED_LINE_REGIMEN': 'MED_LINE_REGIMEN_PLAT'}, inplace=True)
concatenated_dfPLAT.drop(['STOP_DATE_OTHER'], axis=1, inplace = True)

print('\n\n**********************************   MERGE ALL PFS data      ****************************\n\n\n\n')

#Merge to PatientCDE
PatientCDE = PatientCDE.merge(concatenated_dfPLAT, how='left', on='PATIENT_ID')
PatientCDE.drop_duplicates(subset=['PATIENT_ID'], keep='first', inplace=True)
PatientCDE.drop(['STOP_DATE_y', 'STOP_DATE_x'], axis=1, inplace = True)
PatientCDE.rename(columns={'PLATINUM_AGENT_y': 'PLATINUM_AGENT'}, inplace=True)



######################  FINAL HEADER CORRECTION #################
#add cBio specific header meta data
PatientCDEmeta = PatientCDE.dtypes
PatientCDE = PatientCDE.append(PatientCDEmeta, ignore_index=True)
cols = PatientCDE.columns.tolist()
PatientCDE.loc[-1] = cols
PatientCDE = PatientCDE.reindex(np.roll(PatientCDE.index, shift=1))
PatientCDE.loc[-1] = ''  # adding a row
PatientCDE = PatientCDE.reindex(np.roll(PatientCDE.index, shift=1))
PatientCDE.to_csv('PatientCDE.txt', sep='\t', index=False)