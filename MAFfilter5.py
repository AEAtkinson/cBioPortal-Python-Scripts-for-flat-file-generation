#!/usr/bin/env python3

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan  3 08:29:37 2019 through 2/24/2019

@author: aaronatkinson
"""
#import packages
import glob, pandas as pd, numpy as np

#Assign directory to change or us within directory
path = "/Users/aaronatkinson/Desktop/cBioSomaticMafs072020"


#create filelist
read_files = glob.glob('*.Anno.Hg38.anno.filt.vcf.maf')

#Compile master MAF file of all MAFs. 
with open("AvatarSomaticACMG_Master.maf", "wb") as outfile:
    for f in read_files:
        with open(f, "rb") as infile:
            outfile.write(infile.read())
outfile.close()

#Strip: headers, varaints not accepted by cBioPortal, and genes know to give a large number of false positives - i.e. olfactory receptors and large genes. 
bad_orfs = ['#version', 'Intron', 'Silent', "3'UTR", "5'UTR", "3'Flank", 'Splice_Region', 'CallFreq', 'Unknown', 'unprocessed_pseudogene', 'LowEVS', 'MUC16', 'rs753994746', 'rs868981506', 'rs780538719', 'rs768409846', 'rs199941915', 'rs866674389', 'rs528635259', 'rs1384446739', 'rs201899861', 'rs200253954', 'rs113202486', 'rs1993829',
            'HLA-DQB2', 'COL3A1', 'rs751112872', 'rs1472236985', 'rs1366091859', 'rs561234742', 'rs532845932', 'rs11322783', 'rs113322110', 'rs11322783', 'rs4782272', 'rs878943546', 'ZMYM3', 'rs1375638824', 'rs1554844974', 'rs1553961403', 'rs1555889552', 'rs748260316', 'rs1422380205', 'rs745817378', 'rs763997493', 'rs767064362', 'rs1557586047', 'rs756366019', 'rs756349688', 'rs755706336', 'rs369823368', 'rs745614033', 'HLA-A', 'HLA-B', 'HLA-C', 'IGFN1', 'NBPF9', 'NXNL1', 'NPIPB6', 'AC008686.1', 'KIAA0355', 'RAPGEF2', 'ZXDA', 'GJC2', 'rs868981506', 'rs575085', 'OR10A2', 'OR10A3', 'OR10A4', 'OR10A5', 'OR10A6', 'OR10A7', 'OR10AA1P', 'OR10AB1P', 'OR10AC1', 'OR10AD1', 'OR10AE1P',
            'OR10AE3P', 'OR10AF1P', 'OR10AG1', 'OR10AH1P', 'OR10AK1P', 'OR10B1P', 'OR10C1', 'OR10D1P', 'OR10D3', 'OR10D4P', 'OR10D5P', 'OR10G1P', 'OR10G2', 'OR10G3', 'OR10G4', 'OR10G5P', 'OR10G6', 'rs765819981', 'rs200684350', 'rs149815285', 'rs781577243', 'rs779610049', 'rs12627379', 'rs778867790', 'rs751186374', 'rs200075071', 'rs1240508430', 'rs773963042', 'rs765451626',
            'OR10G7', 'OR10G8', 'OR10G9', 'OR10H1', 'OR10H2', 'OR10H3', 'OR10H4', 'OR10H5', 'OR10J1', 'OR10J2P', 'OR10J3', 'OR10J4', 'OR10J5', 'OR10J6P', 'OR10J7P', 'OR10J8P', 'OR10J9P', 'OR10K1', 'GOLGA6L20', 'ZNF575', 'BIRC6', 'ANKRD36', 'LRP1B', 'SCN7A', 'FSIP2', 'BAGE2', 'TOPAZ1', 'SLMAP', 'FRG2C', 'KIAA1109', 'DNAH8', 'RBM33', 'PKHD1L1', 'VPS13A', 'VPS13D',
            'OR10K2', 'OR10N1P', 'OR10P1', 'OR10Q1', 'OR10Q2P', 'OR10R1P', 'OR10R2', 'OR10R3P', 'OR10S1', 'OR10T1P', 'OR10T2', 'OR10U1P', 'OR10V1', 'OR10V2P', 'OR10V3P', 'OR10V7P', 'OR10W1', 'OR10X1',
            'OR10Y1P', 'OR10Z1', 'OR11A1', 'OR11G1P', 'OR11G2', 'OR11H1', 'OR11H12', 'OR11H13P', 'OR11H2', 'OR11H3P', 'OR11H4', 'OR11H5P', 'OR11H6', 'OR11H7', 'OR11I1P', 'OR11J1P', 'OR11J2P',
            'OR11J5P', 'OR11K1P', 'OR11K2P', 'OR11L1', 'OR11M1P', 'OR11N1P', 'OR11P1P', 'OR11Q1P', 'OR12D1', 'OR12D2', 'OR12D3', 'OR13A1', 'OR13C1P', 'OR13C2', 'OR13C3', 'OR13C4', 'OR13C5', 'OR13C6P',
            'OR13C7', 'OR13C8', 'OR13C9', 'OR13D1', 'OR13D2P', 'OR13D3P', 'OR13E1P', 'OR13F1', 'OR13G1', 'OR13H1', 'OR13I1P', 'OR13J1', 'OR13K1P', 'OR13Z1P', 'OR13Z2P', 'OR13Z3P', 'OR14A16',
            'OR14A2', 'OR14C36', 'OR14I1', 'OR14J1', 'OR14K1', 'OR14L1P', 'OR1A1', 'OR1A2', 'OR1AA1P', 'OR1AB1P', 'OR1AC1P', 'OR1B1', 'OR1C1', 'OR1D2', 'OR1D3P', 'OR1D4', 'OR1D5', 'OR1E1', 'OR1E2',
            'OR1E3', 'OR1F1', 'OR1F12', 'OR1F2P', 'OR1G1', 'OR1H1P', 'OR1I1', 'OR1J1', 'OR1J2', 'OR1J4', 'OR1K1', 'OR1L1', 'OR1L3', 'OR1L4', 'OR1L6', 'OR1L8', 'OR1M1', 'OR1M4P', 'OR1N1', 'OR1N2',
            'OR1P1', 'OR1Q1', 'OR1R1P', 'OR1S1', 'OR1S2', 'OR1X1P', 'OR1X5P', 'OR2A1', 'OR2A12', 'OR2A13P', 'OR2A14', 'OR2A15P', 'OR2A2', 'OR2A20P', 'OR2A25', 'OR2A3P', 'OR2A4', 'OR2A41P', 'OR2A42',
            'OR2A5', 'OR2A7', 'OR2A9P', 'OR2AD1P', 'OR2AE1', 'OR2AF1P', 'OR2AG1', 'OR2AG2', 'OR2AH1P', 'OR2AI1P', 'OR2AJ1', 'OR2AK2', 'OR2AL1P', 'OR2AM1P', 'OR2AO1P', 'OR2AP1', 'OR2AQ1P', 'OR2AS1P',
            'OR2AS2P', 'OR2AT1P', 'OR2AT2P', 'OR2AT4', 'OR2B11', 'OR2B2', 'OR2B3', 'OR2B4P', 'OR2B6', 'OR2B7P', 'OR2B8P', 'OR2BH1P', 'OR2C1', 'OR2C3', 'OR2D2', 'OR2D3', 'OR2E1P', 'OR2F1', 'OR2F2',
            'OR2G1P', 'OR2G2', 'OR2G3', 'OR2G6', 'OR2H1', 'OR2H2', 'OR2H4P', 'OR2H5P', 'OR2I1P', 'OR2J1', 'OR2J2', 'OR2J3', 'OR2J4P', 'OR2K2', 'OR2L13', 'OR2L1P', 'OR2L2', 'OR2L3', 'OR2L5', 'OR2L6P',
            'OR2L8', 'OR2L9P', 'OR2M1P', 'OR2M2', 'OR2M3', 'OR2M4', 'OR2M5', 'OR2M7', 'OR2N1P', 'OR2P1P', 'OR2Q1P', 'OR2R1P', 'OR2S1P', 'OR2S2', 'OR2T1', 'OR2T10', 'OR2T11', 'OR2T12', 'OR2T2',
            'OR2T27', 'OR2T29', 'OR2T3', 'OR2T32P', 'OR2T33', 'OR2T34', 'OR2T35', 'OR2T4', 'OR2T5', 'OR2T6', 'OR2T7', 'OR2T8', 'OR2U1P', 'OR2U2P', 'OR2V1', 'OR2V2', 'OR2W1', 'OR2W2P', 'OR2W3',
            'OR2W4P', 'OR2W5', 'OR2W6P', 'OR2X1P', 'OR2Y1', 'OR2Z1', 'OR3A1', 'OR3A2', 'OR3A3', 'OR3A4P', 'OR3B1P', 'OR3D1P', 'OR4A10P', 'OR4A11P', 'OR4A12P', 'OR4A13P', 'OR4A14P', 'OR4A15', 'OR4A16',
            'OR4A17P', 'OR4A18P', 'OR4A19P', 'OR4A1P', 'OR4A21P', 'OR4A2P', 'OR4A3P', 'OR4A40P', 'OR4A41P', 'OR4A42P', 'OR4A43P', 'OR4A44P', 'OR4A45P', 'OR4A46P', 'OR4A47', 'OR4A48P', 'OR4A49P',
            'OR4A4P', 'OR4A5', 'OR4A50P', 'OR4A6P', 'OR4A7P', 'OR4A8', 'OR4A9P', 'OR4B1', 'OR4B2P', 'OR4C10P', 'OR4C11', 'OR4C12', 'OR4C13', 'OR4C14P', 'OR4C15', 'OR4C16', 'OR4C1P', 'OR4C2P', 'OR4C3',
            'OR4C45', 'OR4C46', 'OR4C48P', 'OR4C49P', 'OR4C4P', 'OR4C5', 'OR4C50P', 'OR4C6', 'OR4C7P', 'OR4C9P', 'OR4D1', 'OR4D10', 'OR4D11', 'OR4D12P', 'OR4D2', 'OR4D5', 'OR4D6', 'OR4D7P', 'OR4D8P',
            'OR4D9', 'OR4E1', 'OR4E2', 'OR4F13P', 'OR4F14P', 'OR4F15', 'OR4F16', 'OR4F17', 'OR4F1P', 'OR4F21', 'OR4F28P', 'OR4F29', 'OR4F2P', 'OR4F3', 'OR4F4', 'OR4F5', 'OR4F6', 'OR4F7P', 'OR4F8P',
            'OR4G11P', 'OR4G1P', 'OR4G2P', 'OR4G3P', 'OR4G4P', 'OR4G6P', 'OR4H12P', 'OR4H6P', 'OR4K1', 'OR4K11P', 'OR4K12P', 'OR4K13', 'OR4K14', 'OR4K15', 'OR4K16P', 'OR4K17', 'OR4K2', 'OR4K3',
            'OR4K4P', 'OR4K5', 'OR4K6P', 'OR4K7P', 'OR4K8P', 'OR4L1', 'OR4M1', 'OR4M2', 'OR4N1P', 'OR4N2', 'OR4N3P', 'OR4N4', 'OR4N5', 'OR4P1P', 'OR4P4', 'OR4Q1P', 'OR4Q2', 'OR4Q3', 'OR4R1P',
            'OR4R2P', 'OR4R3P', 'OR4S1', 'OR4S2', 'OR4T1P', 'OR4U1P', 'OR4V1P', 'OR4W1P', 'OR4X1', 'OR4X2', 'OR4X7P', 'OR51A10P', 'OR51A1P', 'OR51A2', 'OR51A3P', 'OR51A4', 'OR51A5P', 'OR51A6P',
            'OR51A7', 'OR51A8P', 'OR51A9P', 'OR51AB1P', 'OR51B2', 'OR51B3P', 'OR51B4', 'OR51B5', 'OR51B6', 'OR51B8P', 'OR51C1P', 'OR51C4P', 'OR51D1', 'OR51E1', 'OR51E2', 'OR51F1', 'OR51F2',
            'OR51F3P', 'OR51F4P', 'OR51F5P', 'OR51G1', 'OR51G2', 'OR51H1', 'OR51H2P', 'OR51I1', 'OR51I2', 'OR51J1', 'OR51K1P', 'OR51L1', 'OR51M1', 'OR51N1P', 'OR51P1P', 'OR51Q1', 'OR51R1P', 'OR51S1',
            'OR51T1', 'OR51V1', 'OR52A1', 'OR52A4P', 'OR52A5', 'OR52B1P', 'OR52B2', 'OR52B3P', 'OR52B4', 'OR52B5P', 'OR52B6', 'OR52D1', 'OR52E1', 'OR52E2', 'OR52E3P', 'OR52E4', 'OR52E5', 'OR52E6',
            'OR52E7P', 'OR52E8', 'OR52H1', 'OR52H2P', 'OR52I1', 'OR52I2', 'OR52J1P', 'OR52J2P', 'OR52J3', 'OR52K1', 'OR52K2', 'OR52K3P', 'OR52L1', 'OR52L2P', 'OR52M1', 'OR52M2P', 'OR52N1', 'OR52N2',
            'OR52N3P', 'OR52N4', 'OR52N5', 'OR52P1P', 'OR52P2P', 'OR52Q1P', 'OR52R1', 'OR52S1P', 'OR52T1P', 'OR52U1P', 'OR52V1P', 'OR52W1', 'OR52X1P', 'OR52Y1P', 'OR52Z1', 'OR55B1P', 'OR56A1',
            'OR56A3', 'OR56A4', 'OR56A5', 'OR56A7P', 'OR56B1', 'OR56B2P', 'OR56B3P', 'OR56B4', 'OR5A1', 'OR5A2', 'OR5AC1', 'OR5AC2', 'OR5AC4P', 'OR5AH1P', 'OR5AK1P', 'OR5AK2', 'OR5AK3P', 'OR5AK4P',
            'OR5AL1', 'OR5AL2P', 'OR5AM1P', 'OR5AN1', 'OR5AN2P', 'OR5AO1P', 'OR5AP1P', 'OR5AP2', 'OR5AQ1P', 'OR5AR1', 'OR5AS1', 'OR5AU1', 'OR5AW1P', 'OR5AZ1P', 'OR5B10P', 'OR5B12', 'OR5B15P',
            'OR5B17', 'OR5B19P', 'OR5B1P', 'OR5B2', 'OR5B21', 'OR5B3', 'OR5BA1P', 'OR5BB1P', 'OR5BC1P', 'OR5BD1P', 'OR5BE1P', 'OR5BH1P', 'OR5BJ1P', 'OR5BK1P', 'OR5BL1P', 'OR5BM1P', 'OR5BN1P',
            'OR5BN2P', 'OR5BP1P', 'OR5BQ1P', 'OR5BR1P', 'OR5BS1P', 'OR5BT1P', 'OR5C1', 'OR5D13', 'OR5D14', 'OR5D15P', 'OR5D16', 'OR5D17P', 'OR5D18', 'OR5D2P', 'OR5D3P', 'OR5E1P', 'OR5F1', 'OR5F2P',
            'OR5G1P', 'OR5G3', 'OR5G4P', 'OR5G5P', 'OR5H1', 'OR5H14', 'OR5H15', 'OR5H2', 'OR5H3P', 'OR5H4P', 'OR5H5P', 'OR5H6', 'OR5H7P', 'OR5H8', 'OR5I1', 'OR5J1P', 'OR5J2', 'OR5J7P', 'OR5K1',
            'OR5K2', 'OR5K3', 'OR5K4', 'OR5L1', 'OR5L2', 'OR5M1', 'OR5M10', 'OR5M11', 'OR5M12P', 'OR5M13P', 'OR5M14P', 'OR5M2P', 'OR5M3', 'OR5M4P', 'OR5M5P', 'OR5M6P', 'OR5M7P', 'OR5M8', 'OR5M9',
            'OR5P1P', 'OR5P2', 'OR5P3', 'OR5P4P', 'OR5R1', 'OR5S1P', 'OR5T1', 'OR5T2', 'OR5T3', 'OR5V1', 'OR5W1P', 'OR5W2', 'OR6A2', 'OR6B1', 'OR6B2', 'OR6B3', 'OR6C1', 'OR6C2', 'OR6C3', 'OR6C4',
            'OR6C5P', 'OR6C6', 'OR6C64P', 'OR6C65', 'OR6C66P', 'OR6C68', 'OR6C69P', 'OR6C70', 'OR6C71P', 'OR6C72P', 'OR6C73P', 'OR6C74', 'OR6C75', 'OR6C76', 'OR6C7P', 'OR6D1P', 'OR6E1P', 'OR6F1',
            'OR6J1', 'OR6K1P', 'OR6K2', 'OR6K3', 'OR6K4P', 'OR6K5P', 'OR6K6', 'OR6L1P', 'OR6L2P', 'OR6M1', 'OR6M2P', 'OR6M3P', 'OR6N1', 'OR6N2', 'OR6P1', 'OR6Q1', 'OR6R1P', 'OR6R2P', 'OR6S1',
            'OR6T1', 'OR6U2P', 'OR6V1', 'OR6W1P', 'OR6X1', 'OR6Y1', 'OR7A10', 'OR7A11P', 'OR7A15P', 'OR7A17', 'OR7A18P', 'OR7A19P', 'OR7A1P', 'OR7A2P', 'OR7A3P', 'OR7A5', 'OR7A8P', 'OR7C1', 'OR7C2',
            'OR7D11P', 'OR7D1P', 'OR7D2', 'OR7D4', 'OR7E100P', 'OR7E101P', 'OR7E102P', 'OR7E104P', 'OR7E105P', 'OR7E106P', 'OR7E108P', 'OR7E109P', 'OR7E10P', 'OR7E110P', 'OR7E111P', 'OR7E115P',
            'OR7E116P', 'OR7E117P', 'OR7E11P', 'OR7E121P', 'OR7E122P', 'OR7E125P', 'OR7E126P', 'OR7E128P', 'OR7E129P', 'OR7E12P', 'OR7E130P', 'OR7E136P', 'OR7E13P', 'OR7E140P', 'OR7E145P', 'OR7E148P',
            'OR7E149P', 'OR7E14P', 'OR7E154P', 'OR7E155P', 'OR7E156P', 'OR7E157P', 'OR7E158P', 'OR7E159P', 'OR7E15P', 'OR7E160P', 'OR7E161P', 'OR7E162P', 'OR7E163P', 'OR7E16P', 'OR7E18P', 'OR7E19P',
            'OR7E1P', 'OR7E21P', 'OR7E22P', 'OR7E23P', 'OR7E24', 'OR7E25P', 'OR7E26P', 'OR7E28P', 'OR7E29P', 'OR7E2P', 'OR7E31P', 'OR7E33P', 'OR7E35P', 'OR7E36P', 'OR7E37P', 'OR7E38P', 'OR7E39P',
            'OR7E41P', 'OR7E43P', 'OR7E46P', 'OR7E47P', 'OR7E4P', 'OR7E53P', 'OR7E55P', 'OR7E59P', 'OR7E5P', 'OR7E62P', 'OR7E66P', 'OR7E7P', 'OR7E83P', 'OR7E84P', 'OR7E85P', 'OR7E86P', 'OR7E87P',
            'OR7E89P', 'OR7E8P', 'OR7E90P', 'OR7E91P', 'OR7E93P', 'OR7E94P', 'OR7E96P', 'OR7E97P', 'OR7E99P', 'OR7G1', 'OR7G15P', 'OR7G2', 'OR7G3', 'OR7H1P', 'OR7H2P', 'OR7K1P', 'OR7L1P', 'OR7M1P',
            'OR8A1', 'OR8A2P', 'OR8A3P', 'OR8B10P', 'OR8B12', 'OR8B1P', 'OR8B2', 'OR8B3', 'OR8B4', 'OR8B5P', 'OR8B6P', 'OR8B7P', 'OR8B8', 'OR8B9P', 'OR8C1P', 'OR8D1', 'OR8D2', 'OR8D4', 'OR8F1P',
            'OR8G1', 'OR8G2P', 'OR8G3P', 'OR8G5', 'OR8G7P', 'OR8H1', 'OR8H2', 'OR8H3', 'OR8I1P', 'OR8I2', 'OR8I4P', 'OR8J1', 'OR8J2', 'OR8J3', 'OR8K1', 'OR8K2P', 'OR8K3', 'OR8K4P', 'OR8K5', 'OR8L1P',
            'OR8Q1P', 'OR8R1P', 'OR8S1', 'OR8S21P', 'OR8T1P', 'OR8U1', 'OR8U8', 'OR8U9', 'OR8V1P', 'OR8X1P', 'OR9A1P', 'OR9A2', 'OR9A3P', 'OR9A4', 'OR9G1', 'OR9G2P', 'OR9G3P', 'OR9G4', 'OR9G9',
            'OR9H1P', 'OR9I1', 'OR9I2P', 'OR9I3P', 'OR9K1P', 'OR9K2', 'OR9L1P', 'OR9M1P', 'OR9N1P', 'OR9P1P', 'OR9Q1', 'OR9Q2', 'OR9R1P', 'OR9S24P', 'rs202003805', 'GSTT4', 'rs3825942',
            'MUC4', 'MUC6', 'MUC12', 'HLA-DRB5', 'MUC17', 'IGHV4-31', 'LRRC37A3', 'MUC3A', 'PRSS1', 'MUC5B', 'TRBV3-1', 'GOLGA6L10', 'TAS2R43','IGHV2-5', 'GOLGA6L2', 'IGHV4-59', 'GOLGA6L2', 'SLC25A5', 'RP1L1', 'TRBV4-2', 'TAS2R19', 'AC118281.1', 'IGLV5-45', 'TAS2R31', 'HRNR', 'MADCAM1', 'PABPC3', 'NPIPB15', 'EVL', 'AHNAK2', 'PRAMEF33', 'TAS2R30', 'FLG', 'MAGEC1', 'HLA-DQA2', 'GXYLT1', 'IGHV4-61', 'OGFR', 'MUC21', 'PLIN4',
            'rs73126218' 'rs750338758', 'rs533172496', ] 
#Strip all ACMG, CPIC, and Foundation ORFs:
good_orfs = ['common_variant', 'ABL1', 'ABL2', 'ACTA2', 'ACTC1', 'ACVR1B', 'AIP', 'AKT1', 'AKT2', 'AKT3', 'ALK', 'AMER1', 'ANKRD26', 'APC', 'APOB', 'AR', 'ARAF', 'ARFRP1', 'ARID1A', 'ARID1B', 'ARID2', 'ASXL1', 
             'ATM', 'ATP7B', 'ATR', 'ATRX', 'AURKA', 'AURKB', 'AXIN1', 'AXIN2', 'AXL', 'BAP1', 'BARD1', 'BCL2', 'BCL2L1', 'BCL2L2', 'BCL6', 'BCOR', 'BCORL1', 'BLM', 'BMPR1A', 'BRAF', 'BRCA1', 
             'BRCA2', 'BRD4', 'BRIP1', 'BTG1', 'BTK', 'BUB1B', 'CACA1C', 'CACNA1C', 'CACNA1S', 'CARD11', 'CASR', 'CBFB', 'CBL', 'CCND1', 'CCND2', 'CCND3', 'CCNE1', 'CD274', 'CD79A', 'CD79B', 
             'CDC73', 'CDH1', 'CDK12', 'CDK4', 'CDK6', 'CDK8', 'CDKN1A', 'CDKN1B', 'CDKN1C', 'CDKN2A', 'CDKN2B', 'CDKN2C', 'CEBPA', 'CEP57', 'CFTR', 'CHD2', 'CHD4', 'CHEK1', 'CHEK2', 'CIC', 
             'COL3A1', 'CPA1', 'CREBBP', 'CRKL', 'CRLF2', 'CSF1R', 'CTC1', 'CTCF', 'CTNNA1', 'CTNNB1', 'CTRC', 'CUL3', 'CYLD', 'CYP2C19', 'CYP2C9', 'CYP2D6', 'CYP3A5', 'CYP4F2', 'DAXX', 
             'DDR2', 'DDX41', 'DES', 'DICER1', 'DIS3L2', 'DKC1', 'DNAJC21', 'DNMT3A', 'DOT1L',  'DSC2', 'DSG2', 'DSP', 'EFL1', 'EGFR', 'EMSY', 'EP300', 'EPCAM', 'EPHA3', 'EPHA5', 
             'EPHA7', 'EPHB1', 'ERBB2', 'ERBB3', 'ERBB4', 'ERCC4', 'ERG', 'ERRFI1', 'ESR1', 'ETV6', 'EZH2', 'FAM175A', 'FAM46C', 'FANCA', 'FANCB', 'FANCC', 'FANCD2', 'FANCE', 'FANCF', 
             'FANCG', 'FANCI', 'FANCL', 'FANCM', 'FAS', 'FAT1', 'FBN1', 'FBXW7', 'FGF10', 'FGF14', 'FGF19', 'FGF23', 'FGF3', 'FGF4', 'FGF6', 'FGFR1', 'FGFR2', 'FGFR3', 'FGFR4', 'FH', 'FLCN', 
             'FLT1', 'FLT3', 'FLT4', 'FOXL2', 'FOXP1', 'FRS2', 'FUBP1', 'G6PD', 'GABRA6', 'GALNT12', 'GATA1', 'GATA2', 'GATA3', 'GATA4', 'GATA6', 'GID4', 'GLA', 'GLI1', 'GNA11', 'GNA13', 
             'GNAQ', 'GNAS', 'GPC3', 'GPR124', 'GREM1', 'GRIN2A', 'GRM3', 'GSK3B', 'H3F3A', 'HGF', 'HNF1A', 'HOXB13', 'HRAS', 'HSD3B1', 'HSP90AA1', 'IDH1', 'IDH2', 'IFNL3', 'IGF1R', 'IGF2', 
             'IKBKE', 'IKZF1', 'IL7R', 'INHBA', 'INPP4B', 'IRF2', 'IRF4', 'IRS2', 'JAK1', 'JAK2', 'JAK3', 'JUN', 'KAT6A', 'KCNH2', 'KCNQ1', 'KDM5A', 'KDM5C', 'KDM6A', 'KDR', 'KEAP1', 'KEL', 
             'KIT', 'KLHL6', 'KMT2A', 'KMT2D', 'KRAS', 'LAMP2', 'LDLR', 'LMNA', 'LMO1', 'LRP1B', 'LYN', 'LZTR1', 'MAD2L2', 'MAGI2', 'MAP2K1', 'MAP2K2', 'MAP2K4', 'MAP3K1', 'MAX', 
             'MC1R', 'MCL1', 'MDM2', 'MDM4', 'MED12', 'MEF2B', 'MEN1', 'MET', 'MITF', 'MLH1', 'MPL', 'MRE11A', 'MSH2', 'MSH3', 'MSH6', 'MTOR', 'MUTYH', 'MYBPC3', 'MYC', 'MYCL', 'MYCL1', 'MYCN', 
             'MYD88', 'MYH11', 'MYH7', 'MYL2', 'MYL3', 'NBN', 'NDE1', 'NF1', 'NF2', 'NFE2L2', 'NFKBIA', 'NKX2-1', 'NOTCH1', 'NOTCH2', 'NOTCH3', 'NPM1', 'NRAS', 'NSD1', 'NTHL1', 'NTRK1', 'NTRK2', 
             'NTRK3', 'NUP93', 'OTC', 'PAK3', 'PALB2', 'PARK2', 'PAX5', 'PBRM1', 'PCSK9', 'PDCD1LG2', 'PDGFRA', 'PDGFRB', 'PDK1', 'PHOX2B', 'PIK3C2B', 'PIK3CA', 'PIK3CB', 'PIK3CG', 'PIK3R1', 
             'PIK3R2', 'PKP2', 'PLCG2', 'PMS2', 'POLD1', 'POLE', 'POT1', 'PPP2R1A', 'PRDM1', 'PREX2', 'PRKAG2', 'PRKAR1A', 'PRKCI', 'PRKDC', 'PRSS1', 'PRSS8', 'PTCH1', 'PTEN', 'PTPN11', 'QKI', 
             'RAC1', 'RAD50', 'RAD51', 'RAD51C', 'RAD51D', 'RAF1', 'RANBP2', 'RARA', 'RB1', 'RBM10', 'RECQL4', 'RET', 'RFWD3', 'RICTOR', 'RINT1 ', 'RNF43', 'ROS1', 'RPL11', 'RPL35A', 'RPL5', 
             'RPS10', 'RPS17', 'RPS19', 'RPS20', 'RPS24', 'RPS26', 'RPTOR', 'RTEL1', 'RUNX1', 'RUNX1T1', 'RYR1', 'RYR2', 'SBDS', 'SCN5A', 'SDHA', 'SDHAF2', 'SDHB', 'SDHC', 'SDHD', 'SETD2', 
             'SF3B1', 'SLCO1B1', 'SLIT2', 'SLX4', 'SMAD2', 'SMAD3', 'SMAD4', 'SMARCA4', 'SMARCB1', 'SMARCE1', 'SMO', 'SNCAIP', 'SOCS1', 'SOX10', 'SOX2', 'SOX9', 'SPEN', 'SPINK1', 'SPOP', 
             'SPTA1', 'SRC', 'SRP54', 'SRP72', 'STAG2', 'STAT3', 'STAT4', 'STK11', 'SUFU', 'SYK', 'TAF1', 'TBX3', 'TERC', 'TERT', 'TET2', 'TG', 'TGFB2', 'TGFB3', 'TGFBR1', 'TGFBR2', 'TINF2', 
             'TMEM127', 'TMEM43', 'TNFAIP3', 'TNFRSF14', 'TNN13', 'TNNI3', 'TNNT2', 'TOP1', 'TOP2A', 'TP53', 'TPM1', 'TPMT', 'TSC1', 'TSC2', 'TSHR', 'U2AF1', 'UBE2T', 'UGT1A1', 'VEGFA', 'VHL', 
             'VKORC1', 'WISP3', 'WRN', 'WT1', 'XPO1', 'XRCC2', 'ZBTB2', 'ZNF217', 'ZNF703']
#Grab all ACMG, CPIC, and Foundation ORFs used for dataframe1
acmg_orfs = ['ABL1', 'ABL2', 'ACTA2', 'ACTC1', 'ACVR1B', 'AIP', 'AKT1', 'AKT2', 'AKT3', 'ALK', 'AMER1', 'ANKRD26', 'APC', 'APOB', 'AR', 'ARAF', 'ARFRP1', 'ARID1A', 'ARID1B', 'ARID2', 'ASXL1', 
             'ATM', 'ATP7B', 'ATR', 'ATRX', 'AURKA', 'AURKB', 'AXIN1', 'AXIN2', 'AXL', 'BAP1', 'BARD1', 'BCL2', 'BCL2L1', 'BCL2L2', 'BCL6', 'BCOR', 'BCORL1', 'BLM', 'BMPR1A', 'BRAF', 'BRCA1', 
             'BRCA2', 'BRD4', 'BRIP1', 'BTG1', 'BTK', 'BUB1B', 'CACA1C', 'CACNA1C', 'CACNA1S', 'CARD11', 'CASR', 'CBFB', 'CBL', 'CCND1', 'CCND2', 'CCND3', 'CCNE1', 'CD274', 'CD79A', 'CD79B', 
             'CDC73', 'CDH1', 'CDK12', 'CDK4', 'CDK6', 'CDK8', 'CDKN1A', 'CDKN1B', 'CDKN1C', 'CDKN2A', 'CDKN2B', 'CDKN2C', 'CEBPA', 'CEP57', 'CFTR', 'CHD2', 'CHD4', 'CHEK1', 'CHEK2', 'CIC', 
             'COL3A1', 'CPA1', 'CREBBP', 'CRKL', 'CRLF2', 'CSF1R', 'CTC1', 'CTCF', 'CTNNA1', 'CTNNB1', 'CTRC', 'CUL3', 'CYLD', 'CYP2C19', 'CYP2C9', 'CYP2D6', 'CYP3A5', 'CYP4F2', 'DAXX', 
             'DDR2', 'DDX41', 'DES', 'DICER1', 'DIS3L2', 'DKC1', 'DNAJC21', 'DNMT3A', 'DOT1L',  'DSC2', 'DSG2', 'DSP', 'EFL1', 'EGFR', 'EMSY', 'EP300', 'EPCAM', 'EPHA3', 'EPHA5', 
             'EPHA7', 'EPHB1', 'ERBB2', 'ERBB3', 'ERBB4', 'ERCC4', 'ERG', 'ERRFI1', 'ESR1', 'ETV6', 'EZH2', 'FAM175A', 'FAM46C', 'FANCA', 'FANCB', 'FANCC', 'FANCD2', 'FANCE', 'FANCF', 
             'FANCG', 'FANCI', 'FANCL', 'FANCM', 'FAS', 'FAT1', 'FBN1', 'FBXW7', 'FGF10', 'FGF14', 'FGF19', 'FGF23', 'FGF3', 'FGF4', 'FGF6', 'FGFR1', 'FGFR2', 'FGFR3', 'FGFR4', 'FH', 'FLCN', 
             'FLT1', 'FLT3', 'FLT4', 'FOXL2', 'FOXP1', 'FRS2', 'FUBP1', 'G6PD', 'GABRA6', 'GALNT12', 'GATA1', 'GATA2', 'GATA3', 'GATA4', 'GATA6', 'GID4', 'GLA', 'GLI1', 'GNA11', 'GNA13', 
             'GNAQ', 'GNAS', 'GPC3', 'GPR124', 'GREM1', 'GRIN2A', 'GRM3', 'GSK3B', 'H3F3A', 'HGF', 'HNF1A', 'HOXB13', 'HRAS', 'HSD3B1', 'HSP90AA1', 'IDH1', 'IDH2', 'IFNL3', 'IGF1R', 'IGF2', 
             'IKBKE', 'IKZF1', 'IL7R', 'INHBA', 'INPP4B', 'IRF2', 'IRF4', 'IRS2', 'JAK1', 'JAK2', 'JAK3', 'JUN', 'KAT6A', 'KCNH2', 'KCNQ1', 'KDM5A', 'KDM5C', 'KDM6A', 'KDR', 'KEAP1', 'KEL', 
             'KIT', 'KLHL6', 'KMT2A', 'KMT2D', 'KRAS', 'LAMP2', 'LDLR', 'LMNA', 'LMO1', 'LRP1B', 'LYN', 'LZTR1', 'MAD2L2', 'MAGI2', 'MAP2K1', 'MAP2K2', 'MAP2K4', 'MAP3K1', 'MAX', 
             'MC1R', 'MCL1', 'MDM2', 'MDM4', 'MED12', 'MEF2B', 'MEN1', 'MET', 'MITF', 'MLH1', 'MPL', 'MRE11A', 'MSH2', 'MSH3', 'MSH6', 'MTOR', 'MUTYH', 'MYBPC3', 'MYC', 'MYCL', 'MYCL1', 'MYCN', 
             'MYD88', 'MYH11', 'MYH7', 'MYL2', 'MYL3', 'NBN', 'NDE1', 'NF1', 'NF2', 'NFE2L2', 'NFKBIA', 'NKX2-1', 'NOTCH1', 'NOTCH2', 'NOTCH3', 'NPM1', 'NRAS', 'NSD1', 'NTHL1', 'NTRK1', 'NTRK2', 
             'NTRK3', 'NUP93', 'OTC', 'PAK3', 'PALB2', 'PARK2', 'PAX5', 'PBRM1', 'PCSK9', 'PDCD1LG2', 'PDGFRA', 'PDGFRB', 'PDK1', 'PHOX2B', 'PIK3C2B', 'PIK3CA', 'PIK3CB', 'PIK3CG', 'PIK3R1', 
             'PIK3R2', 'PKP2', 'PLCG2', 'PMS2', 'POLD1', 'POLE', 'POT1', 'PPP2R1A', 'PRDM1', 'PREX2', 'PRKAG2', 'PRKAR1A', 'PRKCI', 'PRKDC', 'PRSS1', 'PRSS8', 'PTCH1', 'PTEN', 'PTPN11', 'QKI', 
             'RAC1', 'RAD50', 'RAD51', 'RAD51C', 'RAD51D', 'RAF1', 'RANBP2', 'RARA', 'RB1', 'RBM10', 'RECQL4', 'RET', 'RFWD3', 'RICTOR', 'RINT1 ', 'RNF43', 'ROS1', 'RPL11', 'RPL35A', 'RPL5', 
             'RPS10', 'RPS17', 'RPS19', 'RPS20', 'RPS24', 'RPS26', 'RPTOR', 'RTEL1', 'RUNX1', 'RUNX1T1', 'RYR1', 'RYR2', 'SBDS', 'SCN5A', 'SDHA', 'SDHAF2', 'SDHB', 'SDHC', 'SDHD', 'SETD2', 
             'SF3B1', 'SLCO1B1', 'SLIT2', 'SLX4', 'SMAD2', 'SMAD3', 'SMAD4', 'SMARCA4', 'SMARCB1', 'SMARCE1', 'SMO', 'SNCAIP', 'SOCS1', 'SOX10', 'SOX2', 'SOX9', 'SPEN', 'SPINK1', 'SPOP', 
             'SPTA1', 'SRC', 'SRP54', 'SRP72', 'STAG2', 'STAT3', 'STAT4', 'STK11', 'SUFU', 'SYK', 'TAF1', 'TBX3', 'TERC', 'TERT', 'TET2', 'TG', 'TGFB2', 'TGFB3', 'TGFBR1', 'TGFBR2', 'TINF2', 
             'TMEM127', 'TMEM43', 'TNFAIP3', 'TNFRSF14', 'TNN13', 'TNNI3', 'TNNT2', 'TOP1', 'TOP2A', 'TP53', 'TPM1', 'TPMT', 'TSC1', 'TSC2', 'TSHR', 'U2AF1', 'UBE2T', 'UGT1A1', 'VEGFA', 'VHL', 
             'VKORC1', 'WISP3', 'WRN', 'WT1', 'XPO1', 'XRCC2', 'ZBTB2', 'ZNF217', 'ZNF703']

with open("AvatarSomaticACMG_Master.maf") as infile, open ("AvatarSomaticMasterFilteredTemp.maf", "w") as outfile:
#strip headers, flags, and olfactory receptors from list above 
        for line in infile:
                if not any(bad_orf in line for bad_orf in bad_orfs):
                    outfile.write(line)
                    
#Strip ALL acmg, cpic and foundation orfs
with open("AvatarSomaticMasterFilteredTemp.maf") as infile2, open ("AvatarSomaticNoACMG_MasterFiltered2Temp.maf", "w") as outfile2:
        for line in infile2:
                if not any(good_orf in line for good_orf in good_orfs):
                        outfile2.write(line)
                        
#Grab ALL ACMG, CPIC, and Foundation variants even if a "common_variant"              
columns = ['Hugo_Symbol', 'Entrez_Gene_Id', 'Center', 'NCBI_Build', 'Chromosome', 'Start_Position', 'End_Position', 'Strand', 'Variant_Classification', 'Variant_Type', 'Reference_Allele', 'Tumor_Seq_Allele1', 'Tumor_Seq_Allele2', 'dbSNP_RS', 'dbSNP_Val_Status', 'Tumor_Sample_Barcode', 'Matched_Norm_Sample_Barcode', 'Match_Norm_Seq_Allele1', 'Match_Norm_Seq_Allele2', 'Tumor_Validation_Allele1', 'Tumor_Validation_Allele2', 'Match_Norm_Validation_Allele1', 'Match_Norm_Validation_Allele2', 'Verification_Status', 'Validation_Status', 'Mutation_Status', 'Sequencing_Phase', 'Sequence_Source', 'Validation_Method', 'Score', 'BAM_File', 'Sequencer', 'Tumor_Sample_UUID', 'Matched_Norm_Sample_UUID', 'HGVSc', 'HGVSp', 'HGVSp_Short', 'Transcript_ID', 'Exon_Number', 't_depth', 't_ref_count', 't_alt_count', 'n_depth', 'n_ref_count', 'n_alt_count', 'all_effects', 'Allele', 'Gene', 'Feature', 'Feature_type', 'Consequence', 'cDNA_position', 'CDS_position', 'Protein_position', 'Amino_acids', 'Codons', 'Existing_variation', 'ALLELE_NUM', 'DISTANCE', 'STRAND_VEP', 'SYMBOL', 'SYMBOL_SOURCE', 'HGNC_ID', 'BIOTYPE', 'CANONICAL', 'CCDS', 'ENSP', 'SWISSPROT', 'TREMBL', 'UNIPARC', 'RefSeq', 'SIFT', 'PolyPhen', 'EXON', 'INTRON', 'DOMAINS', 'AF', 'AFR_AF', 'AMR_AF', 'ASN_AF', 'EAS_AF', 'EUR_AF', 'SAS_AF', 'AA_AF', 'EA_AF', 'CLIN_SIG', 'SOMATIC', 'PUBMED', 'MOTIF_NAME', 'MOTIF_POS', 'HIGH_INF_POS', 'MOTIF_SCORE_CHANGE', 'IMPACT', 'PICK', 'VARIANT_CLASS', 'TSL', 'HGVS_OFFSET', 'PHENO', 'MINIMISED', 'GENE_PHENO', 'FILTER', 'flanking_bps', 'vcf_id', 'vcf_qual', 'gnomAD_AF', 'gnomAD_AFR_AF', 'gnomAD_AMR_AF', 'gnomAD_ASJ_AF', 'gnomAD_EAS_AF', 'gnomAD_FIN_AF', 'gnomAD_NFE_AF', 'gnomAD_OTH_AF', 'gnomAD_SAS_AF', 'vcf_pos', 'T_DP', 'T_AF', 'N_DP', 'N_AF', 'BKZ']
df = (pd.read_table("AvatarSomaticMasterFilteredTemp.maf", usecols=columns))
#df = df.drop(list(df)[121:124], axis=1)
df = df[df['FILTER'].isin(['PASS'])]
df1 = df[df['Hugo_Symbol'].isin(['ABL1', 'ABL2', 'ACTA2', 'ACTC1', 'ACVR1B', 'AIP', 'AKT1', 'AKT2', 'AKT3', 'ALK', 'AMER1', 'ANKRD26', 'APC', 'APOB', 'AR', 'ARAF', 'ARFRP1', 'ARID1A', 'ARID1B', 'ARID2', 'ASXL1', 
                                 'ATM', 'ATP7B', 'ATR', 'ATRX', 'AURKA', 'AURKB', 'AXIN1', 'AXIN2', 'AXL', 'BAP1', 'BARD1', 'BCL2', 'BCL2L1', 'BCL2L2', 'BCL6', 'BCOR', 'BCORL1', 'BLM', 'BMPR1A', 'BRAF', 'BRCA1', 
                                 'BRCA2', 'BRD4', 'BRIP1', 'BTG1', 'BTK', 'BUB1B', 'CACA1C', 'CACNA1C', 'CACNA1S', 'CARD11', 'CASR', 'CBFB', 'CBL', 'CCND1', 'CCND2', 'CCND3', 'CCNE1', 'CD274', 'CD79A', 'CD79B', 
                                 'CDC73', 'CDH1', 'CDK12', 'CDK4', 'CDK6', 'CDK8', 'CDKN1A', 'CDKN1B', 'CDKN1C', 'CDKN2A', 'CDKN2B', 'CDKN2C', 'CEBPA', 'CEP57', 'CFTR', 'CHD2', 'CHD4', 'CHEK1', 'CHEK2', 'CIC', 
                                 'COL3A1', 'CPA1', 'CREBBP', 'CRKL', 'CRLF2', 'CSF1R', 'CTC1', 'CTCF', 'CTNNA1', 'CTNNB1', 'CTRC', 'CUL3', 'CYLD', 'CYP2C19', 'CYP2C9', 'CYP2D6', 'CYP3A5', 'CYP4F2', 'DAXX', 
                                 'DDR2', 'DDX41', 'DES', 'DICER1', 'DIS3L2', 'DKC1', 'DNAJC21', 'DNMT3A', 'DOT1L',  'DSC2', 'DSG2', 'DSP', 'EFL1', 'EGFR', 'EMSY', 'EP300', 'EPCAM', 'EPHA3', 'EPHA5', 
                                 'EPHA7', 'EPHB1', 'ERBB2', 'ERBB3', 'ERBB4', 'ERCC4', 'ERG', 'ERRFI1', 'ESR1', 'ETV6', 'EZH2', 'FAM175A', 'FAM46C', 'FANCA', 'FANCB', 'FANCC', 'FANCD2', 'FANCE', 'FANCF', 
                                 'FANCG', 'FANCI', 'FANCL', 'FANCM', 'FAS', 'FAT1', 'FBN1', 'FBXW7', 'FGF10', 'FGF14', 'FGF19', 'FGF23', 'FGF3', 'FGF4', 'FGF6', 'FGFR1', 'FGFR2', 'FGFR3', 'FGFR4', 'FH', 'FLCN', 
                                 'FLT1', 'FLT3', 'FLT4', 'FOXL2', 'FOXP1', 'FRS2', 'FUBP1', 'G6PD', 'GABRA6', 'GALNT12', 'GATA1', 'GATA2', 'GATA3', 'GATA4', 'GATA6', 'GID4', 'GLA', 'GLI1', 'GNA11', 'GNA13', 
                                 'GNAQ', 'GNAS', 'GPC3', 'GPR124', 'GREM1', 'GRIN2A', 'GRM3', 'GSK3B', 'H3F3A', 'HGF', 'HNF1A', 'HOXB13', 'HRAS', 'HSD3B1', 'HSP90AA1', 'IDH1', 'IDH2', 'IFNL3', 'IGF1R', 'IGF2', 
                                 'IKBKE', 'IKZF1', 'IL7R', 'INHBA', 'INPP4B', 'IRF2', 'IRF4', 'IRS2', 'JAK1', 'JAK2', 'JAK3', 'JUN', 'KAT6A', 'KCNH2', 'KCNQ1', 'KDM5A', 'KDM5C', 'KDM6A', 'KDR', 'KEAP1', 'KEL', 
                                 'KIT', 'KLHL6', 'KMT2A', 'KMT2D', 'KRAS', 'LAMP2', 'LDLR', 'LMNA', 'LMO1', 'LRP1B', 'LYN', 'LZTR1', 'MAD2L2', 'MAGI2', 'MAP2K1', 'MAP2K2', 'MAP2K4', 'MAP3K1', 'MAX', 
                                 'MC1R', 'MCL1', 'MDM2', 'MDM4', 'MED12', 'MEF2B', 'MEN1', 'MET', 'MITF', 'MLH1', 'MPL', 'MRE11A', 'MSH2', 'MSH3', 'MSH6', 'MTOR', 'MUTYH', 'MYBPC3', 'MYC', 'MYCL', 'MYCL1', 'MYCN', 
                                 'MYD88', 'MYH11', 'MYH7', 'MYL2', 'MYL3', 'NBN', 'NDE1', 'NF1', 'NF2', 'NFE2L2', 'NFKBIA', 'NKX2-1', 'NOTCH1', 'NOTCH2', 'NOTCH3', 'NPM1', 'NRAS', 'NSD1', 'NTHL1', 'NTRK1', 'NTRK2', 
                                 'NTRK3', 'NUP93', 'OTC', 'PAK3', 'PALB2', 'PARK2', 'PAX5', 'PBRM1', 'PCSK9', 'PDCD1LG2', 'PDGFRA', 'PDGFRB', 'PDK1', 'PHOX2B', 'PIK3C2B', 'PIK3CA', 'PIK3CB', 'PIK3CG', 'PIK3R1', 
                                 'PIK3R2', 'PKP2', 'PLCG2', 'PMS2', 'POLD1', 'POLE', 'POT1', 'PPP2R1A', 'PRDM1', 'PREX2', 'PRKAG2', 'PRKAR1A', 'PRKCI', 'PRKDC', 'PRSS1', 'PRSS8', 'PTCH1', 'PTEN', 'PTPN11', 'QKI', 
                                 'RAC1', 'RAD50', 'RAD51', 'RAD51C', 'RAD51D', 'RAF1', 'RANBP2', 'RARA', 'RB1', 'RBM10', 'RECQL4', 'RET', 'RFWD3', 'RICTOR', 'RINT1 ', 'RNF43', 'ROS1', 'RPL11', 'RPL35A', 'RPL5', 
                                 'RPS10', 'RPS17', 'RPS19', 'RPS20', 'RPS24', 'RPS26', 'RPTOR', 'RTEL1', 'RUNX1', 'RUNX1T1', 'RYR1', 'RYR2', 'SBDS', 'SCN5A', 'SDHA', 'SDHAF2', 'SDHB', 'SDHC', 'SDHD', 'SETD2', 
                                 'SF3B1', 'SLCO1B1', 'SLIT2', 'SLX4', 'SMAD2', 'SMAD3', 'SMAD4', 'SMARCA4', 'SMARCB1', 'SMARCE1', 'SMO', 'SNCAIP', 'SOCS1', 'SOX10', 'SOX2', 'SOX9', 'SPEN', 'SPINK1', 'SPOP', 
                                 'SPTA1', 'SRC', 'SRP54', 'SRP72', 'STAG2', 'STAT3', 'STAT4', 'STK11', 'SUFU', 'SYK', 'TAF1', 'TBX3', 'TERC', 'TERT', 'TET2', 'TG', 'TGFB2', 'TGFB3', 'TGFBR1', 'TGFBR2', 'TINF2', 
                                 'TMEM127', 'TMEM43', 'TNFAIP3', 'TNFRSF14', 'TNN13', 'TNNI3', 'TNNT2', 'TOP1', 'TOP2A', 'TP53', 'TPM1', 'TPMT', 'TSC1', 'TSC2', 'TSHR', 'U2AF1', 'UBE2T', 'UGT1A1', 'VEGFA', 'VHL', 
                                 'VKORC1', 'WISP3', 'WRN', 'WT1', 'XPO1', 'XRCC2', 'ZBTB2', 'ZNF217', 'ZNF703', 'POLQ', 'CCNH', 'RAD17', 'EXO1', 'POLN', 'POLK', 'ERCC2', 'XPA', 'ERCC3', 'ERCC6', 'PTCH', 'RNF168'])]
df1['Mutation_Status'] = 'Somatic'
print(df1)
df1 = df1[~df1['Variant_Classification'].isin(['Intron'])]
df1 = df1[~df1['Variant_Classification'].isin(['RNA'])]
df1 = df1[~df1['Variant_Classification'].isin(['Translation_Start_Site'])]
df1 = df1[~df1['Variant_Classification'].isin(['IGR'])]
df1.to_csv('AvatarSomaticACMG_MasterFiltered3Temp.maf', sep='\t', index=False)

#make counter table with no ACMG genes or 5'FLANK
df2 = pd.read_table(('AvatarSomaticNoACMG_MasterFiltered2Temp.maf'), header=None)
#df2 = df2.drop(list(df)[121:124], axis=1)
df2.columns = columns
#remove all 5'Flank/promoter mutations as TERT promoter mutations will be caught in the ACMG list among other un-ingestable 5' promoter mutations
df2 = df2[~df2['Variant_Classification'].isin(["5'Flank"])]
#df2 = df2[df2['FILTER'].isin(['PASS'])]
df2['Mutation_Status'] = 'Somatic'
print(df2)
#remerge ACMG and non-ACMG dfs
df3 = pd.concat([df1, df2], ignore_index=True)
#df3 = df2
df3 = df3.sort_values(by=['Chromosome', 'Start_Position'])
#df3.drop(['t_--tumor-id', 'n_--tumor-id'], axis=1, inplace = True)
#df3 = df3[~df3['IMPACT'].isin(['MODIFIER'])]
df3 = df3[~df3['Variant_Classification'].isin(['Intron'])]
df3 = df3[~df3['Variant_Classification'].isin(['RNA'])]
df3 = df3[~df3['Variant_Classification'].isin(['Translation_Start_Site'])]
df3 = df3[~df3['Variant_Classification'].isin(['IGR'])]
df3 = df3.loc[~df3['CLIN_SIG'].astype(str).isin(['likely_benign,benign/likely_benign', 'association,benign', 'not_provided,benign', 'likely_benign,benign/likely_benign', 'benign,not_provided', 'not_provided,benign', 'benign,not_provided', 'likely_benign,benign/likely_benign', 'benign,association', 'benign,association', 'likely_benign,benign/likely_benign', 'benign,not_provided', 'likely_benign,benign/likely_benign', 'likely_benign,benign/likely_benign', 'likely_benign,benign/likely_benign', 'likely_benign,benign/likely_benign', 'likely_benign,benign/likely_benign', 'benign,not_provided', 'affects,benign', 'not_provided,benign', 'benign,not_provided', 'likely_benign,benign/likely_benign', 'likely_benign,benign/likely_benign', 'likely_benign,benign/likely_benign', 'likely_benign,benign/likely_benign', 'likely_benign,benign/likely_benign', 'likely_benign,benign/likely_benign', 'likely_benign,benign/likely_benign', 'likely_benign,benign/likely_benign', 'likely_benign,benign/likely_benign', 'likely_benign,benign/likely_benign', 'likely_benign,benign/likely_benign', 'benign,not_provided', 'not_provided,benign', 'not_provided,benign', 'likely_benign,benign/likely_benign', 'not_provided,benign', 'not_provided,benign', 'likely_benign,benign/likely_benign,benign', 'likely_benign,benign/likely_benign', 'likely_benign,benign/likely_benign', 'benign,not_provided', 'likely_benign,benign/likely_benign', 'likely_benign,benign/likely_benign', 'likely_benign,benign/likely_benign', 'benign,not_provided', 'likely_benign,benign/likely_benign', 'likely_benign,benign/likely_benign', 'not_provided,benign', 'not_provided,benign', 'benign,not_provided', 'not_provided,benign', 'not_provided,benign', 'not_provided,benign', 'likely_benign,benign/likely_benign'])]
#drop all clin sig benign
#df3 = df3.loc[~df3['CLIN_SIG'].astype(str).str.contains('benign')]
#df3['t_depth'].replace('.', '',  inplace=True)
#df3['t_depth'].replace('.', 0,  inplace=True)
#df3['n_depth'] = df3['t_depth']
#df3['t_depth'] = 0
#df3['n_depth'].replace('.', 0,  inplace=True)
#df3['n_ref_count'].replace('.', 0,  inplace=True)
#df3['n_alt_count'].replace('.', 0,  inplace=True)
#df3['n_ref_count'] = df3['t_ref_count']
#df3['t_ref_count'] = 0
#df3['t_alt_count'].replace('.', 0,  inplace=True)
#df3['n_alt_count'] = df3['t_alt_count']
#df3['t_alt_count'] = 0
#df3.to_csv('Temp.maf', sep='\t', index=False)
#print(df3.dtypes)
#df3[['T_DP']] = df3[['T_DP']].astype(float, errors= 'ignore')
#df3[['T_AF']] = df3[['T_AF']].astype(float, errors= 'ignore')
#df3[['N_AF']] = df3[['N_AF']].astype(float, errors= 'ignore')
df3[['BKZ']] = df3[['BKZ']].astype(float, errors= 'ignore')
#df3['t_ref_count'] = df3['t_ref_count'].astype(float, errors= 'ignore')
#print(df3.dtypes) 
#df3['t_depth'] = df3['T_DP'].astype(float, errors= 'ignore')
#df3['t_alt_count'] = df3['T_DP'] * df3['T_AF']
#df3['t_alt_count'] = np.ceil(df3['t_alt_count']).astype(int)
#df3['t_ref_count'] = df3['t_depth'] - df3['t_alt_count']
#df3['n_depth'] = df3['N_DP'].astype(float, errors= 'ignore')
#df3['n_alt_count'] = df3['N_DP'] * df3['N_AF']
#df3['n_alt_count'] = np.ceil(df3['n_alt_count']).astype(int)
#df3['n_ref_count'] = df3['n_depth'] - df3['n_alt_count']

df3 = df3.astype({"t_depth": float, "t_ref_count": float,"t_alt_count": float})
df3a = df3[(df3['t_depth'] >12) & (df3['t_ref_count'] >8) & (df3['t_alt_count'] >4)& (df3['BKZ'] >4)]
df3a.drop(['T_DP', 'T_AF', 'N_DP', 'N_AF', 'BKZ'], axis=1, inplace = True)
df3a['NCBI_Build'] = 'hg38'
#df3['Center'] = 'HCI'
#df3a['Tumor_Sample_Barcode'] = df3a['Matched_Norm_Sample_Barcode']
#df3a['Matched_Norm_Sample_Barcode'] = df3a['Tumor_Sample_Barcode'].str.split('-').str[0]
df3a.drop_duplicates(subset=['Tumor_Sample_Barcode', 'Hugo_Symbol', 'Start_Position', 'End_Position'], inplace=True)
df3a['Entrez_Gene_Id'] = df3a['Entrez_Gene_Id'].astype(str)
df3a = df3a[~df3a['Entrez_Gene_Id'].isin(['0'])]
#update Sample names
#SL_SAMPLE = pd.read_csv('SL_SAMPLE.txt', sep='\t')
#df3a = df3a.merge(SL_SAMPLE, on='Matched_Norm_Sample_Barcode', how='left')
#df3a['Tumor_Sample_Barcode'] = df3a['SAMPLE_ID']
#df3a.drop(['SAMPLE_ID'], axis=1, inplace = True)
df3a.to_csv('SomaticACMG_MasterFilteredSortedFinal3.maf', sep='\t', index=False)
df4 = df3a['Hugo_Symbol'].value_counts(dropna=False).to_frame()
#print(df4)
df5 = df3a['dbSNP_RS'].value_counts(dropna=False).to_frame()
df4.to_csv('GeneCount.txt', sep='\t')
df5.to_csv('RSCount.txt', sep='\t')