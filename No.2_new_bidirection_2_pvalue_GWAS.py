# -*- coding: utf-8 -*-
"""
Created on Thu Jul 25 10:12:04 2024

@author: Buyobillion
"""

import re
import gzip
import pandas as pd
import numpy as np
import os
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
import multiprocessing
from warnings import simplefilter
simplefilter(action = "ignore", category=FutureWarning())
import warnings
warnings.filterwarnings("ignore", category=pd.errors.SettingWithCopyWarning)
from io import StringIO


SNP = ['SNP', 'SNPID', 'rsID', 'MarkerName', 'rsid', 'variant_id']
CHR = ['CHR', 'Chr', 'chromosome', 'Chromosome']
A1 = ['A1', 'EFFECT_ALLELE', 'Effect_allele', 'alt', 'effect_allele']
A2 = ['A2', 'OTHER_ALLELE', 'Other_allele', 'ref', 'other_allele']
Beta = ['Beta', 'BETA', 'direct_Beta', 'stdBeta', 'beta']
SE = ['SE', 'direct_SE', 'standard_error']
Pval = ['Pval', 'PVALUE', 'P', 'direct_log10_P', 'p_value', 'direct_pval']


dicColExporsure, dicColOutcome = {}, {}
for ColumnsPre, ColumnExposure, ColumnOutcome in zip([SNP, A1, A2, Beta, SE, Pval],
                                ['SNP', 'A1.exposure', 'A2.exposure', 'beta.exposure', 'se.exposure', 'pval.exposure'],
                                ['SNP', 'A1.outcome', 'A2.outcome', 'beta.outcome', 'se.outcome', 'pval.outcome']):
    for col in ColumnsPre:
        dicColExporsure[col], dicColOutcome[col] = ColumnExposure, ColumnOutcome

path_outcome = r'F:\GWAS_brainstructure_1-3935'
path_exposure = r"F:\result\imp\16_tan\exposure\rename"

path_exposure_to_outcome = r'F:\result\imp\16_tan\e2o'
path_outcome_to_exposure = r'F:\result\imp\16_tan\o2e'

def initialize_exposure(path):
    exposure_dict = {}
    for root, dirs, files in os.walk(path):
        for filename in files: 
            file_path = os.path.join(root, filename)
            exposure_dict[filename] = pd.read_csv(file_path, sep = '\t')
    return exposure_dict

exposure_dict = initialize_exposure(path_exposure)
exposure_SNP = []
for _, data in exposure_dict.items():
    exposure_SNP += data[data.Pval < 5e-8]["SNP"].values.tolist()
exposure_SNP = list(set(exposure_SNP))


def initialize_outcome(path, exposure_SNP):
    try:
        with gzip.open(path, 'rt') as f:
            df_brain = pd.read_csv(f, sep= ' ', usecols=['rsid', 'a1', 'a2', 'beta', 'se', 'pval(-log10)'])
    except:
        return os.path.basename(path), 'Failed'
    df_brain["Pval"] = np.power(10, -1 * df_brain["pval(-log10)"])
    df_brain = df_brain.drop('pval(-log10)', axis=1)
    df_brain = df_brain[(df_brain["rsid"].isin(exposure_SNP)) | (df_brain['Pval'] < 5e-8)]
    return os.path.basename(path), df_brain

def Exposure_to_Outcome(dic_Exposure, dic_Outcome, p_threshold, new_path):
    for fileOutcome, dataOutcome in dic_Outcome.items(): 
        dataOutcome = dataOutcome[dataOutcome.Pval > p_threshold]
        for fileExposure, dataExposure in dic_Exposure.items():
            dataExposure = dataExposure[dataExposure.Pval < p_threshold]  
            if len(dataExposure) == 0 or len(dataOutcome) == 0:
                data_new = pd.DataFrame()
                data_new.to_csv(os.path.join(new_path, f'None_{Path(fileExposure).stem}_to_{Path(fileOutcome).stem}.csv'), sep='\t', index=False)
                continue
            
            for column_exposure in [x for x in dataExposure.columns if x in dicColExporsure.keys()]:
                dataExposure.rename(columns={column_exposure: dicColExporsure[column_exposure]}, inplace=True)
            for column_outcome in [x for x in dataOutcome.columns if x in dicColOutcome.keys()]:
                dataOutcome.rename(columns={column_outcome: dicColOutcome[column_outcome]}, inplace=True)

            data_new = pd.merge(dataExposure, dataOutcome, on="SNP")
            
            if data_new.empty:
                data_new.to_csv(os.path.join(new_path, f'None_{Path(fileExposure).stem}_to_{Path(fileOutcome).stem}.csv'), sep='\t', index=False)
            else:
                data_new.to_csv(os.path.join(new_path, f'{Path(fileExposure).stem}_to_{Path(fileOutcome).stem}.csv'), sep='\t', index=False)

def Outcome_to_Exposure(dic_Outcome, dic_Exposure, p_threshold, new_path):
    for fileOutcome, dataOutcome in dic_Outcome.items(): 
        dataOutcome = dataOutcome[dataOutcome.Pval < p_threshold]
        for fileExposure, dataExposure in dic_Exposure.items():
            if isinstance(dataOutcome, str): 
                Fail_brainstructure.append(fileOutcome)
                continue
            
            dataExposure = dataExposure[dataExposure.Pval > p_threshold]  
            if len(dataExposure) == 0 or len(dataOutcome) == 0:
                data_new = pd.DataFrame()
                data_new.to_csv(os.path.join(new_path, f'None_{Path(fileOutcome).stem}_to_{Path(fileExposure).stem}.csv'), sep='\t', index=False)
                continue
            
            for column_exposure in [x for x in dataExposure.columns if x in dicColExporsure.keys()]:
                dataExposure.rename(columns={column_exposure: dicColExporsure[column_exposure]}, inplace=True)
            for column_outcome in [x for x in dataOutcome.columns if x in dicColOutcome.keys()]:
                dataOutcome.rename(columns={column_outcome: dicColOutcome[column_outcome]}, inplace=True)

            data_new = pd.merge(dataExposure, dataOutcome, on="SNP")
            
            if data_new.empty:
                data_new.to_csv(os.path.join(new_path, f'None_{Path(fileOutcome).stem}_to_{Path(fileExposure).stem}.csv'), sep='\t', index=False)
            else:
                data_new.to_csv(os.path.join(new_path, f'{Path(fileOutcome).stem}_to_{Path(fileExposure).stem}.csv'), sep='\t', index=False)

#Exposure_to_Outcome(exposure_dict, dict, 5e-8, path_exposure_to_outcome)

#Outcome_to_Exposure(dict, exposure_dict, 5e-8, path_outcome_to_exposure)



Fail_brainstructure = []
for project_index in [i for i in range(1, 3935, 10)]:

    dict_brain = {}
    brain_num_llm = project_index
    brain_num_ulm = project_index + 10
    print(f"Processing files: {brain_num_llm} to {brain_num_ulm-1}")
     
    paths_brain = [os.path.join(path_outcome, f) for f in os.listdir(path_outcome) if int(f.replace(".txt.gz","")) in range(brain_num_llm, brain_num_ulm)]
    
    with ThreadPoolExecutor(max_workers=10) as executor:
        futures = {executor.submit(initialize_outcome, path, exposure_SNP): path for path in paths_brain}
        for future in as_completed(futures): 
            file_name, data_brain = future.result()
            if isinstance(data_brain, str): 
                Fail_brainstructure.append(file_name)
                continue
            dict_brain[file_name] = data_brain
            print(f"Brain GWAS finish processing {file_name} with {len(data_brain)} SNPs")
    
    Exposure_to_Outcome(exposure_dict, dict_brain, 5e-8, path_exposure_to_outcome)
    Outcome_to_Exposure(dict_brain, exposure_dict, 5e-8, path_outcome_to_exposure)
    
    print(f"Finished files: {brain_num_llm} to {brain_num_ulm-1}")
    del dict_brain

print(Fail_brainstructure)