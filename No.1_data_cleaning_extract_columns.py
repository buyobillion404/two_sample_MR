# -*- coding: utf-8 -*-
"""
Created on Thu Sep 26 22:26:11 2024

@author: Administrator
"""

import os 
import pandas as pd
import gzip


path = r"F:\result\16_tan"
new_path = r"F:\result\imp\16_tan"
def extract_chosen_columns(path): 
    
    SNP = ['SNP', 'SNPID', 'rsID', 'MarkerName', 'rsid', 'variant_id']
    CHR = ['CHR', 'Chr', 'chromosome', 'Chromosome']
    A1 = ['A1', 'EFFECT_ALLELE', 'Effect_allele', 'alt', 'effect_allele']
    A2 = ['A2', 'OTHER_ALLELE', 'Other_allele', 'ref', 'other_allele']
    Beta = ['Beta', 'BETA', 'direct_Beta', 'stdBeta', 'beta']
    SE = ['SE', 'direct_SE', 'standard_error']
    Pval = ['Pval', 'PVALUE', 'P', 'direct_log10_P', 'p_value', 'direct_pval']
    all_columns = SNP + CHR + A1 + A2 + Beta + SE + Pval
    
    for root, dirs, files in os.walk(path):
        for filename in files: 
            file_path = os.path.join(root, filename)
            df = None 
            
            if filename.endswith(".gz"): 
                with gzip.open(file_path, "rt") as f: 
                    df = pd.read_csv(f, sep=r'\s+')
            else: 
                df = pd.read_table(file_path)
            
            extracted_columns = [col for col in all_columns if col in df.columns]
            if extracted_columns: 
                extracted_df = df[extracted_columns]
                
                if not extracted_df.empty:  
                    relative_path = os.path.relpath(root, path)
                    new_dir = os.path.join(new_path, relative_path)
                    os.makedirs(new_dir, exist_ok=True)
                    output_file_path = os.path.join(new_dir, os.path.splitext(filename)[0] + '.csv')
                    extracted_df.to_csv(output_file_path, sep='\t', index=False)

extract_chosen_columns(path)

new_path = r'F:\result\imp\16_tan\rename'
def rename_files(path): 
    SNP = ['SNP', 'SNPID', 'rsID', 'MarkerName', 'rsid', 'variant_id']
    CHR = ['CHR', 'Chr', 'chromosome', 'Chromosome']
    A1 = ['A1', 'EFFECT_ALLELE', 'Effect_allele', 'alt', 'effect_allele']
    A2 = ['A2', 'OTHER_ALLELE', 'Other_allele', 'ref', 'other_allele']
    Beta = ['Beta', 'BETA', 'direct_Beta', 'stdBeta', 'beta']
    SE = ['SE', 'direct_SE', 'standard_error']
    Pval = ['Pval', 'PVALUE', 'P', 'direct_log10_P', 'p_value', 'direct_pval']

    column_mapping = {}

    for new_names in [SNP, CHR, A1, A2, Beta, SE, Pval]:
        new_name = new_names[0]  
        for old_name in new_names:
            column_mapping[old_name] = new_name
            
    for root, dirs, files in os.walk(path): 
        for filename in files: 
            file_path = os.path.join(root, filename)
            df = pd.read_table(file_path)
            df = df.rename(columns = column_mapping)
            relative_path = os.path.relpath(root, path)
            new_dir = os.path.join(new_path, relative_path)
            os.makedirs(new_dir, exist_ok=True)
            output_file_path = os.path.join(new_dir, os.path.splitext(filename)[0] + '.csv')
            df.to_csv(output_file_path, sep='\t', index=False)
rename_files(r'F:\result\imp\16_tan')

def log10_P_to_P(path, var): 
    for root, dirs, files in os.walk(path): 
        for file in files: 
            file_path = os.path.join(root, file)
            df = pd.read_table(file_path)
            df[var] = 10 ** (-df[var])
            file_new = os.path.splitext(file_path)[0] + '_new.csv'
            output_file_path = os.path.join(root, file_new)
            df.to_csv(output_file_path, sep='\t', index=False)
log10_P_to_P(r'F:\exposure_in_format_csv_sep_tab\15_EA&Cognitive ability_Young et al', 'direct_log10_P')

def snp_filter(path, new_path): 
    for root, dirs, files in os.walk(path): 
        for filename in files: 
            file_path = os.path.join(root, filename)
            df = pd.read_table(file_path)
            df = df[df['SNP'].str.contains('rs')]
            relative_path = os.path.relpath(root, path)
            new_dir = os.path.join(new_path, relative_path)
            os.makedirs(new_dir, exist_ok=True)
            output_file_path = os.path.join(new_dir, os.path.splitext(filename)[0] + '.csv')
            df.to_csv(output_file_path, sep='\t', index=False)
snp_filter(r'F:\exposure_in_format_csv_sep_tab', r'F:\exposure_in_format_csv_sep_tab_rs_filtered')

def read_column_names(file_path):
    if file_path.endswith('.gz'):
        with gzip.open(file_path, 'rt') as f:
            header = f.readline().strip()
    else:
        with open(file_path, 'r') as f:
            header = f.readline().strip()
    column_names = header.split()  
    return column_names, len(column_names)

a = 0
for root, dirs, files in os.walk(new_path):
    for filename in files: 
        file_path = os.path.join(root, filename)
        columns, colnum = read_column_names(file_path)
        print(f"{filename}\n{columns}\n{colnum}")
        a += 1
