# -*- coding: utf-8 -*-
"""
Created on Mon Oct 14 19:32:00 2024

@author: Buyobillion
"""

import os
import subprocess
import networkx as nx
import pandas as pd

e2o_path = r"F:\result\imp\imp\e2o"
o2e_path = r"F:\result\imp\imp\o2e"
path = e2o_path
e2o_out = r"F:\result\imp\imp\e2o_set"
o2e_out = r"F:\result\imp\imp\o2e_set"
new_path = e2o_out

snp = {}
for file in os.listdir(path): 
    file_path = os.path.join(path, file)
    if "None" not in file: 
        df = pd.read_table(file_path, sep='\t')
        for name, data in df.groupby('CHR'): 
            if not name in snp: 
                snp[name] = set()
            snp[name].update(data['SNP'].tolist())
for chr_num in range(1, 23): 
    if chr_num not in snp:
        continue
    now_new_path = os.path.join(new_path, f'chr{chr_num}')
    os.makedirs(now_new_path, exist_ok = True)
    mapping_path = rf"F:\phasing_rsid_snpid_mapping\rsid_snpid_table\rsid_snpid_table_chr{chr_num}.csv"
    mapping = pd.read_table(mapping_path)
    set_snp = []
    set_snp = [i for i in snp[chr_num] if i in mapping['rsid'].values]
    missed_snp = [i for i in snp[chr_num] if i not in set_snp]
    missed_df = pd.DataFrame(list(missed_snp), columns = ['missed_snp'])
    missed_df.to_csv(os.path.join(now_new_path, f'chr{chr_num}_missed_snp.csv'), sep='\t')
    set_rsid_snpid = mapping[mapping['rsid'].isin(set_snp)]
    set_rsid_snpid.to_csv(os.path.join(now_new_path, f'chr{chr_num}_kept_snp.csv'), sep='\t')
    set_snpid = []
    set_snpid = set(mapping[mapping['rsid'].isin(set_snp)]['snpid'].tolist())
    set_snpid_df = pd.DataFrame(list(set_snpid), columns = ['snpid'])
    set_snpid_df.to_csv(os.path.join(now_new_path, f'chr{chr_num}_snpid.csv'), sep='\t', index=False, header=False)

data_dir = r"F:\phasing_rsid_snpid_mapping\Grch38_bed"
ref_path = r"F:\phasing_rsid_snpid_mapping\Grch38_bed\EUR.txt"

for chr_num in range(1, 23): 
    chr_folder_path = os.path.join(new_path, f"chr{chr_num}")
    if not os.path.exists(chr_folder_path):
        print(f"chr{chr_num} not found in target folder. Skipping...")
        continue

    snp_path = os.path.join(new_path, rf"chr{chr_num}\chr{chr_num}_snpid.csv")
    output_path = os.path.join(new_path, rf"chr{chr_num}")
    output_file = os.path.join(new_path, rf"chr{chr_num}\chr{chr_num}_set_ld")
    
    data_file = f"Grch38_chr{chr_num}"
    existing_files = {f for f in os.listdir(output_path) if f.endswith('.ld')}    
    input_bfile = os.path.join(data_dir, data_file)
    if f"{output_file}.ld" in existing_files: 
        print(f"Output file {output_file}.ld already exists. Skipping...")
        continue
              
    command = f"plink --bfile {input_bfile} " \
          f"--extract {snp_path} " \
          f"--keep {ref_path} " \
          f"--r2 " \
          f"--ld-window 999999999 " \
          f"--ld-window-kb 999999999 " \
          f"--ld-window-r2 0 " \
          f"--maf 0.01 " \
          f"--out {output_file}"                   
    print(f"Running command: {command}")
    subprocess.run(command, shell=True)
    print("Processing completed")
    
e2o_path = r"F:\result\imp\imp\e2o"
o2e_path = r"F:\result\imp\imp\o2e"
path = o2e_path
e2o_out = r"F:\result\imp\imp\e2o_set"
o2e_out = r"F:\result\imp\imp\o2e_set"
new_path = o2e_out

e2o_out_path = r"F:\result\imp\imp\post_ld\e2o"
o2e_out_path = r"F:\result\imp\imp\post_ld\o2e"
out_path = o2e_out_path

clusters = {}
ld = {}
mapping = {}

for chr_num in range(1, 23): 
    chr_folder_path = os.path.join(new_path, f"chr{chr_num}")
    if not os.path.exists(chr_folder_path):
        print(f"chr{chr_num} not found in target folder. Skipping...")
        continue

    mapping_path = os.path.join(chr_folder_path, f"chr{chr_num}_kept_snp.csv")
    mapping[chr_num] = pd.read_table(mapping_path)

    ld_path = os.path.join(chr_folder_path, f"chr{chr_num}_set_ld.ld")
    if os.path.exists(ld_path):
        ld[chr_num] = pd.read_table(ld_path, sep=r'\s+')
        ld[chr_num] = ld[chr_num][ld[chr_num]['R2'] > 0.2]
        
        G = nx.Graph()
        for index, row in ld[chr_num].iterrows():
            G.add_edge(row['SNP_A'], row['SNP_B'])
        clusters[chr_num] = list(nx.connected_components(G))
    else:
        print(f"{ld_path} not found. Adding single SNP cluster from kept SNP file.")
        single_snp_data = mapping[chr_num].iloc[0]  # Assuming the first row contains the SNP info
        cluster_data = {
            'chr': single_snp_data['chr'],
            'pos': single_snp_data['pos'],
            'rsid': single_snp_data['rsid'],
            'snpid': single_snp_data['snpid'],
            'A1': single_snp_data['A1'],
            'A2': single_snp_data['A2']
        }
        clusters[chr_num] = [{cluster_data['snpid']}]
        
for file in os.listdir(path):
    out_file = os.path.join(out_path, file)
    df = pd.DataFrame()
    if "None" in file: 
        df.to_csv(out_file, sep='\t')
        print(f"Empty file {out_file}. Skipping...")
        continue
    existing_files = {f for f in os.listdir(out_path)}
    if f"{file}" in existing_files: 
        print(f"Output file {out_file} already exists. Skipping...")
        continue
    print(f"Processing file {out_file} ...")
    df = pd.read_table(os.path.join(path, file)) 
    df_merged = pd.DataFrame()
    for chr_num in range(1, 23): 
        chr_folder_path = os.path.join(new_path, f"chr{chr_num}")
        if not os.path.exists(chr_folder_path):
            print(f"chr{chr_num} not found in target folder. Skipping...")
            continue
        df_part = pd.DataFrame()
        df_part = pd.merge(df[df['CHR'] == chr_num], mapping[chr_num][['rsid', 'snpid']], left_on='SNP', right_on='rsid', how='left')
        clus = []
        snps = []
        df_part_snp_list = df_part['snpid'].tolist()
        for cluster in clusters[chr_num]:
            snps = [snp for snp in cluster if snp in df_part_snp_list]
            if snps:  
                min_snp = min(snps, key=lambda snp: df_part.loc[df_part['snpid'] == snp, 'pval.outcome'].values[0])
                clus.append(min_snp) 
        df_part = df_part[df_part['snpid'].isin(clus)]
        df_merged = pd.concat([df_merged, df_part])
        
    df_merged.to_csv(out_file, sep='\t')
    print(f"Successvfully finished processing file {out_file}. ")

df = pd.read_csv(r"F:\result\post_ld\e2o_2401-3935\13_EA1_multi_p5e-8_sumstats_to_2835.txt.csv", sep='\t')
df[df['CHR']==3]['SNP'].to_csv(r"F:\result\post_ld\verify.csv", sep='\t', header=False, index=False)
