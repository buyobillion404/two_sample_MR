# -*- coding: utf-8 -*-
"""
Created on Wed Oct 16 09:52:48 2024

@author: Buyobillion
"""

import os 
import pandas as pd

out_path = r"F:\phasing\rsid_snpid_table"
rec = pd.DataFrame(columns=['chr', 'rs_len', 'bim_len', 'merged_len', 'rm_len'])
for chr_num in range(1, 23): 
    path_rs = rf"F:\phasing\dsSNP_mapping_ref\rsid_table_chr{chr_num}.csv"
    path_bim = rf"F:\phasing\Grch38_bed\Grch38_chr{chr_num}.bim"

    rs = pd.read_table(path_rs, header=None)
    bim = pd.read_table(path_bim, sep=r'\s+', header=None)
    rs.columns = ["chr", "pos", "rsid", "A2", "A1"]
    bim.columns = ["chr_cp", "snpid", "uk", "pos", "A1_cp", "A2_cp"]

    df = pd.merge(rs, bim, on='pos', how='inner')
    df = df[["chr", "pos", "rsid", "snpid", "A1", "A2"]]
    rm_snps = rs[~rs['pos'].isin(df['pos'])]
    new_row = pd.DataFrame([[f'chr{chr_num}', len(rs), len(bim), len(df), len(rm_snps)]], columns=['chr', 'rs_len', 'bim_len', 'merged_len', 'rm_len'])
    rec = pd.concat([rec, new_row], ignore_index=True)
    df.to_csv(os.path.join(out_path, f"rsid_snpid_table_chr{chr_num}.csv"), sep='\t', index=False)
    rm_snps.to_csv(os.path.join(out_path, f"removed_rsid_chr{chr_num}.csv"), sep='\t', index=False)
    
rec.to_csv(os.path.join(out_path, "len_table_rec.csv"), sep='\t', index=False)