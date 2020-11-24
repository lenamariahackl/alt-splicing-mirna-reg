from helper_functions import *
import tarpmir_to_pandas
import translate_pos_to_chrom
import merge_bs_with_miRNA

import sqlite3
import os
import pandas as pd
from pathlib import Path
import numpy as np

def main(path, path_tarp):
    log_path = Path(path/'log')
    tcga_path = Path(path/'KICH')
    path_input = Path(path_tarp/'all_cdna_subsets')
    #TODO
    
    #0) read in TarPmiR output files to one pandas DataFrame
    parquet_file = Path(path/'bs_all.parquet')
    if parquet_file.is_file():
        bs_tarp = pd.read_parquet(parquet_file)
    else:
        bs_tarp = tarpmir_to_pandas.main(path_tarp, path, log_path) #DONE
    #1) translate binding site position from transcript to chromosome
    bs_tarp = translate_pos_to_chrom.main(path, bs_tarp, log_path)
    #2) merge to add chromosome and strand to binding sites using the transcript id
    bs_tarp = merge_bs_with_transcripts.main(path, bs_tarp, log_path)
    
    #1 + 2) get exon and miRNA counts from XENA
    mirna_counts, exon_counts = counts_to_pandas.main(path, log_path)
    #3) outer join exon counts with mirna counts
    bs_counts = merge_exons_with_mirnas.main(path, mirna_counts, exon_counts, log_path)
    #4) Elastic Net Regression
    bs_counts = elastic_net.main(path, bs_counts, log_path)
    
    # Filter bindingsites
    bs = filter_bs.main(path, bs_tarp, bs_counts, log_path)
    # Analyse alternative splicing events using DIGGER
    
    
if __name__ == '__main__':
    path = Path('/nfs/home/users/l.hackl/alt-splicing-mirna-reg/data')
    path_tarp = Path('/nfs/home/users/l.hackl/data')
    main(path, path_tarp)

#TODO put this somewhere

parquet_file = Path(path/'bs_subset.parquet')
if parquet_file.is_file():
    #load file from disk
    bs_new = pd.read_parquet(parquet_file)
else:
    #make small bs subset
    bs_new = bs[['miRNA','mRNA','bs_start','bs_end']].copy()
    bs_new['miRNA'] = pd.Categorical(bs_new.miRNA)
    bs_new['mRNA'] = pd.Categorical(bs_new.mRNA)
    bs_new.info()
    bs_new.to_parquet(parquet_file)

#parse miRNA and mRNA seed from sequence
seperate_mRNA = []
for filename in os.listdir(path_input):
    seperate_mRNA.append(pd.DataFrame.from_dict(parse_seq(path_input/filename)))
mRNA_sequences = pd.concat(seperate_mRNA, axis=0, ignore_index=True)
mRNA_sequences.to_parquet(path/'mRNA_all.parquet')
print('All input mRNAs were read into pandas.')

exon_positions = read_in_exon_pos(path/'exon_positions.csv') # mapping of transcript ids to exon ids, chrom_exon_starts, exon_starts, exon_ends

#split bs in 10 parts
for i in range(1,9):
    bs = pd.read_parquet(path/'bs_all.parquet') #983499270 rows
    bs = bs[(i-1)*100000000:i*100000000]
    print('read bs')
    chrom_pos(bs, path/('bs_chrom_pos'+str(i)+'.parquet'), exon_positions)
    print('saved parquet ',i)
    
#10. part
bs = pd.read_parquet(path/'bs_all.parquet') #983499270 rows
bs = bs[900000000:983499269] #TODO is this right?
print('read bs')
chrom_pos(bs, path/('bs_chrom_pos10.parquet'), exon_positions)
print('saved parquet 10')

