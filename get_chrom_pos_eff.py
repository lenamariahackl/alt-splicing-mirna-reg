from helper_fcts import *
import os
import pandas as pd
from pathlib import Path
import numpy as np
import sys

test = False
path = Path('/nfs/home/users/l.hackl/alt-splicing-mirna-reg/data')
tcga_path = Path(path/'KICH')
path_tarp = Path(path/'tarp-bs') if test else Path('../data')
path_input = Path(path/'tarp-bs/cdna') if test else Path('../data/all_cdna_subsets')
print('start')
sys.stdout.flush()
exon_positions = read_in_exon_pos(path/'exon_positions.csv') # mapping of transcript ids to exon ids, chrom_exon_starts, exon_starts, exon_ends
bs = pd.read_parquet(path/'bs_all.parquet') #983499270 rows

#split bs in 10 parts
i=9
bs = bs[(i-1)*100000000:i*100000000]
print('read bs')
sys.stdout.flush()
chrom_pos(bs, path/('bs_chrom_pos'+str(i)+'.parquet'), exon_positions)
print('saved parquet',i)
sys.stdout.flush()
bs = 0
#10. part
bs = pd.read_parquet(path/'bs_all.parquet') #983499270 rows
bs = bs[900000000:983499269]
print('read bs')
sys.stdout.flush()
chrom_pos(bs, path/('bs_chrom_pos10.parquet'), exon_positions)
print('saved parquet 10')
sys.stdout.flush()