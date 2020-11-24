import os
import pandas as pd
from pathlib import Path
import numpy as np

#TODO in case bs chrom pos trans has finished
#run in script for bs already
all_chrom = []
for i in range(1,10):
    all_chrom.append(pd.read_parquet(path/'bs_chrom_pos'+str(i)+'.parquet'))
chrom = pd.concat(all_chrom, axis=1)
chrom.to_parquet(path/'bs_all_chrom.parquet')
chrom = chrom[['miRNA','mRNA','chrom_bs_start','chrom_bs_end','exon_id']]
chrom['miRNA'] = pd.Categorical(chrom.miRNA)
chrom['mRNA'] = pd.Categorical(chrom.mRNA)
chrom.to_parquet(path/'bs_subset_chrom.parquet')