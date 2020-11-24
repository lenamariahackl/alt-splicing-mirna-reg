import csv
import gc
# read exon position dictionary from file
# returns format {transcript_id : {exon_id : (exon_start, exon_end, chromosome_exon_start)}}
def read_in_exon_pos(file):
    exon_positions = {}
    with open(file,'r') as data: 
        for line in csv.reader(data): 
            tid, eid, exon_start, exon_end, chrom_exon_start = line
            exon_positions[tid] = exon_positions.get(tid,{})
            exon_positions[tid][eid] = (int(exon_start),int(exon_end),int(chrom_exon_start))
    sorted_dict1 = {}
    for t in exon_positions:
        dict1 = exon_positions[t]
        sorted_values = sorted(dict1.values()) # Sort the values
        sorted_dict = {}
        for i in sorted_values:
            for k in dict1.keys():
                if dict1[k] == i:
                    sorted_dict[k] = dict1[k]
                    break
        sorted_dict1[t] = sorted_dict1.get(t,{})
        sorted_dict1[t].update(sorted_dict)
    return sorted_dict1

def calc_chrom_pos(start, end, exon_positions):
    start_chrom_pos = False
    eids = []
    for eid, data in exon_positions.items():
        exon_start = data[0]
        exon_end = data[1]
        chrom_exon_start = data[2]
        if (start >= exon_start) and (start <= exon_end):
            start_chrom_pos = chrom_exon_start + (start - exon_start)
            if (end >= exon_start) and (end <= exon_end):
                end_chrom_pos = chrom_exon_start + (end - exon_start)
                return start_chrom_pos, end_chrom_pos, eid, 1
            else: eids.append(eid)   
        elif (end >= exon_start) and (end <= exon_end):
            end_chrom_pos = chrom_exon_start + (end - exon_start)
            eids.append(eid)
            return start_chrom_pos, end_chrom_pos, ','.join(eids), len(eids)
        elif start_chrom_pos: #in case of bindingsites spanning > 2 exons
            eids.append(eid)
    return None, None, None, None

def chrom_pos(df, file, exon_positions):
    v1s, v2s, v3s, v4s = [], [], [], []
    i = 0
    for _, row in df.iterrows():
        v1, v2, v3, v4 = calc_chrom_pos(row.bs_start, row.bs_end, exon_positions[row.mRNA])
        if i % 5000 == 0:
            gc.collect()
        v1s.append(v1)
        v2s.append(v2)
        v3s.append(v3)
        v4s.append(v4)
        i += 1
    print('Done computation')
    df_result = pd.DataFrame({'chrom_bs_start': v1s,
                            'chrom_bs_end': v2s,
                            'exon_id': v3s,
                            'in_exon': v4s})
    df.reset_index(drop=True, inplace=True)
    df_result.reset_index(drop=True, inplace=True)
    df = pd.concat([df,df_result], axis=1)
    df.to_parquet(file)
    gc.collect()
    return df

chrom = pd.read_parquet(path/'chrom'/'bs_na_chrom.parquet')[1:1000] #19 GB
exon_positions = read_in_exon_pos(path/'exon_positions.csv') # mapping of transcript ids to exon ids, chrom_exon_starts, exon_starts, exon_ends
a = chrom_pos(chrom, path/'chrom'/'bs_tarp_multi.parquet', exon_positions) #look at mem usage 