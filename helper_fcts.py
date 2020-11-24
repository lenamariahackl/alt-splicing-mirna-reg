import re
import csv
import pandas as pd 
from pathlib import Path
import pickle
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
	return exon_positions

# returns format {transcript_id : {exon_id : data}}
def read_in_fasta(file):
	file = open(file,'r').readlines()
	dict_ = {}
	for id_, line in enumerate(file):
		if line[0] == '>': 	# data is in next line
			transcript_id = re.search('ENST[^_]*',line).group(0)
			exon_id = re.search('ENSE[^\n]*',line).group(0)
			data = file[id_+1][:-1]
			data = int(data) if data != 'NA' else None
			dict_[transcript_id] = dict_.get(transcript_id, {})
			dict_[transcript_id][exon_id] = data
	return dict_

# returns format {transcript_id : data} or {transcript_id : {exon_id : data}}
def read_in_miRNA(file):
	pickle1 = Path(str(file).split('.', 1)[0]+".pkl")
	if pickle1.is_file():
		#load from disk
		pkl_file1 = open(pickle1, 'rb')
		dict_ = pickle.load(pkl_file1)
		pkl_file1.close()
	else:
		file = open(file,'r').readlines()
		dict_ = {}
		for id_, line in enumerate(file):
			if line[0] == '>':
				miRNA_id = line[1:-1]
				# data is in next line
				data = file[id_+1][:-1]
				dict_[miRNA_id] = data
		#write to disk
		f = open(pickle1,"wb")
		pickle.dump(dict_,f)
		f.close()
	return dict_

# add header to beginning of file
def append_header(fasta):
	header = 'miRNA	mRNA	binding_site	binding_probability	energy	seed	accessibility	AU_content	PhyloP_Stem	PyloP_Flanking	m/e 	number_of_pairings	binding_region_length	longest_consecutive_pairings	position_of_longest_consecutive_pairings	pairings_in_3prime_end	difference_of_pairings_between_seed_and_3prime_end\n'
	f = open(fasta,'r+')
	lines = f.readlines()
	f.seek(0)
	f.write(header)
	for line in lines:
	    f.write(line)
	f.close()

# read in binding site file ( file that is produced by tarpmir )
# returns pandas DataFrame
def parse_tarp_bs(file):
	if (open(file).readline()[0:5] != 'miRNA'):
		append_header(file)
	binding_sites = pd.read_csv(file, delimiter='	')
	# 1 column binding_sites to 2 new columns bs_start, bs_end
	new = binding_sites['binding_site'].str.split(',', 1, expand=True)
	binding_sites['bs_start'] = new[0].astype(int)
	binding_sites['bs_end'] = new[1].astype(int)
	binding_sites.drop(columns =["binding_site"], inplace = True) 
	return binding_sites

#returns format {name : sequence}
def parse_seq(file):
	file = open(file,'r').readlines()
	seqs = {}
	for i, line in enumerate(file):
		if line[0] == '>':
			name = line[1:-1]
			sequence = file[i+1][:-1]
			seqs[name] = sequence
	return seqs

#returns format {transcript_id : {exon_id : chromosome_start}}, {transcript_id : {exon_id : start}}, {transcript_id : {exon_id : end}}, pandas DataFrame [exon_id, transcript_id, chrom_exon_start, chrom_exon_end]
def calc_exon_data(path):
	feather_file = Path(path/'exon_info.feather')
	if feather_file.is_file():
		#load file from disk
		exon_info = pd.read_feather(feather_file)
	else:
		chrom_exon_starts = read_in_fasta(path/'all_exon_start.fasta')
		chrom_exon_ends = read_in_fasta(path/'all_exon_end.fasta')
		exon_info = {}
		for tid in chrom_exon_starts:
			sorted_ = {k: v for k, v in sorted(chrom_exon_starts[tid].items(), key=lambda item: item[1])}
			for eid in sorted_:
				chrom_exon_start = chrom_exon_starts[tid][eid]
				chrom_exon_end = chrom_exon_ends[tid][eid]
				exon_info[eid] = [tid, chrom_exon_start, chrom_exon_end]
				#write to disk
		exon_info = pd.DataFrame.from_dict(exon_info, orient='index', columns=['transcript_id','chrom_exon_start', 'chrom_exon_end']).reset_index()
		exon_info = exon_info.rename(columns={'index': "exon_id"})
		exon_info.to_feather(feather_file)
	return exon_info

def get_chrom_pos(start, end, tid, exon_positions):#old
	for eid, data in exon_positions[tid].items():
		exon_start = data[0]
		exon_end = data[1]
		chrom_exon_start = data[2]
		if (start >= exon_start) and (start <= exon_end) and (end >= exon_start) and (end <= exon_end):
			start_chrom_pos = chrom_exon_start + (start - exon_start)
			end_chrom_pos = chrom_exon_start + (end - exon_start)
			return start_chrom_pos, end_chrom_pos, eid
	return None, None, None

def compute_chrom_pos(df, file, exon_positions):#old
	v1s, v2s, v3s = [], [], []
	i = 0
	for _, row in df.iterrows():
		v1, v2, v3 = get_chrom_pos(row.bs_start, row.bs_end, row.mRNA, exon_positions)
		if i % 5000 == 0:
			gc.collect()
		v1s.append(v1)
		v2s.append(v2)
		v3s.append(v3)
		i += 1
	print('Done computation')
	df_result = pd.DataFrame({'chrom_bs_start': v1s,
							'chrom_bs_end': v2s,
							'exon_id': v3s})
	df = pd.concat([df,df_result], axis=1)
	df.to_parquet(file)
	gc.collect()
#	df.head()
#	return df

def calc_chrom_pos(start, end, exon_positions):
	for eid, data in exon_positions.items():
		exon_start = data[0]
		exon_end = data[1]
		chrom_exon_start = data[2]
		if (start >= exon_start) and (start <= exon_end) and (end >= exon_start) and (end <= exon_end):
			start_chrom_pos = chrom_exon_start + (start - exon_start)
			end_chrom_pos = chrom_exon_start + (end - exon_start)
			return start_chrom_pos, end_chrom_pos, eid
	return None, None, None

def chrom_pos(df, file, exon_positions):
	v1s, v2s, v3s = [], [], []
	i = 0
	for _, row in df.iterrows():
		v1, v2, v3 = calc_chrom_pos(row.bs_start, row.bs_end, exon_positions[row.mRNA])
		if i % 5000 == 0:
			gc.collect()
		v1s.append(v1)
		v2s.append(v2)
		v3s.append(v3)
		i += 1
	print('Done computation')
	df_result = pd.DataFrame({'chrom_bs_start': v1s,
							'chrom_bs_end': v2s,
							'exon_id': v3s})
	df.reset_index(drop=True, inplace=True)
	df_result.reset_index(drop=True, inplace=True)
	df = pd.concat([df,df_result], axis=1)
	df.to_parquet(file)
	gc.collect()