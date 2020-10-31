import re
import pandas as pd 
from pathlib import Path
import pickle

# returns format {gene_id : {transcript_id : exon_id}} 
def read_in_ids(file):
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
				transcript_id = re.search('ENST[^_]*',line).group(0)
				exon_id = file[id_+1][:-1]
				dict_[transcript_id] = dict_.get(transcript_id, [])
				dict_[transcript_id].append(exon_id)
		#write to disk
		f = open(pickle1,"wb")
		pickle.dump(dict_,f)
		f.close()
	return dict_

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
	pickle0 = Path(path/"chrom_exon_starts.pkl")
	pickle1 = Path(path/"exon_starts.pkl")
	pickle2 = Path(path/"exon_ends.pkl")
	feather_file = Path(path/'exon_info.feather')
	if pickle0.is_file() and pickle1.is_file() and pickle2.is_file() and feather_file.is_file():
		#load files from disk
		pkl_file0 = open(pickle0, 'rb')
		chrom_exon_starts = pickle.load(pkl_file0)
		pkl_file0.close()
		pkl_file1 = open(pickle1, 'rb')
		exon_starts = pickle.load(pkl_file1)
		pkl_file1.close()
		pkl_file2 = open(pickle2, 'rb')
		exon_ends = pickle.load(pkl_file2)
		pkl_file2.close()
		exon_info = pd.read_feather(feather_file)
	else:
		chrom_exon_starts = read_in_fasta(path/'all_exon_start.fasta')
		chrom_exon_ends = read_in_fasta(path/'all_exon_end.fasta')
		exon_starts = {}
		exon_ends = {}
		exon_info = {}
		for tid in chrom_exon_starts:
			last_exon_end = 0
			sorted_ = {k: v for k, v in sorted(chrom_exon_starts[tid].items(), key=lambda item: item[1])}
			for eid in sorted_:
				chrom_exon_start = chrom_exon_starts[tid][eid]
				chrom_exon_end = chrom_exon_ends[tid][eid]
				exon_info[eid] = [tid, chrom_exon_start, chrom_exon_end]
				if chrom_exon_end - chrom_exon_start < 0:
					print('Error: difference is negative! start of',tid,'is bigger than end!')
				exon_starts[tid] = exon_starts.get(tid,{})
				exon_ends[tid] = exon_ends.get(tid,{})
				exon_start = last_exon_end + 1
				exon_end = exon_start + (chrom_exon_end - chrom_exon_start)
				last_exon_end = exon_end
				exon_starts[tid][eid] = exon_start
				exon_ends[tid][eid] = exon_end
				#write to disk
		exon_info = pd.DataFrame.from_dict(exon_info, orient='index', columns=['transcript_id','chrom_exon_start', 'chrom_exon_end']).reset_index()
		exon_info = exon_info.rename(columns={'index': "exon_id"})
		f = open(pickle0,"wb")
		pickle.dump(chrom_exon_starts,f)
		f.close()
		f = open(pickle1,"wb")
		pickle.dump(exon_starts,f)
		f.close()
		f = open(pickle2,"wb")
		pickle.dump(exon_ends,f)
		f.close()
		exon_info.to_feather(feather_file)
	return chrom_exon_starts, exon_starts, exon_ends, exon_info