import pandas as pd 
from pathlib import Path
import os
import argparse

def bp_file_in_path(path):
    for file in os.listdir(path):
        if file.endswith('.bp'): 
            return True
    return False

def check_args(input_path, output_path, log_path):
    
    if not input_path or not os.path.exists(input_path) or not bp_file_in_path(input_path):
        print(f"Error: You need to pass a path to the folder containing .bp files as argument.")
        exit()

    if log_path:
        logging = True
        log_path = args.log
        if not os.path.exists(log_path):
            os.makedirs(log_path)
    else:
        logging = False
        log_path = ""
        
    return input_path, output_path, log_path, logging

# read in binding site files (all files with .bs ending in path) produced by TarPmiR
# returns pandas DataFrame
def parse_tarpmir_bs(input_path, output_path):
    if not output_path: print(f"Bindingsites won't be written to disk as no output path was specified."")
    seperate_bs = []
    for filename in os.listdir(input_path):
        if filename.endswith('.bp'):
            file = input_path/filename
            if (open(file).readline()[0:5] != 'miRNA'):
                append_header(file)
            binding_sites = pd.read_csv(file, delimiter='	')
            # 1 column binding_sites to 2 new columns bs_start, bs_end
            new = binding_sites['binding_site'].str.split(',', 1, expand=True)
            binding_sites['bs_start'] = new[0].astype(int)
            binding_sites['bs_end'] = new[1].astype(int)
            binding_sites.drop(columns =["binding_site"], inplace = True)
            seperate_bs.append(binding_sites)
    all_bs = pd.concat(seperate_bs, axis=0, ignore_index=True) #14 388 + 10 326 = 24 714 rows
    print('All predicted bindingsites were read into pandas.')
    if output_path: 
        all_bs.to_parquet(output_path/'bs_all.parquet')
        print(f'All bindingsites were written into bs_all.parquet in {output_path}.')
    return all_bs

def main(input_path, output_path=False, log_path=""):
    input_path, output_path, log_path, logging = check_args(input_path, output_path, log_path)
    binding_sites = parse_tarpmir_bs(input_path, output_path)
    return binding_sites
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", help="path to the input data folder containing TarPmiR output files")
    parser.add_argument("-o", "--output", help="path to save the bs_all.parquet into")
    parser.add_argument("-l", "--log", help="save log files to folder")
    args = parser.parse_args()
    main(args.input, args.output, args.log)