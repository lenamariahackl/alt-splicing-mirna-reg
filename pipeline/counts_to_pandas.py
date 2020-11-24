import pandas as pd 
from pathlib import Path
import os
import argparse

def check_args(input_path, bindingsites, log_path):
    
    if not bindingsites:
        if not input_path or not os.path.exists(input_path):
            print(f"Error: You need to pass either a DataFrame bindingsites or a folder containing bs_subset.parquet.")
            exit()
        parquet_file = Path(input_path/'bs_subset.parquet') #TODO
        if parquet_file.is_file():
            print(f"No bindingsites DataFrame was passed as argument. Defaulted back to using bs_subset.parquet from folder {input_path}.")
            bindingsites = pd.read_parquet(parquet_file)
        else:
            print(f"Error: No DataFrame bindingsites was passed as argument and bs_subset.parquet couldn't be found in folder {input_path}. Please first execute the prior steps.")
            exit()

    if log_path:
        logging = True
        log_path = args.log
        if not os.path.exists(log_path):
            os.makedirs(log_path)
    else:
        logging = False
        log_path = ""
        
    return input_path, bindingsites, log_path, logging


#TODO add logging
def main(input_path="", bindingsites=False, log_path=""):
    input_path, bindingsites, log_path, logging = check_args(input_path, bindingsites, log_path)
    
        
    
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", help="path to the input data folder")
    parser.add_argument("-b", "--bindingsites", help="load calculated bindingsites DataFrame from parquet file")
    parser.add_argument("-l", "--log", help="save log files to folder")
    args = parser.parse_args()
    main(args.input, args.bindingsites, args.log)