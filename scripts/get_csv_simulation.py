#script that inside a directory, get the csv files and merge them into one csv file
import os
import pandas as pd
import numpy as np
import glob
import sys
import argparse

#get the arguments
parser = argparse.ArgumentParser()
parser.add_argument('-d', '--directory', help='directory where the csv files are', required=True)
parser.add_argument('-o', '--output', help='output file name', required=True)
args = parser.parse_args()
#get the csv header with ; as separator

csv_header = pd.read_csv("scripts/params/header-csv.csv", sep=";")
#add vcs and routing fields to csv_header
csv_header = csv_header.append(pd.DataFrame([['vcs', 'vcs'], ['routing', 'routing']], columns=csv_header.columns))
print(csv_header)

"""
#iter the files in the directory recived as argument and ends with *_sim_out.csv
files = glob.glob(args.directory + '/*_sim_out.csv')
#iter 
for file in files:
    #read the file
    df = pd.read_csv(file)
    #split the file name by "_"
    file_name = file.split("_")
    vcs = file_name[0]
    routing = file_name[2]
    #append the vcs and routing to the csv file, two lasts columns
    df['vcs'] = vcs
    df['routing'] = routing
    df.to_csv(args.output, mode='a', header=False, index=False)
"""
