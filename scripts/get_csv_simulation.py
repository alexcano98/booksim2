#script that inside a directory, get the csv files and merge them into one csv file
import os
import pandas as pd
import numpy as np
import glob
import sys
import argparse
import re

#get the arguments
parser = argparse.ArgumentParser()
parser.add_argument('-d', '--directory', help='directory where the csv files are', required=True)
parser.add_argument('-o', '--output', help='output file name', required=True)
args = parser.parse_args()
#get the csv header with ; as separator

csv_header = pd.read_csv("scripts/params/header-csv.csv", sep=";")
print(csv_header)

#read the trafic files in traffic file as plain text
traffics = []
# Check if traffics exists
if not os.path.exists(args.directory + "traffics.txt"):
    print(args.directory + "traffics.txt")
    print("Error: traffics file does not exist")
    exit(1)
# Read traffics from file as a float separated by space
with open(args.directory + "traffics.txt", "r") as f:
    traffics = [x for x in f.read().split()]
print("traffics: ", traffics)

for t in traffics:
    #iter the files in the directory recived as argument and ends with _sim_out.csv
    files = glob.glob(args.directory + t +"/*_sim_out.csv")

    print("files in " + args.directory + t +'/*_sim_out.csv' + " files: " + str(files))

    for file in files:
        #read the file
        df = pd.read_csv(file)

        #split the file name by "_"
        file_name = file.split("_")

        vcs = file_name[1]
        routing = file_name[3]
        #print(str(file_name))
        #append the vcs and routing to the csv file, two lasts columns
        df['vcs'] = vcs
        df['routing'] = routing
        print(df)
        df.to_csv(args.output, mode='a', header=False, index=False) # mode a means that the file will be appended