#script that reads the path of a json file as argument and launches the jobs
#read the json file
import json
import os
import sys
import subprocess
import time
import datetime
import shutil
import glob
import re
import argparse
import logging
import logging.handlers
import traceback

json = sys.argv[1]

data = json.load(open(json))

# read the fields topology, traffic, routing_function, num_vcs, injection rate and config file
traffic = data["traffic"]
routing_function = data["routing_function"]
num_vcs = data["num_vcs"]
injection_rate = data["injection_rate"]
config_file = data["config_file"]
topology = data["topology"]
results_dir = data["results_dir"]

sbatch_string = "\#!/bin/bash
        \#SBATCH --job-name=booksim
        \#SBATCH --cpus-per-task=1
        \#SBATCH --ntasks=1
        \#SBATCH --time=UNLIMITED;"

#check that the config file exists
if not os.path.isfile(config_file):
    print("Config file not found")
    sys.exit(1)

#check that results_dir exists
if not os.path.isdir(results_dir):
    print("Results directory not found")
    #create the directory
    os.makedirs(results_dir)
    print("Results directory created")



#for each traffic, for each topology, for each injection rate, for each routing function, launch the job
for traffic in data["traffic"]:
    for num_vcs in data["num_vcs"]:
        for injection_rate in data["injection_rate"]:
            for routing_function in data["routing_function"]:
                #create a file to launch the job
                file_name = "sbatch_" + topology + "_" + num_vcs + "_" + traffic + "_" + injection_rate + "_" + routing_function + ".sbatch"
                file = open(file_name, "w")
                file.write(sbatch_string)
                file.write("\n")
                file.write("\#SBATCH --output=logs/%j.out")
                file.write("\n")
                file.write("\#SBATCH --error=logs/%j.err")
                file.write("\n")
                file.write("./src/booksim" + config_file + " traffic=" + traffic + " num_vcs=" + num_vcs + " injection_rate=" + injection_rate + " routing_function=" + routing_function)
                file.close()
                #launch the job
                subprocess.call(["sbatch", file_name])
                #move the file to sbatch_archive dir
                shutil.move(file_name, "sbatch_archive")
                print("Launched job for traffic: " + traffic + " topology: " + topology + " num_vcs: " + num_vcs + " injection_rate: " + injection_rate + " routing_function: " + routing_function)



               

