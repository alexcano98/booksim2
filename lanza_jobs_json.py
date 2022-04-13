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

#check that theres an argument
if len(sys.argv) < 2:
    print("Error: no json file specified")
    sys.exit(1)
json_file = sys.argv[1]


print(open(json_file).read())
data = json.load(open(json_file))

# read the fields topology, traffic, routing_function, num_vcs, injection rate and config file
traffic = data["traffic"]
routing_function = data["routing_function"]
num_vcs = data["num_vcs"]
injection_rate = data["injection_rate"]
config_file = data["config_file"]
topology = data["topology"]
results_dir = data["results_dir"]    
string_to_append = data["string"]



#check that exists the allocator field in the json file or if its empty
if "allocator" in data:
    allocator = str(data["allocator"])

    #if its a list return an error
    if isinstance(allocator, list):
        print("Error: allocator field is a list and not a string")
        sys.exit(1)

    #check if the string is not empty
    if allocator != "":
        string_to_append += " sw_allocator=" + allocator + " vc_allocator=" + allocator


#check that exists the packet size in the json file and if not, set it to default
if "packet_size" in data:
    packet_size = str(data["packet_size"])
    #check that packet size is not a list
    if isinstance(packet_size, list):
        print("Error: packet_size field is a list and not a string")
        sys.exit(1)
    string_to_append += " packet_size=" + packet_size


#check that theres a second argument and that is equal to 1 to remove the whole results directory
if len(sys.argv) == 3 and sys.argv[2] == "1":
    print("Removing previous results")
    shutil.rmtree(results_dir)

sbatch_string = "#!/bin/bash\n\
#SBATCH --job-name=booksim\n\
#SBATCH --cpus-per-task=1\n\
#SBATCH --ntasks=1\n\
#SBATCH --time=UNLIMITED"

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


#results_dir + topology variable
topology_dir = results_dir + "/" + topology

#check that the topology exists inside reulsts_dir
if not os.path.isdir(topology_dir):
    print("Topology directory not found")
    #create the directory
    os.makedirs(topology_dir)
    print("Topology directory created")

#save the config file inside the topology directory
shutil.copy(config_file, topology_dir)

#save the injection rates inside a file in the topology directory if it doesnt exist
if not os.path.isfile(topology_dir + "/injection_rates.txt"):
    with open(topology_dir + "/injection_rates.txt", "w") as f:
        for i in injection_rate:
            f.write(str(i) + " ")

#save the routing_functions inside a file  in the topology directory if it doesnt exist
if not os.path.isfile(topology_dir + "/routing_functions.txt"):
    with open(topology_dir + "/routing_functions.txt", "w") as f:
        for i in routing_function:
            f.write(str(i) + " ")
else: #if the file exists, add the new routing functions to the file
    with open(topology_dir + "/routing_functions.txt", "a") as f:
        for i in routing_function:
            f.write(str(i) + " ")

#save the num_vcs inside a file  in the topology directory if it doesnt exist
if not os.path.isfile(topology_dir + "/num_vcs.txt"):
    with open(topology_dir + "/num_vcs.txt", "w") as f:
        for i in num_vcs:
            f.write(str(i) + " ")

else: #if the file already exists, add the new num_vcs to the file
    with open(topology_dir + "/num_vcs.txt", "a") as f:
        for i in num_vcs:
            f.write(str(i) + " ")

#save the traffics inside a file  in the topology directory if it doesnt exist
if not os.path.isfile(topology_dir + "/traffics.txt"):
    with open(topology_dir + "/traffics.txt", "w") as f:
        for i in traffic:
            f.write(str(i) + " ")
else: #add the new traffic to the file
    with open(topology_dir + "/traffics.txt", "a") as f:
        for i in traffic:
            f.write(str(i) + " ");

#for each traffic, for each topology, for each injection rate, for each routing function, launch the job
for traffic in data["traffic"]:
    #create the directory for the traffic if does not exist
    traffic_dir = topology_dir + "/" + traffic
    if not os.path.isdir(traffic_dir):
        os.makedirs(traffic_dir)

    for num_vcs in data["num_vcs"]:
        for injection_rate in data["injection_rate"]:
            for routing_function in data["routing_function"]:
                #create a file to launch the job
                file_name = "sbatch_" + topology + "_" + str(num_vcs) + "_" + traffic + "_" + str(injection_rate) + "_" + routing_function + ".sbatch"
                file = open(file_name, "w")
                file.write(sbatch_string)
                file.write("\n")
                file.write("#SBATCH --output=" + topology_dir + "/" + traffic + "/" + str(num_vcs) + "_" + str(injection_rate) + "_" + routing_function +  ".out")
                
                file.write("\n")
                file.write("#SBATCH --error=" + topology_dir + "/" + traffic + "/" + str(num_vcs) + "_" + str(injection_rate) + "_" + routing_function + ".err")
                file.write("\n")
                file.write("./src/booksim " + config_file + " traffic=" + traffic + " num_vcs=" + str(num_vcs) + " injection_rate=" + str(injection_rate) + " routing_function=" + routing_function + string_to_append)
                file.close()
                #launch the job
                subprocess.call(["sbatch", file_name])
                #move the file to sbatch_archive dir and override it if it exists
                shutil.move(file_name, "sbatch_archive/" + file_name)
                print("Launched job for traffic: " + traffic + " topology: " + topology + " num_vcs: " + str(num_vcs) + " injection_rate: " + str(injection_rate) + " routing_function: " + routing_function)