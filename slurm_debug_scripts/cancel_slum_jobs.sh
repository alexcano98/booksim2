#!/bin/bash

# By Alejandro Cano Cos and Github Copilot
# This script tails the output file of the jobs RUNNING of the user

# 12/04/2022 first version

user="$USER"

#save first argumentn as the name of the job, if no argument default to "all"
if [ -z "$1" ]
then
    # exit the script if no argument is given
    echo "No argument given"
    exit 1
else
    job_name="$1"
fi

#ask the user if he really wanna cancel the jobs with the name $job_name
echo "Are you sure you wanna cancel the jobs with the name $job_name? (y/n)"
read answer

# if the user answers "y" or "Y" then cancel the jobs
if [ "$answer" != "y" ] && [ "$answer" != "Y" ]
then
    # if the user answers "n" or "N" then exit the script
    echo "Exiting the script"
    exit 1
fi

job_ids=$(squeue -o "%i %j" -h -u $user)

#get the first column of all the lines with the name $job_name
job_ids=$(echo "$job_ids" | egrep "$job_name" | cut -d " " -f 1)

#get the details of the jobs
for job_id in $job_ids; do
    echo "Cancelling job with ID: $job_id"
    scancel $job_id
    echo "=============================="

done
