#!/bin/bash

# By Alejandro Cano Cos and Github Copilot
# This script tails the output file of the jobs RUNNING of the user

# 12/04/2022 first version

state="RUNNING"
user="$USER"

#save first argumentn as the name of the job, if no argument default to "all"
if [ -z "$1" ]
then
    job_name="*"
else
    job_name="$1"
fi


#get the job ids of the jobs with the name $job_name

#job_ids=$(squeue -o "%i" -h -u $user -t $state)

job_ids=$(squeue -o "%i %j" -h -u $user -t $state)

#get the first column of all the lines with the name $job_name
job_ids=$(echo "$job_ids" | egrep "$job_name" | cut -d " " -f 1)

#get the details of the jobs
for job_id in $job_ids; do
    echo "Job ID: $job_id"
    #get the output directory of the jobs 
    scontrol show jobid -dd $job_id | grep StdOut | cut -d "=" -f 2

    scontrol show jobid -dd $job_id | grep StdOut | cut -d "=" -f 2 | xargs tail

    echo "=============================="

done
