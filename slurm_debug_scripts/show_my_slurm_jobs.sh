#!/bin/bash

# By Alejandro Cano Cos and Github Copilot
# This script shows the details of the jobs of the user, by default the selected state is running

# First argument is the state of the job, by default is running
# Second argument is the user, by default is the user executing the script

# 12/04/2022 first version


#get the first argument and set it as the state
state=$1

#get the second argument and set it as the user
user=$2

#if state is not set, set it to RUNNING
if [ -z "$state" ]; then
    state="RUNNING"
fi

#if user is not set, set it to $USER
if [ -z "$user" ]; then
    user="$USER"
fi

#get the job ids of the jobs in the state
job_ids=$(squeue -o "%i" -h -u $user -t $state)

#get the details of the jobs
for job_id in $job_ids; do
    echo "Job ID: $job_id"
    scontrol show jobid -dd $job_id
    echo "=============================="
done
