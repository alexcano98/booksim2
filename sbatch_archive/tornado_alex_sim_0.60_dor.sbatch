#!/bin/bash
        #SBATCH --job-name=booksim
        #SBATCH -D .
        #SBATCH --output= results/flat_fly/flatfly_topology_2vcs_sim/tornado_alex/sim_0.60_dor.out
        #SBATCH --error= results/flat_fly/flatfly_topology_2vcs_sim/tornado_alex/sim_0.60_dor.err
        #SBATCH --cpus-per-task=1
        #SBATCH --ntasks=1
        #SBATCH --time=1-11:59:59
        #SBATCH --mem=5GB
 /home/alejandro/booksim2/src/booksim injection_rate=0.60
         traffic=tornado_alex routing_function=dor
