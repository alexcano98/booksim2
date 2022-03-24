#!/bin/bash
        #SBATCH --job-name=booksim
        #SBATCH -D .
        #SBATCH --output= results/flat_fly/flatfly_topology_2vcs_sim/uniform/sim_0.20_valiant.out
        #SBATCH --error= results/flat_fly/flatfly_topology_2vcs_sim/uniform/sim_0.20_valiant.err
        #SBATCH --cpus-per-task=1
        #SBATCH --ntasks=1
        #SBATCH --time=1-11:59:59
        #SBATCH --mem=5GB
 /home/alejandro/booksim2/src/booksim injection_rate=0.20
         traffic=uniform routing_function=valiant
