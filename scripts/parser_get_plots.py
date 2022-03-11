import csv
import sys
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os

P_LAT = 5
N_LAT = 8
F_LAT = 11
FRAGM = 14
P_INJ = 17
P_ACC = 20
F_INJ = 23
F_ACC = 26


def main():
    """
    Main function
    """
    flits_latency = []
    accepted_flits = []
    injected_rate = []

    routings = []
    traffics = []
#   injected_rate = [0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 1.0]

    out_dir = sys.argv[1]
    sim_out = "sim_out.csv"
    cycles = "cycles.csv"
    run_time = "run_time.csv"


    # check if the output directory exists
    if not os.path.exists(out_dir):
        print("Error: Output directory does not exist")
        exit(1)
    # Delete / from the end of the directory
    if out_dir[-1] == '/':
        out_dir = out_dir[:-1]


    # Check if injection_rate exists
    if not os.path.exists(out_dir + "/injection_rates"):
        print("Error: injection_rates file does not exist")
        exit(1)
    # Read injection rates from file as a float separated by space
    with open(out_dir + "/injection_rates", "r") as f:
        injected_rate = [float(x) for x in f.read().split()]
    print("Injection rates: ", injected_rate)


    # Check if routings exists
    if not os.path.exists(out_dir + "/routings"):
        print("Error: routings file does not exist")
        exit(1)
    # Read routings from file as a float separated by space
    with open(out_dir + "/routings", "r") as f:
        routings = [float(x) for x in f.read().split()]
    print("routings ", routings)


    # Check if traffics exists
    if not os.path.exists(out_dir + "/traffics"):
        print("Error: traffics file does not exist")
        exit(1)
    # Read traffics from file as a float separated by space
    with open(out_dir + "/traffics", "r") as f:
        traffics = [float(x) for x in f.read().split()]
    print("traffics: ", traffics)


    for t in traffics:

        # move to output directory
        os.chdir(out_dir+t)

        # check if csv files exist
        if not os.path.exists(sim_out):
            print("Error: sim_out.csv does not exist")
            exit(1)

        if not os.path.exists(cycles):
            print("Error: cycles.csv does not exist")
            exit(1)
        if not os.path.exists(run_time):
            print("Error: run_time.csv does not exist")
            exit(1)


        # create plots directory
        print("Creating plots directory...")
        if not os.path.exists('plots'):
            os.makedirs('plots')

        df_sim = pd.read_csv(sim_out, sep=',', header=None)
        flits_latency = df_sim[F_LAT]
        accepted_flits = df_sim[F_ACC]
        traffic_pattern = df_sim[1][0]

        
        x = pd.DataFrame({'injected_rate': injected_rate, 'flits_latency': flits_latency, 'accepted_flits': accepted_flits})
        x.plot(x='injected_rate', y='flits_latency', marker="o", kind='line', color='red', label='flits_latency', legend=True, title='Flits latency')
        plt.ylabel('Latency (cycles)')
        plt.xlabel('Injected rate (Flits/cycle/node)')
        plt.title('Flits latency (cycles) [TP={}]'.format(traffic_pattern))
        plt.xticks(np.arange(0, 1.1, 0.1))
        plt.grid()
        plt.savefig('./plots/flits_latency.png')
        print("Flits latency plot created")

        x.plot(x='injected_rate', y='accepted_flits', marker="o" ,kind='line', color='red', label='accepted_flits', legend=True, title='Accepted flits', ylim=(0, 1))
        plt.ylabel('Accepted flits')
        plt.xlabel('Injected Rate (Flits/cycle/node)')
        plt.title('Throughput [TP={}]'.format(traffic_pattern))
        plt.xticks(np.arange(0, 1.1, 0.1))
        plt.yticks(np.arange(0, 1.1, 0.1))
        plt.grid()
        plt.savefig('./plots/throughput.png')
        print("Throughput plot created")

        df_cycles = pd.read_csv(cycles, sep=',', header=None)
        x = pd.DataFrame({'injected_rate': injected_rate, 'cycles' : df_cycles[0]})
        x.plot(x='injected_rate', y='cycles', marker="o",kind='line', color='red', label='cycles', legend=True)
        plt.xlabel('Injected Rate (Flits/cycle/node)')
        plt.ylabel('Simulation Cycles (inc. WARMUP')
        plt.title('Simulation Cycles [TP={}]'.format(traffic_pattern))
        plt.xticks(np.arange(0, 1.1, 0.1))
        plt.grid()
        plt.savefig('./plots/cycles.png')
        print("Cycles plot created")

        df_run = pd.read_csv(run_time, sep=',', header=None)
        x = pd.DataFrame({'injected_rate': injected_rate, 'run time' : df_run[0]})
        x.plot(x='injected_rate', y='run time', marker="o", kind='line', color='red', label='run time', legend=True, title='Run time')
        plt.ylabel('Time (s)')
        plt.xlabel('Injected Rate (Flits/cycle/node)')
        plt.title('Run time [TP={}]'.format(traffic_pattern))
        plt.xticks(np.arange(0, 1.1, 0.1))
        plt.grid()
        plt.savefig('./plots/run_time.png')
        print("Run time plot created")

    # Add usage message
if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python parser_2.py <output_directory>")
        exit(1)
    main()
