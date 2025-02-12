import csv
import sys
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os

P_LAT = 5
P_W_LAT = 6
N_LAT = 8
N_W_LAT = 9
F_LAT = 11
FRAGM = 14
P_INJ = 17
P_ACC = 20
F_INJ = 23
F_ACC = 26
F_W_ACC = 25
H_AVG = 30


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
    vc_utilization = "vc_util.csv"
    fairness_injected = "jain_injected.csv"
    fairness_accepted = "jain_accepted.csv"

    topology  = sys.argv[2]

    if topology == "":
        print("Topology not specified")
        exit(1)

    # check if the output directory exists
    if not os.path.exists(out_dir):
        print("Error: Output directory does not exist")
        exit(1)
    # Delete / from the end of the directory
    if out_dir[-1] == '/':
        out_dir = out_dir[:-1]


    # Check if injection_rate exists
    if not os.path.exists(out_dir + "/injection_rates.txt"):
        print("Error: injection_rates file does not exist")
        exit(1)
    # Read injection rates from file as a float separated by space
    with open(out_dir + "/injection_rates.txt", "r") as f:
        injected_rate = [float(x) for x in f.read().split()]
    print("Injection rates: ", injected_rate)


    # Check if routings exists
    if not os.path.exists(out_dir + "/routing_functions.txt"):
        print("Error: routing_functions.txt file does not exist")
        exit(1)
    # Read routings from file as a float separated by space
    with open(out_dir + "/routing_functions.txt", "r") as f:
        routings = [x for x in f.read().split()]
    print("routing_functions.txt ", routings)


    # Check if traffics exists
    if not os.path.exists(out_dir + "/traffics.txt"):
        print("Error: traffics file does not exist")
        exit(1)
    # Read traffics from file as a float separated by space
    with open(out_dir + "/traffics.txt", "r") as f:
        traffics = [x for x in f.read().split()]
    print("traffics: ", traffics)

    #check if num_vcs exists
    if not os.path.exists(out_dir + "/num_vcs.txt"):
        print("Error: num_vcs file does not exist")
        exit(1)
    # Read num_vcs from file as a float separated by space
    with open(out_dir + "/num_vcs.txt", "r") as f:
        num_vcs = [x for x in f.read().split()]
    print("num_vcs: ", num_vcs)


    #check if allocators exists
    if not os.path.exists(out_dir + "/allocators.txt"):
        print("Error: allocators file does not exist")
        exit(1)
    # Read allocators from file as a float separated by space
    with open(out_dir + "/allocators.txt", "r") as f:
        allocators = [x for x in f.read().split()]

    
    #check if packet_sizes exists
    if not os.path.exists(out_dir + "/packet_sizes.txt"):
        print("Error: packet_sizes file does not exist")
        exit(1)
    # Read packet_sizes from file as a float separated by space
    with open(out_dir + "/packet_sizes.txt", "r") as f:
        packet_sizes = [x for x in f.read().split()]


    flits_latency_list = {} #'injected_rate': injected_rate
    accepted_flits_list = {}
    traffic_pattern_list = {}
    hops_avg_list = {}
    worst_flits_latency_list = {}
    worst_accepted_flits_list = {}
    vc_utilization_values_list = {}

    fairness_injected_values_list = {}
    fairness_accepted_values_list = {}

    os.chdir(out_dir)
    for t in traffics:
        # move to output directory
        os.chdir(t)
        for a in allocators:
            for n in num_vcs:
                for p in packet_sizes:
                    sim_file_plot = n +"_" + a + "_" + "_"+ p
                    error = 0
                    for routing in routings:
                        # check if csv files exist

                        sim_file= n +"_" + a + "_" + routing  + "_"+ p

                        
                        if not os.path.exists(sim_file + "_" + sim_out):
                            print("sim_file: \n", sim_file)
                            print("current directory: \n", os.getcwd())
                            print("Error: sim_out.csv does not exist")
                            exit(1)

                        if not os.path.exists(sim_file + "_" + cycles):
                            print("sim_file: \n", sim_file)
                            print("current directory: \n", os.getcwd())
                            print("Error: cycles.csv does not exist")
                            exit(1)

                        if not os.path.exists(sim_file + "_" + run_time):
                            print("sim_file: \n", sim_file)
                            print("current directory: \n", os.getcwd())
                            print("Error: run_time.csv does not exist")
                            exit(1)

                        # Check if vc_utilization exists
                        if not os.path.exists(sim_file + "_" + vc_utilization):
                            print("sim_file: \n", sim_file)
                            print("current directory: \n", os.getcwd())
                            print("Error: vc_util.csv does not exist")
                            exit(1)

                        # Check if fairness_injected exists
                        if not os.path.exists(sim_file + "_" + fairness_injected):
                            print("sim_file: \n", sim_file)
                            print("current directory: \n", os.getcwd())
                            print("Error: jain_injected.csv does not exist")
                            exit(1)
                        
                        # Check if fairness_accepted exists
                        if not os.path.exists(sim_file + "_" + fairness_accepted):
                            print("sim_file: \n", sim_file)
                            print("current directory: \n", os.getcwd())
                            print("Error: jain_accepted.csv does not exist")
                            exit(1)

                        # create plots directory
                        print("Creating plots directory...")
                        if not os.path.exists('plots'):
                            os.makedirs('plots')

                        #catch any exception and throw it again
                        try:
                            df_sim = pd.read_csv(sim_file + "_" +sim_out, sep=',', header=None)
                            #read df_vc_util file that contains in each line the vc utilization
                            df_vc_util = pd.read_csv(sim_file + "_" +vc_utilization, sep='\t', header=None)
                            df_fairness_injected = pd.read_csv(sim_file + "_" + fairness_injected, sep='\t', header=None)
                            df_fairness_accepted = pd.read_csv(sim_file + "_" + fairness_accepted, sep='\t', header=None)
                            #print("df_vc_util: \n", df_vc_util)


                        except Exception as e:
                            print("==============================ERROR===============================")
                            #Show filename and current directory
                            print("sim_file: \n", sim_file)
                            print("current directory: \n", os.getcwd())

                            #Show exception
                            print("Exception: \n", e)
                            break

                        accepted_flits = df_sim[F_ACC]
                        hops_avg = df_sim[H_AVG]
                        traffic_pattern = df_sim[1][0]
                        worst_accepted_flits = df_sim[F_W_ACC]
                        # get df_vc_util values of the first column in a list
                        vc_utilization_values = df_vc_util[0].tolist()
                        # get df_fairness_injected values of the first column in a list
                        fairness_injected_values = df_fairness_injected[0].tolist()
                        # get df_fairness_accepted values of the first column in a list
                        fairness_accepted_values = df_fairness_accepted[0].tolist()


                        accepted_flits_list[routing] = round(100 * accepted_flits.mean())
                        #print("Accepted flits: ", accepted_flits.mean())
                        traffic_pattern_list[routing] = traffic_pattern
                        hops_avg_list[routing] = float(hops_avg)
                        worst_accepted_flits_list[routing] = worst_accepted_flits
                        vc_utilization_values_list[routing] = vc_utilization_values

                        fairness_injected_values_list[routing] = fairness_injected_values[0]
                        fairness_accepted_values_list[routing] = fairness_accepted_values[0]


                    
                    markers = ['v', '^', '*', '<', '>',  '*', 'p', '*', '*', '*', '*', '*']
                    if error == 1:
                        continue
                    #Throughput
                    #print(accepted_flits_list)
                    d = {"routings": list(accepted_flits_list.keys()), "accepted_flits": list(accepted_flits_list.values())}

                    data = pd.DataFrame(data=d) 
                    
                    ax =  data.plot.bar(y="accepted_flits", x="routings", legend=False, title='Throughput')
                    ax.set_xlabel('Routing Functions')
                    ax.set_ylabel('Accepted Flits')
                    ax.set_title('Accepted Flits')
                    ax.set_xticklabels(routings, rotation=30)
                    #show numbers on bars
                    for p in ax.patches:
                        ax.annotate(str(p.get_height()), (p.get_x() * 1.005, p.get_height() * 1.005))
                    #show the grid
                    ax.grid(True)
                    plt.title('Throughput [TP={}]'.format(traffic_pattern))
                    plt.tight_layout()

                    plt.savefig('./plots/'+topology+"_"+sim_file_plot+'_throughput_'+t+'.png')
                    print("Throughput plot created")

                    #HOPS
                    d = {"routings": list(hops_avg_list.keys()), "hops_avg": list(hops_avg_list.values())}
                    data = pd.DataFrame(data=d)
                    ax =  data.plot.bar(y="hops_avg", x="routings", legend=False, title='Hops')
                    ax.set_xlabel('Routing Functions')
                    ax.set_ylabel('Hops')
                    ax.set_title('Hops')
                    ax.set_xticklabels(routings, rotation=30)
                    #show numbers on bars
                    for p in ax.patches:
                        ax.annotate(str(p.get_height()), (p.get_x() * 1.005, p.get_height() * 1.005))
                    #show the grid
                    ax.grid(True)
                    plt.tight_layout()

                    plt.title('Hops [TP={}]'.format(traffic_pattern))
                    plt.savefig('./plots/'+topology+"_"+sim_file_plot+'_hops_'+t+'.png')
                    print("Hops plot created")


                    #VC Utilization
                    data = pd.DataFrame(vc_utilization_values_list)
                    data.plot(y=routings, legend=True, title='% VC Utilization')

                    ax = data.plot()
                    for i, line in enumerate(ax.get_lines()):
                        if line.get_label() == 'injected_rate':
                            print("pasando")
                            continue
                        line.set_marker(markers[i])
                        #marker also the legend

                    ax.legend(loc='upper left')

                     #show only x range between 0 and 1.1
                    #plt.xlim(0, num_vcs)
                    plt.ylabel('% VC Utilization')
                    plt.xlabel('VCs')
                    plt.title('VCs util [TP={}]'.format(traffic_pattern))

                    #show grid 0.1
                    plt.grid(True, which='both', alpha=0.1)
                    plt.tight_layout()
                    #plt.show()
                    plt.savefig('./plots/'+topology+"_"+sim_file_plot+'_VC_util_'+t+'.png')

                    #Fairness injected
                    d = {"routings": list(fairness_injected_values_list.keys()), "fairness_injected": fairness_injected_values_list.values()}
                    #print(d)
                    data = pd.DataFrame(data=d)
                    #print(data)
                    ax =  data.plot.bar(x="routings", y="fairness_injected", legend=False, title='Jain index injected')
                    ax.set_xlabel('Routing Functions')
                    ax.set_ylabel('Index')
                    ax.set_title('Jain index injection')
                    #show numbers on bars
                    for p in ax.patches:
                        ax.annotate(str(p.get_height()), (p.get_x() * 1.005, p.get_height() * 1.005))
                    #show the grid
                    ax.grid(True)

                    plt.tight_layout()
                    plt.title('Fairness [TP={}]'.format(traffic_pattern))
                    plt.savefig('./plots/'+topology+"_"+sim_file_plot+'_fairness_injected_'+t+'.png')
                    print("Fairness plot created")

                    #Fairness accepted
                    d = {"routings": list(fairness_accepted_values_list.keys()), "fairness_accepted": list(fairness_accepted_values_list.values())}
                    data = pd.DataFrame(data=d)
                    ax =  data.plot.bar(y="fairness_accepted", x="routings", legend=False, title='Jain index accepted')
                    ax.set_xlabel('Routing Functions')
                    ax.set_ylabel('Index')
                    ax.set_title('Jain index accepted')
                    ax.set_xticklabels(routings, rotation=30)
                    #show numbers on bars
                    for p in ax.patches:
                        ax.annotate(str(p.get_height()), (p.get_x() * 1.005, p.get_height() * 1.005))
                    #show the grid
                    ax.grid(True)
                    plt.title('Fairness [TP={}]'.format(traffic_pattern))
                    plt.tight_layout()
                    plt.savefig('./plots/'+topology+"_"+sim_file_plot+'_fairness_accepted_'+t+'.png')
                    print("Fairness plot created")





                    """
                    #PLOT 4
                    data = pd.DataFrame(worst_accepted_flits_list)
                    data.plot(y=routings, legend=True, title='Worst throughput')
                    
                    ax = data.plot(kind='line')
                    for i, line in enumerate(ax.get_lines()):
                        if line.get_label() == 'injected_rate':
                            print("pasando")
                            continue
                        line.set_marker(markers[i])
                        line.set_xdata(injected_rate)
                        #marker also the legend

                    ax.legend(loc='upper left')

                    #show only x range between 0 and 1.1
                    plt.xlim(0, 1.1)
                    plt.ylim(0, 1.1)

                    #show grid 0.1 separation between lines in the grid
                    
                    plt.grid(True, which='both', alpha=0.1)

                    plt.ylabel('Accepted flits')
                    plt.xlabel('Injected Rate (Flits/cycle/node)')
                    plt.title('Throughput [TP={}]'.format(traffic_pattern))

                    #plt.grid()
                    #plt.show()
                    plt.savefig('./plots/'+topology+"_"+sim_file_plot+'_worst_throughput_'+t+'.png')
                    print("Throughput plot created")

                    #PLOT 5
                    data = pd.DataFrame(hops_avg_list)
                    data.plot(y=routings, legend=True, title='Hops average')

                    ax = data.plot(kind='line')
                    for i, line in enumerate(ax.get_lines()):
                        if line.get_label() == 'injected_rate':
                            print("pasando")
                            continue
                        line.set_marker(markers[i])
                        line.set_xdata(injected_rate)
                        #marker also the legend

                    ax.legend(loc='upper left')

                    #show only x range between 0 and 1.1
                    plt.xlim(0, 1.1)
                    plt.ylabel('Avg hops')
                    plt.xlabel('Injected Rate (Flits/cycle/node)')
                    plt.title('Hops [TP={}]'.format(traffic_pattern))

                    #show grid 0.1
                    plt.grid(True, which='both', alpha=0.1)
                    #plt.show()
                    plt.savefig('./plots/'+topology+"_"+sim_file_plot+'_hops_'+t+'.png')
                    print("Hops plot created") """

        os.chdir("../") # volvemos atras 




        """
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
        """

    # Add usage message
if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python parser_2.py <output_directory>")
        exit(1)
    main()
