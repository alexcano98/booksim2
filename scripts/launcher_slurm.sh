#!/bin/bash

###################################################################
#Script Name	:   launcher.sh
#Description	:   Launcher for Booksim2 sims in a Slurm cluster
#Args           :   <CONFIG_FILE> <PARAMS> <PURGE>
#Author       	:   Daniel Postigo Diaz
#Email         	:   daniel.postigo@alumnos.unican.es
###################################################################

LC_NUMERIC="en_US.UTF-8"
# Booksim binary
BOOKSIM_HOME=/home/alejandro/booksim2/src
# Config file for the simulation
CONFIG_FILE=$1
# Injection rates
PARAMS=$2

PURGE=0
PURGE="${3:-0}"
echo "${PURGE}"

# Check the number of arguments
if [ $# -ne 3 ]; then
    echo "Usage: $0 <CONFIG_FILE> <PARAMS> <PURGE>"
    exit 1
fi

echo "Config file: $CONFIG_FILE"

# Check if the config file exists
if ! [ -f "$CONFIG_FILE" ]; then
    echo "Error: $CONFIG_FILE does not exist"
    exit 1
fi


TO_REMOVE=runfiles\/
OUT_DIR=results/${CONFIG_FILE##$TO_REMOVE}_sim
echo "OUT DIR: $OUT_DIR"

# Check if the OUT_DIR directory already exists
uno=1
if [ -d "$OUT_DIR" ] && [ "${PURGE}" -eq "1" ]; then
    echo "Removing: "$OUT_DIR" "
    rm -rf $OUT_DIR
    #exit 1
fi

mkdir -p ${OUT_DIR}
echo "hola"
# Copy the config file to simulation directory
cp ${CONFIG_FILE} ${OUT_DIR}/config
cp ${PARAMS}/* ${OUT_DIR}/

for traffic in `cat ${PARAMS}/traffics` #extremely bbuggy
do
  if ! [ -d "${OUT_DIR}/${traffic}" ]; then
    mkdir -p ${OUT_DIR}/${traffic}
    for inj_rate in `cat ${PARAMS}/injection_rates` #extremely bbuggy
    do

      for routing in `cat ${PARAMS}/routings` #extremely bbuggy
      do

        echo "#!/bin/bash
        #SBATCH --job-name=booksim
        #SBATCH -D .
        #SBATCH --output= ${OUT_DIR}/${traffic}/sim_${inj_rate}_${routing}.out
        #SBATCH --error= ${OUT_DIR}/${traffic}/sim_${inj_rate}_${routing}.err
        #SBATCH --cpus-per-task=1
        #SBATCH --ntasks=1
        #SBATCH --time=1-11:59:59
        #SBATCH --mem=5GB" > ${traffic}_sim_${inj_rate}_${routing}.sbatch


        # Create the sbatch body
        echo "/usr/bin/time --verbose -o ./time_results/sim_${inj_rate}.time ${BOOKSIM_HOME}/booksim injection_rate=${inj_rate}
         traffic=${traffic} routing_function=${routing}" >> ${traffic}_sim_${inj_rate}_${routing}.sbatch

        # Submit the sbatch
        echo "Submitting sbatch ${traffic}_sim_${inj_rate}_${routing}.sbatch"
        sbatch ${traffic}_sim_${inj_rate}_${routing}.sbatch

        # Archive the sbatch
        mv ${traffic}_sim_${inj_rate}_${routing}.sbatch ./sbatch_archive/

        #echo $inj_rate
        ./src/booksim ${CONFIG_FILE} injection_rate="${inj_rate}" traffic="${traffic}" routing_function="${routing}" >> ${OUT_DIR}/${traffic}/sim_${inj_rate}_${routing}.out
        echo "Done: ${inj_rate}, ${routing}, ${traffic}"

      done
    done
  fi
done
