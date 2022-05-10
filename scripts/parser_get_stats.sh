#!/bin/bash

###################################################################
#Script Name	:   parser_1.sh
#Description	:   Parser Booksim2 sims outputs
#Args           :   <OUT_DIR>
#Author       	:   Daniel Postigo Diaz
#Email         	:   daniel.postigo@alumnos.unican.es
###################################################################

OUT_DIR=$1

# Check the number of arguments
if [ $# -ne 1 ]; then
    echo "Usage: $0 <PARSED_DIR>"
    exit 1
fi

# delete / from the end of the string
OUT_DIR=${OUT_DIR%/}

# Check if the OUT_DIR directory exists
if [ ! -d "$OUT_DIR" ]; then
    echo "Error: $OUT_DIR does not exist"
    exit 1
fi

PARSED_DIR=${OUT_DIR}/parsed
echo "Parsed directory: $PARSED_DIR"
# Check if the parsed directory already exists
if [ -d "$PARSED_DIR" ]; then
    echo "Error: "$PARSED_DIR" already exists"
    rm -rf $PARSED_DIR
    #exit 1
fi
mkdir $PARSED_DIR


# Check if the injection rates file exists
if [ ! -f "${OUT_DIR}/injection_rates.txt" ]; then
    echo "Error: "${OUT_DIR}/injection_rates" does not exist"
    exit 1
else
    INJECTION_RATES=$(cat ${OUT_DIR}/injection_rates.txt)
fi

# Check if the routings file exists
if [ ! -f "${OUT_DIR}/routing_functions.txt" ]; then
    echo "Error: "${OUT_DIR}/routings" does not exist"
    exit 1
else
    ROUTINGS=$(cat ${OUT_DIR}/routing_functions.txt)
fi

# Check if the traffics file exists
if [ ! -f "${OUT_DIR}/traffics.txt" ]; then
    echo "Error: "${OUT_DIR}/traffics" does not exist"
    exit 1
else
    TRAFFICS=$(cat ${OUT_DIR}/traffics.txt)
fi

#check if the num_vcs file exists
if [ ! -f "${OUT_DIR}/num_vcs.txt" ]; then
    echo "Error: "${OUT_DIR}/num_vcs" does not exist"
    exit 1
else
    NUM_VCS=$(cat ${OUT_DIR}/num_vcs.txt)
fi

#check if the packet_sizes file exists
if [ ! -f "${OUT_DIR}/packet_sizes.txt" ]; then
    echo "Error: "${OUT_DIR}/packet_sizes" does not exist"
    exit 1
else
    PACKET_SIZES=$(cat ${OUT_DIR}/packet_sizes.txt)
fi

#check if the allocators file exists
if [ ! -f "${OUT_DIR}/allocators.txt" ]; then
    echo "Error: "${OUT_DIR}/allocators" does not exist"
    exit 1
else
    ALLOCATORS=$(cat ${OUT_DIR}/allocators.txt)
fi


for traffic in $TRAFFICS;
do
    mkdir -p ${PARSED_DIR}/${traffic}
    for routing in $ROUTINGS;
    do
        for allocator in $ALLOCATORS; 
        do
            for num_vcs in $NUM_VCS;
            do
                for packet_size in $PACKET_SIZES;
                do
                    for inj_rate in $INJECTION_RATES;
                    do
                    
                    # Parse stats from sim
                    grep "results:" ${OUT_DIR}/${traffic}/${num_vcs}_${inj_rate}_${routing}_${allocator}_${packet_size}.out | sed 's/results://g' >> ${PARSED_DIR}/${traffic}/${num_vcs}_${allocator}_${routing}_${packet_size}_sim_out.csv
                    grep "Time taken is" ${OUT_DIR}/${traffic}/${num_vcs}_${inj_rate}_${routing}_${allocator}_${packet_size}.out| awk '{print $4}' >> ${PARSED_DIR}/${traffic}/${num_vcs}_${allocator}_${routing}_${packet_size}_cycles.csv
                    grep "Total run time" ${OUT_DIR}/${traffic}/${num_vcs}_${inj_rate}_${routing}_${allocator}_${packet_size}.out | awk '{print $4}' >> ${PARSED_DIR}/${traffic}/${num_vcs}_${allocator}_${routing}_${packet_size}_run_time.csv
                    #grep "Total hops average" ${OUT_DIR}/${traffic}/${num_vcs}_${inj_rate}_${routing}_${allocator}_${packet_size}.out | awk '{print $4}' >> ${PARSED_DIR}/${traffic}/${num_vcs}_${allocator}_${routing}_${packet_size}_hops_avg.csv
                    grep "utilization" ${OUT_DIR}/${traffic}/${num_vcs}_${inj_rate}_${routing}_${allocator}_${packet_size}.out | awk '{print $4}' >> ${PARSED_DIR}/${traffic}/${num_vcs}_${allocator}_${routing}_${packet_size}_vc_util.csv

                    done
                    # Delete "results:" from the beginning of the lines
                    sed -i 's/results://g' ${PARSED_DIR}/${traffic}/${num_vcs}_${allocator}_${routing}_${packet_size}_sim_out.csv            
                done  
            done
        done
    done
done



# Delete .err files if they are empty
for file in $(find ${OUT_DIR} -name "*.err");
do
    if [ ! -s $file ]; then
        rm $file
    fi
done

# Copy the injection rates file to the parsed directory
cp ${OUT_DIR}/injection_rates.txt ${PARSED_DIR}/injection_rates.txt
cp ${OUT_DIR}/routing_functions.txt ${PARSED_DIR}/routing_functions.txt
cp ${OUT_DIR}/traffics.txt ${PARSED_DIR}/traffics.txt
cp ${OUT_DIR}/num_vcs.txt ${PARSED_DIR}/num_vcs.txt
cp ${OUT_DIR}/packet_sizes.txt ${PARSED_DIR}/packet_sizes.txt
cp ${OUT_DIR}/allocators.txt ${PARSED_DIR}/allocators.txt

