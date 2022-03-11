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
    exit 1
else
    mkdir $PARSED_DIR
fi


# Check if the injection rates file exists
if [ ! -f "${OUT_DIR}/injection_rates" ]; then
    echo "Error: "${OUT_DIR}/injection_rates" does not exist"
    exit 1
else
    INJECTION_RATES=$(cat ${OUT_DIR}/injection_rates)
fi

# Check if the routings file exists
if [ ! -f "${OUT_DIR}/routings" ]; then
    echo "Error: "${OUT_DIR}/routings" does not exist"
    exit 1
else
    ROUTINGS=$(cat ${OUT_DIR}/routings)
fi

# Check if the traffics file exists
if [ ! -f "${OUT_DIR}/traffics" ]; then
    echo "Error: "${OUT_DIR}/traffics" does not exist"
    exit 1
else
    TRAFFICS=$(cat ${OUT_DIR}/traffics)
fi




for traffic in $TRAFFICS;
do
  mkdir -p ${PARSED_DIR}/${traffic}
  for routing in $ROUTINGS;
  do
    for inj_rate in $INJECTION_RATES;
    do
        # Parse stats from sim
        grep "results:" ${OUT_DIR}/${traffic}/sim_${inj_rate}_${routing}.out | sed 's/results://g' >> ${PARSED_DIR}/${traffic}/${routing}_sim_out.csv
        grep "Time taken is" ${OUT_DIR}/${traffic}/sim_${inj_rate}_${routing}.out| awk '{print $4}' >> ${PARSED_DIR}/${traffic}/${routing}_cycles.csv
        grep "Total run time" ${OUT_DIR}/${traffic}/sim_${inj_rate}_${routing}.out | awk '{print $4}' >> ${PARSED_DIR}/${traffic}/${routing}_run_time.csv

    done
    # Delete "results:" from the beginning of the lines
    sed -i 's/results://g' ${PARSED_DIR}/${traffic}/${routing}_sim_out.csv
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
cp ${OUT_DIR}/injection_rates ${PARSED_DIR}/injection_rates
cp ${OUT_DIR}/routings ${PARSED_DIR}/routings
cp ${OUT_DIR}/traffics ${PARSED_DIR}/traffics
