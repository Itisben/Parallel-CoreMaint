#!/bin/bash

CORE="../core-maint/core2"

TIMEOUT_TIME="1800s" # 200 minutes timeout
#WORKER_LIST="1 2 4 6 8 10 12 14 16 32"
WORKER_LIST="1"

BATCH_SIZE="100000" # the number of inserted or removed edges. 

# misc fields
BENCHDIR=`pwd`
TMPFILE="${BENCHDIR}/test.out" # Create a temporary file to store the output
FAILFOLDER="${BENCHDIR}/failures"

# input graphs folders
#OUR_TEST_GRAPH="/media/guob15/WIN/experiment-data/SCC-test-data/TestTrim/large"
#OUR_TEST_GRAPH="/home/guob15/TestTrim/test"
OUR_TEST_GRAPH="/home/guo/Documents/test"

# results
OUR_RESULT="${BENCHDIR}/results/results-9-15-2021-hp-core2.csv"

#result column
COLUMN="model,N,M,workers,alg,init-mstime,insert-num,remove-num,mstime,V*,V+,S,tag,assert,date"

#repeat times
REPEAT_TIME="15"

trap "exit" INT #ensures that we can exit this bash program

# create new output file, or append to the exisiting one
create_results() {
    output_file=${1}
    if [ ! -e ${output_file} ]
    then 
        touch ${output_file}
        #add column in output csv.
        echo "${COLUMN}" > ${output_file}
    fi
}

# create necessary folders and files
init() {
    if [ ! -d "${FAILFOLDER}" ]; then
      mkdir "${FAILFOLDER}"
    fi
    if [ ! -d "${BENCHDIR}/results" ]; then
      mkdir "${BENCHDIR}/results"
    fi
    touch ${TMPFILE}
    create_results ${OUR_RESULT}
}

# test_real_offline MODEL ALG WORKERS
test_real_offline() {
    name=${1}
    alg=${2}
    IR=${3}
    workers=${4} #use in future

    # -T 3 is CSR file -T 4 is CSR file without sampling.
    
    cmd="${TIMEOUT_TIME} ${CORE} -p ${OUR_TEST_GRAPH}/${name} ${IR} ${BATCH_SIZE} -m ${alg} -T 3"
    echo ${cmd}
    timeout ${cmd} &> ${TMPFILE}
    #timeout ${cmd}

    base=${name%.bin} # without the .bin#
    python3 parse-output.py "${base}" "${alg}" "${workers}" "${FAILFOLDER}" "${TMPFILE}" "${OUR_RESULT}" 
}


# test_all_real_offline ALG WORKERS
test_all_real_offline() {
    alg=${1}
    IR=${2}
    workers=${3}
    if [ ! "$(ls -A ${OUR_TEST_GRAPH})" ]; then
        echo "graphs/real-offline folder is empty"
    else
        for file in ${OUR_TEST_GRAPH}/*.bin
        do
            name=${file##*/} # remove preceding path
            test_real_offline ${name} ${alg} ${IR} ${workers}
        done
    fi
}

do_tests() {
    
    for workers in `echo ${WORKER_LIST}`
    do
        #test_all_real_offline "1" "-I" "${workers}" #traversal,  
        test_all_real_offline "2" "-I" "${workers}" #glist
        test_all_real_offline "3" "-I" "${workers}" #ours

        #test_all_real_offline "1" "-R" "${workers}" #traversal,  
        test_all_real_offline "2" "-R" "${workers}" #glist
        test_all_real_offline "3" "-R" "${workers}" #ours
    done
}


# do_all_experiments ITERATIONS
do_all_experiments() {
    iterations=${1}
    for experiments in `seq 1 ${iterations}`
    do
        echo ""
        echo "iteraton ${experiments}"
        do_tests
    done
}

############################################################

# initialize
init

# performs all experiments N times
do_all_experiments ${REPEAT_TIME}

# cleanup
#rm "${TMPFILE}"
