#!/bin/bash

CORE="../core-maint/core"
#compared parallel methods: fast core maintenance 
KCORE="../FasterCoreMaintenance/kcore/kcore"
#compared parallel methods: 
MYKCORE=""

TIMEOUT_TIME="600s" # 3 minutes timeout
#WORKER_LIST="0 1 2 4 8 16 24 32 40 48 56 64 128"
#WORKER_LIST="1 2 4 8 16 32 64"
WORKER_LIST="16"

#sample a fixed number of edges
#EDGE_SIZE="100000 200000 400000 800000 1600000 3200000 6400000 12800000 25600000 51200000 102400000"
#EDGE_SIZE="100000 200000 300000 400000 500000 600000 700000 800000 900000 1000000"
#EDGE_SIZE="500000 600000 700000 800000 900000 1000000"
EDGE_SIZE="100000"

# misc fields
BENCHDIR=`pwd`
TMPFILE="${BENCHDIR}/test.out" # Create a temporary file to store the output
FAILFOLDER="${BENCHDIR}/failures"

# input graphs folders
#OUR_TEST_GRAPH="/media/guob15/WIN/experiment-data/SCC-test-data/TestTrim/large"
#OUR_TEST_GRAPH="/home/guob15/TestTrim/test"
#OUR_TEST_GRAPH="/home/guob15/test-graph/test"
OUR_TEST_GRAPH=${1}

#GRAPHT="5" #repeated random csr
GRAPHT=${2}

# results
#OUR_RESULT="${BENCHDIR}/results/results-10-9-2021-T3-hp.csv"
#GRAPHT="3" #random csr

OUR_RESULT="${BENCHDIR}/results/results-10-7-2022-stability.csv"


#result column
#total number of relabel tags. 
COLUMN="model,N,M,workers,alg,init-mstime,insert-num,remove-num,mstime,V*,V+,S,relabel,tag,subtag,assert,Insert,Delete,date"

#repeat times
REPEAT_TIME="100"

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
    base=${name%.bin} # without the .bin#

    for edge_size in `echo ${EDGE_SIZE}`
    do
        # -T 3 is CSR file -T 4 is CSR file without sampling.
        cmd="${TIMEOUT_TIME} ${CORE} -p ${OUR_TEST_GRAPH}/${name} -I ${edge_size} -m ${alg} -T ${GRAPHT} -w ${workers}"
        echo ${cmd}

        timeout ${cmd} &> ${TMPFILE}
        cat ${TMPFILE}

        python3 ParCM-parse-output.py "${base}" "${alg}" "${workers}" "${FAILFOLDER}" "${TMPFILE}" "${OUR_RESULT}" 

        echo ""

        cmd="${TIMEOUT_TIME} ${CORE} -p ${OUR_TEST_GRAPH}/${name} -R ${edge_size} -m ${alg} -T ${GRAPHT} -w ${workers}"
        echo ${cmd}

        timeout ${cmd} &> ${TMPFILE}
        cat ${TMPFILE}

        python3 ParCM-parse-output.py "${base}" "${alg}" "${workers}" "${FAILFOLDER}" "${TMPFILE}" "${OUR_RESULT}" 


        echo ""
        echo ""
        echo ""

        # generate the graph.
        cmd="${TIMEOUT_TIME} ${CORE} -t 31 -p ${OUR_TEST_GRAPH}/${name} -T ${GRAPHT} -I ${edge_size}" 
        timeout ${cmd}

        cmd="${TIMEOUT_TIME} ${KCORE} ${OUR_TEST_GRAPH}/${name}".edge" ${OUR_TEST_GRAPH}/${name}".edge-${edge_size}" ${workers}"
        echo ${cmd}

        timeout ${cmd} &> ${TMPFILE}
        cat ${TMPFILE}

        #first compared methods
        python3 ParCM-parse-output.py "${base}" "41" "${workers}" "${FAILFOLDER}" "${TMPFILE}" "${OUR_RESULT}" 

        echo ""
        echo ""
        echo ""
        echo ""
    done 

    
}


# test_all_real_offline ALG WORKERS
test_all_real_offline() {
    alg=${1}
    IR=${2}
    workers=${3}
    other=${4}
    if [ ! "$(ls -A ${OUR_TEST_GRAPH})" ]; then
        echo ${OUR_TEST_GRAPH} "is empty"
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
  
        test_all_real_offline "4" "-IorR" "${workers}" "-l 0 -t 0" #ours parallel remove

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
