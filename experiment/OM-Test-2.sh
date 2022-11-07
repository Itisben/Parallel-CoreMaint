#!/bin/bash
# test the scalability

CORE="../core-maint/core2"

TIMEOUT_TIME="1200s" # timeout
WORKER_LIST="0 1 2 4 8 16 24 32 40 48 56 64 128"
INSERT_SIZE="10000000 20000000 30000000 40000000 50000000 60000000 70000000 80000000 90000000 100000000"
#WORKER_LIST="0"

#repeat times
REPEAT_TIME="100"



# misc fields
BENCHDIR=`pwd`
TMPFILE="${BENCHDIR}/test.out" # Create a temporary file to store the output
FAILFOLDER="${BENCHDIR}/failures"

# input graphs folders
OUR_TEST_GRAPH="/home/guob15/test-graph/test"


# results
OUR_RESULT="${BENCHDIR}/results/results-om-scalability-7-22-2022.csv"

#result column I=10,000,000, T 3 (repeated random test), t 21 22 23, 
#l 0 omplock 1 cas lock. w worker# 
COLUMN="I,T,t,l,w,Itime,Otime,Dtime,Mtime,tag#,subtag#,orepeat#,relabel#,maxtag#,morepeat#"

trap "exit" INT #ensures that we can exit this bash program

# create new output file, or append to the exisiting one
create_results() {
    output_file=${1}
    if [ ! -e ${output_file} ]
    then 
        touch ${output_file}
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
    I=${1}
    T=${2}
    t=${3}
    l=${4}
    w=${5} #use in future

    # -T 3 is CSR file -T 4 is CSR file without sampling.
    
    cmd="${TIMEOUT_TIME} ${CORE} -p ${OUR_TEST_GRAPH}/ -I ${I} -T ${T} -t ${t} -l ${l} -w ${w}"
    echo ${cmd}
    
    #timeout ${cmd}

    timeout ${cmd} &> ${TMPFILE}
    python3 OM-parse-output.py  "${I}" "${T}" "${t}" "${FAILFOLDER}" "${TMPFILE}" "${OUR_RESULT}" 
}

do_tests() {
    
    # for workers in `echo ${WORKER_LIST}`
    # do 
    #     #test three cases, cas lock
    #     test_real_offline "${INSERT_SIZE}" "2" "20" "0" "${workers}"
    #     test_real_offline "${INSERT_SIZE}" "2" "21" "0" "${workers}"
    #     test_real_offline "${INSERT_SIZE}" "2" "22" "0" "${workers}"
    #     test_real_offline "${INSERT_SIZE}" "2" "23" "0" "${workers}"

    #     #test three cases, omplock
    #     test_real_offline "${INSERT_SIZE}" "2" "20" "1" "${workers}"
    #     test_real_offline "${INSERT_SIZE}" "2" "21" "1" "${workers}"
    #     test_real_offline "${INSERT_SIZE}" "2" "22" "1" "${workers}"
    #     test_real_offline "${INSERT_SIZE}" "2" "23" "1" "${workers}"
    # done
    
    workers="32"
    for size in `echo ${INSERT_SIZE}`
    do 
        #test three cases, cas lock
        test_real_offline "${size}" "2" "20" "0" "${workers}"
        test_real_offline "${size}" "2" "21" "0" "${workers}"
        test_real_offline "${size}" "2" "22" "0" "${workers}"
        test_real_offline "${size}" "2" "23" "0" "${workers}"

        #test three cases, omplock
        # test_real_offline "${size}" "2" "20" "1" "${workers}"
        # test_real_offline "${size}" "2" "21" "1" "${workers}"
        # test_real_offline "${size}" "2" "22" "1" "${workers}"
        # test_real_offline "${size}" "2" "23" "1" "${workers}"
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
