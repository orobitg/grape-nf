#!/bin/bash

#
# A generic wrapper for generating the index file
#
# It takes as input a genome file and generates the index file
#
# CLI Parameters:
# - $1: The mapping tool used
# - $2: The genome FASTA file to format
# - $3: The output index name
# - $4: The number of threads

[ ! -z $NF_DEBUG_SCRIPT ] && set -x
set -e
set -u

GENOME=${2}
INDEX=${3}
CPUS=${4}

case "$1" in
'gem')
gemtools index -i ${GENOME} -o index.gem -t ${CPUS} --no-hash
mv index.gem ${INDEX}
;;

'tophat2')
bowtie2-build ${GENOME} ${INDEX}
;;

'hpg')
mkdir $INDEX
hpg-aligner build-index -g $GENOME -r 10  -i $INDEX
;;

*) echo "Not a valid indexer strategy: $1"; exit 1
;;

esac 
