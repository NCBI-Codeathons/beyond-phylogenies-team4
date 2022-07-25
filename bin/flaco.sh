#!/bin/bash

# Florida COVID Sequencing Pipeline
# A. Riva, B. Rife-Magalis, University of Florida

# Configuration
BASESPACE=/home/ariva/basespace.icbr/
CLEAN_NFRAC=30
CLEAN_MLEN=0

# Script variables
SCRIPT_HOME=$(dirname $(readlink -f $0))
PROJ=$1
CONSENSUS=${PROJ}.consensus.fasta
CLEANCONS=${PROJ}.consensus.clean.fasta
MSA_OUTPUT="msa_output"
MSA_EMAIL="ariva@ufl.edu"

function get_consensus_sequences() {
  rm -f $CONSENSUS
  PROJDIR=${BASESPACE}/Projects/${PROJ}/AppResults/
  echo Copying consensus fastq files from $PROJDIR
  for smp in ${PROJDIR}/*/; do
    for fasta in ${smp}/Files/consensus/*.fasta; do
      echo $fasta
      cat $fasta >> $CONSENSUS
    done
  done
}

function clean_consensus_sequences() {
  rm -f $CLEANCONS
  $SCRIPT_HOME/seq_cleaner.py -f $CLEAN_NFRAC -l $CLEAN_MLEN $CONSENSUS $CLEANCONS
}

function run_msa() {
  module load viralmsa
  ViralMSA.py -s $CLEANCONS -r SARS-CoV-2 -o $MSA_OUTPUT -e $MSA_EMAIL
}

get_consensus_sequences
NS=$(grep -c ^ $CONSENSUS)
echo "Sequences in consensus: $NS"
echo "Cleaning sequences..."
clean_consensus_sequences
NC=$(grep -c ^ $CLEANCONS)
echo "Cleaned sequences: $NC"
echo "Running ViralMSA..."
run_msa
