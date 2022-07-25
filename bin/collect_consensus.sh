#!/bin/bash

SAMPLES=$1
OUTFILE=$2

rm -f $OUTFILE

for smp in $(cat $SAMPLES); do 
  cat ${smp}/consensus/*.fasta >> $OUTFILE
done
