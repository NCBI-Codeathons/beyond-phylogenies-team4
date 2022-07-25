#!/bin/bash

BAM=$1
POSITIONS=$2

if [[ -z $BAM ]];
then
  echo "Usage: get-mq.sh BAMFILE POSITIONS"
  echo
  echo "where POSITIONS is a file containing genomic positions, one per line."
  echo 
  echo "Output is to standard output, redirect it to a file with > output.txt."
  echo
  exit 1
fi

for pos in $(cat $POSITIONS);
do
  echo -e -n "${pos}\t"
  samtools view $BAM MN908947.3:${pos}-${pos} | cut -f 5 | tr '\n' '\t'
  echo
done

