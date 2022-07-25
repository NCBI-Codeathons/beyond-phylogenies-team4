#!/bin/bash

while read smp path;
do
  f1=$(ls ${path}/*_R1_*.fastq.gz)
  f2=$(ls ${path}/*_R2_*.fastq.gz)
  echo -e "${smp}\t${PWD}/${f1}\t${PWD}/${f2}"
done < $1
