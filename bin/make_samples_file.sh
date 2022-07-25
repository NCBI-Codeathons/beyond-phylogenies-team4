#!/bin/bash

# With no args, sample ID is second field of directory name
# Otherwise, specify field as first argument.

if [[ "$1" == "" ]];
then
  FIELD=2
else
  FIELD=$1
fi

for dir in */;
do
  smp=$(echo $dir | cut -f $FIELD -d _)
  echo -e "${smp}\t${dir}"
done | sort
