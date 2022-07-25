#!/bin/bash

OUT="/dev/stdout"
PLOT=""

while getopts "o:p:h" opt; do
    case $opt in
	h)
	    echo "Usage: covdepth [-o outfile] [-p plotfile] bamfile"
	    exit 0
	    ;;
	o)
	    OUT="$OPTARG"
	    ;;
	p)
	    PLOT="$OPTARG"
	    ;;
    esac
done
shift $((OPTIND-1))

BAM=$1

samtools depth -d 0 -a -r MN908947.3 $BAM > $OUT

if [[ "$PLOT" ]];
then
  Plots.py plot -o $PLOT -xc 2 -yc 3 -xs 16 $OUT
fi
