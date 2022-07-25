#!/bin/bash

IN=$1
TMP=${IN}.tmp

grep -v ^submissions $IN | sed 's/; /\t/' | \
while read SEQID GISAID;
do
  SMPID=$(echo $SEQID | cut -d / -f 3)
  echo -e "$SMPID\t$GISAID"
done > $TMP

sqlite3 $FLACO_DB <<EOF
.separator "\t"
.import $TMP SampleGISAID
EOF
