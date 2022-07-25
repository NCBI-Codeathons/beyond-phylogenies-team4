#!/bin/bash

module purge
module load python/3.8 nextflow

SCRIPT_HOME=$(dirname $(readlink -f $0))
shopt -s nullglob
set -e

# Number of sequences in each BLAST job
MAX_BLAST=25

CMD=$1

if [[ "$CMD" == "makedb" ]];
then
  FASTA=$2
  shift 2
  REMOVE=$*

  if [[ "$FASTA" == "" ]];
  then
    echo "Usage: $0 makedb msa.fasta [patterns]"
    exit 1
  fi

  mkdir -p split-by-month
  FNAME=$(basename $FASTA)
  CLEAN=${FNAME%.*}.clean.fa
  EXCLUDED=${FNAME%.*}.target.fa
  echo "Cleaning alignment $FASTA"
  if [[ "$REMOVE" != "" ]]; then
    echo "Removing sequences containing \"$REMOVE\""
  fi
  ${SCRIPT_HOME}/flacotools.py splitdb $FASTA $CLEAN split-by-month $EXCLUDED $REMOVE

  echo "Generating BLAST database"
  find split-by-month -name DB > split-by-month/DBLIST
  nextflow run ${SCRIPT_HOME}/flacoblast.nf --cmd makedb --dbfile split-by-month/DBLIST
  nextflow clean
  echo "Database creation complete. Please use file $CLEAN in the run step."
  exit 0
fi

if [[ "$CMD" == "split" ]];
then
    shift
    mkdir -p split-by-month
    find split-by-month -name \*.fa | xargs rm -f
    for FASTA in $*;
    do
	echo "Splitting file $FASTA by month..."
	${SCRIPT_HOME}/flacotools.py split $FASTA split-by-month
    done
    exit 0
fi

if [[ "$CMD" == "blast" ]];
then
  INDEX=DB
  echo "Starting BLAST jobs..."
  rm -f blast.cmd
  for spdir in split-by-month/2*/;
  do
      find $spdir -name \*.fa | xargs -n $MAX_BLAST echo -e "$spdir\t$spdir/${INDEX}\t" >> blast.cmd
  done
  nextflow run ${SCRIPT_HOME}/flacoblast.nf --cmd blast --blastcmd blast.cmd
  #nextflow clean
  exit 0
fi

if [[ "$CMD" == "parse" ]];
then
  OUTFILE=$2
  MATCHES=$3
  ${SCRIPT_HOME}/flacotools.py parse split-by-month $OUTFILE $MATCHES
  exit 0
fi

if [[ "$CMD" == "extract" ]];
then
  FASTA=$2
  MATCHES=$3
  OUTFILE=$4
  touch ALLMATCHES.txt
  ${SCRIPT_HOME}/flacotools.py extract -a ALLMATCHES.txt -o $OUTFILE $FASTA $MATCHES
  exit 0
fi

if [[ "$CMD" == "run" ]];
then
  R_TARGETFASTA=$2
  R_DBFASTA=$3

  if [[ "$R_TARGETFASTA" == "" || "$R_DBFASTA" == "" ]];
  then
    echo "Usage: flaco_blast.sh run target.fa database.clean.fa"
    exit 1
  fi
  R_OUTFILE=${R_TARGETFASTA%.*}.out.txt
  R_MATCHES=${R_TARGETFASTA%.*}.matches.txt
  R_OUTFASTA=${R_TARGETFASTA%.*}.out.fa
  $0 split $R_TARGETFASTA
  $0 blast
  $0 parse $R_OUTFILE $R_MATCHES
  $0 extract $R_DBFASTA $R_MATCHES $R_OUTFASTA
  exit 0
fi

echo "flaco_blast.sh - BLAST pipeline for GISAID sequences."
echo
echo "Usage: flaco_blast.sh command [options...]"
echo
echo "where command is one of: run, makedb, split, blast, parse, extract."
echo
