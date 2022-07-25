#!/bin/bash

#SBATCH --time=12:00:00
#SBATCH --mem=5G

SESSIONID=$1
OUTDIR=$2

mkdir -p $OUTDIR

pushd $OUTDIR

bs await appsession ${SESSIONID} -f csv -F Name > SESSION_SAMPLES

for smp in $(cat SESSION_SAMPLES);
do
  bs download appresults --name $smp -o $smp
done

popd
