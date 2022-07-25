#!/bin/bash

# Samples should be a tab-delimited file with two columns:
# SampleName     DirName

RUN=$1
SAMPLES=$2
OUTDIR=$3

if [[ -z $FLACO_DB ]];
then
  echo FLACO environment not initialized - please execute \"source FLACO\"
  exit 1
fi

if [[ "$RUN" == "" ]];
then
  echo Usage: collect-stats.sh RUN SAMPLES OUTDIR
  echo
  echo "SAMPLES should be a tab-delimited file with two columns:"
  echo " - Sample name"
  echo " - Directory"
  exit 1
fi

mkdir -p $OUTDIR

# Run the following in a subshell because samtools module messes with python
(module load samtools;
while read sample smpdir;
do
  echo $sample
  smpname=${smpdir%_*}
  BAM=${smpdir}/Files/${smpname}_tumor.bam
  if [[ -f $BAM ]];
  then
     samtools idxstats $BAM | bam_stats_summary.py > $OUTDIR/${sample}_mapping.txt
     mapdata=$(samtools idxstats $BAM | bam_stats_summary.py -s)
     app_stats.py ${RUN} ${sample} ${smpdir}/Files $mapdata >> $OUTDIR/all-stats.txt
  else
    echo Warning: no BAM file found for sample $sample
  fi

  F1=$OUTDIR/${sample}.${RUN}.plotcov.png
  F2=$OUTDIR/${sample}.pathogen-coverage-report.tsv 
  cp ${smpdir}/Files/*.png $F1 2>/dev/null
  cp ${smpdir}/Files/*.pathogen-coverage-report.tsv $F2 2>/dev/null
  chmod 664 $F1 $F2
done < $SAMPLES
)

flacodb.py stats -r ${RUN} -f w3 -o $OUTDIR/${RUN}-stats.html
flacodb.py stats -r ${RUN} -f csv -o $OUTDIR/${RUN}-stats.txt

flacodb.py cov -r ${RUN} -f w3 -o $OUTDIR/${RUN}-cov.html
flacodb.py cov -r ${RUN} -f csv -o $OUTDIR/${RUN}-cov.txt

echo "HTML reports written to:"
echo "  $OUTDIR/${RUN}-stats.html"
echo "  $OUTDIR/${RUN}-cov.html"
echo "Tab-delimited reports written to:"
echo "  $OUTDIR/${RUN}-stats.txt"
echo "  $OUTDIR/${RUN}-cov.txt"
