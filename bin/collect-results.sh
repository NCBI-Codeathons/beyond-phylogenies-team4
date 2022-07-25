#!/bin/bash

if [[ -z $FLACO_DB ]];
then
  echo FLACO environment not initialized - please execute \"source FLACO\"
  exit 1
fi

MINCOVERAGE=70			# Required % genome coverage
FIELD=HiDepthFrac20

# This script is meant to be run in the top-level directory produced
# by the RNA Pathogen Detection App.

# First arg is the project name
# Second arg is a file containing the mapping between sample names and directory names

NAME=$1
SAMPLES=$2

if [[ "$NAME" == "" ]];
then
  echo "Usage: collect-results.sh RUN SAMPLES"
  echo
  echo "SAMPLES should be a tab-delimited file with two columns:"
  echo " - Sample name"
  echo " - Directory"
  exit 1
fi

source /apps/dibig_tools/1.0/lib/sh/utils.sh

# Figure out if result dirs have a "Files" subdir or not

dir1=$(cut -f 2 $SAMPLES | head -1)
if [[ -d $dir1/Files ]];
then
  FDIR="Files/"
else
  FDIR=""
fi

# Begin by copying all consensus fasta files into consensus/ dir

RAW_CONSENSUS="${NAME}-consensus.raw.fa"
GOOD_CONSENSUS="${NAME}-consensus.good.fa"
BAD_CONSENSUS="${NAME}-consensus.bad.fa"
ALN_CONSENSUS="${NAME}-consensus.good.aln.fa"
MASKED_CONSENSUS="${NAME}-consensus.good.aln.mask.fa"
LINEAGES="${NAME}-lineages.csv"

mkdir -p consensus
rm -f consensus/$RAW_CONSENSUS

echo Collecting consensus sequences
consfastas=""
while read sample smpdir;
do
  cons=$(ls ${smpdir}/$FDIR/consensus/*.fasta | head -1)
  if [[ -f $cons ]];
  then
      consfastas="$consfastas $sample:$cons"
  fi
done < $SAMPLES

fix_fasta_names.py $METADATA $consfastas > consensus/$RAW_CONSENSUS
echo Consensus sequences written to consensus/$RAW_CONSENSUS

# Split consensus sequences into good and bad :)

pushd consensus

flacodb.py -r $NAME -o GOOD_SAMPLES.txt query "${FIELD}>${MINCOVERAGE}" 

#seq_cleaner.py -f ${MINCOVERAGE} -r ${NAME}-report.txt -b ${BAD_CONSENSUS} $RAW_CONSENSUS ${GOOD_CONSENSUS}
seq_cleaner.py -w GOOD_SAMPLES.txt -r ${NAME}-report.txt -b ${BAD_CONSENSUS} $RAW_CONSENSUS ${GOOD_CONSENSUS}

# Produce multiple sequence alignment (in a subshell beause of module load)
#echo Aligning consensus sequences with ViralMSA
#(module load viralmsa
# rm -fr alignment
# ViralMSA.py -s ${GOOD_CONSENSUS} -t 4 -e ariva@ufl.edu -o alignment -r MN908947
# cp alignment/*.fa.aln ${ALN_CONSENSUS}
#)

# Produce multiple sequence alignment with mafft
echo Aligning consensus sequences with MAFFT
tmp=$(mktemp -p .)
submit -done mafft.@.done mafft.qsub ${GOOD_CONSENSUS} $tmp ${FLACO_DATA}/MN908947.3.fa
wait_for mafft 1
remgaps.py -r $tmp > ${ALN_CONSENSUS}
rm -f $tmp

# Mask problematic sites (in a subshell beause of module load)
#echo Masking problematic sites
#(module purge;
# module load python3;
# python ${FLACO_HOME}/bin/mask_aln_using_vcf.py -i ${ALN_CONSENSUS} -o ${MASKED_CONSENSUS} -v ${FLACO_DATA}/problematic_sites_sarsCov2.vcf
#)

# Run Pangolin on masked consensus
echo Running Pangolin on ${ALN_CONSENSUS}
(module purge;
 module load pangolin;
 pangolin --outfile ${LINEAGES} ${ALN_CONSENSUS}
)

# Parse pangolin output
parse-pangolin.py ${NAME} ${LINEAGES}

# Write pangolin reports
flacodb.py results -r ${NAME} -f w3 -o ${NAME}-results.html
flacodb.py results -r ${NAME} -f csv -o ${NAME}-results.txt

flacodb.py lintable -r ${NAME} -f w3 -o ${NAME}-lintable.html
flacodb.py lintable -r ${NAME} -f csv -o ${NAME}-lintable.txt

popd
echo "Clean consensus file: consensus/${GOOD_CONSENSUS}"
echo "Masked consensus file: consensus/${MASKED_CONSENSUS}"
echo "Lineages: consensus/${LINEAGES}"
echo "Report (html): ${NAME}-results.html"
echo "Report (txt): ${NAME}-results.txt"
