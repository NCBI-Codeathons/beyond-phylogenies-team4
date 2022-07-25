#!/bin/bash
#SBATCH --job-name=mut_analysis
#SBATCH --mail-user=brittany.rife@epi.ufl.edu	# Doesn't have to be uf email
#SBATCH --mail-type=END,FAIL					# Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --account=salemi
#SBATCH --qos=salemi							# salemi-b for burst
#SBATCH --nodes=1								# Number of physical computers (can leave out as long as ntasks=1)
#SBTACH --ntasks=1								# Really only need more than one for mpi
#SBATCH --cpus-per-task=4 ## Needs to be N+1	# Number of cores/processors per task
#SBATCH --output=%mut_analysis_%j.out
#SBATCH --mem=120gb
#SBATCH --time=10:00:00						# days-hrs:min:sec
##SBATCH --array=1-1000

#####################################################################################

date; pwd; hostname # Prints some info

date1=$(date +"%Y%m%d")

today=$(tail -n 1 dates.txt)
aln=$(ls *${today}_masked.aln)
meta=$(ls *lineages_${today}.csv)

ml R/4.0

Rscript mut_analysis.R -q ${aln} -m ${meta} -d 2021-12-01


#####################################################################################

date
date2=$(date +"%Y%m%d")
diff=$(($date2-$date1))
echo "$(($diff / 60)) minutes and $(($diff % 60)) seconds elapsed."


