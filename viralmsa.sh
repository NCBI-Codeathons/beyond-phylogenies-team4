#!/bin/bash
#SBATCH --job-name=viralmsa
#SBATCH --mail-user=brittany.rife@epi.ufl.edu	# Doesn't have to be uf email
#SBATCH --mail-type=END,FAIL					# Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --account=salemi
#SBATCH --qos=salemi							# salemi-b for burst
#SBATCH --nodes=1								# Number of physical computers (can leave out as long as ntasks=1)
#SBTACH --ntasks=1								# Really only need more than one for mpi
#SBATCH --cpus-per-task=4 ## Needs to be N+1	# Number of cores/processors per task
#SBATCH --output=%viralmsa_%j.out
#SBATCH --mem=12gb
#SBATCH --time=10:00:00						# days-hrs:min:sec
##SBATCH --array=1-1000

#####################################################################################

date; pwd; hostname # Prints some info

date1=$(date +"%s")

file=$(ls *gapstripped.fa)


module load viralmsa

ViralMSA.py -s $file -t 4 -e brittany.rife@ufl.edu -o ${file%_gapstripped.fa}_aln -r ../cov_reference/*fasta

#####################################################################################

date
date2=$(date +"%s")
diff=$(($date2-$date1))
echo "$(($diff / 60)) minutes and $(($diff % 60)) seconds elapsed."


