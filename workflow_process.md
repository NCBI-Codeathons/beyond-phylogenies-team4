# Workflow for addition of new sequence data to existing cluster phylogenies

## Dependencies
* [python3](https://www.python.org)
* [R v4.2.1](https://cloud.r-project.org/)
* [ViralmMSA](https://github.com/niemasd/ViralMSA)
* [mafft (optional)](https://mafft.cbrc.jp/alignment/software/)
* [Pangolin](https://cov-lineages.org/resources/pangolin/installation.html)
* [FastTreeMP](http://www.microbesonline.org/fasttree/#Install)
* [UShER](https://usher-wiki.readthedocs.io/en/latest/Installation.html)


## Step-by-step (user)


1. Create directory for initial sequencing run to initiate sub-project (**Must be in the format of <NAME_DATE>, where the name is the name of the sequencing run and the date is today's date in yyy-mm-dd**). Then copy FASTA with generated consensus sequences into this directory (**this must be done manually, cannot be automated**):
	
	```	
	mkdir <NAME_DATE>;
	cd <NAME_DATE>;
	cp /path/to/consensus/sequences/FASTA ./<NAME_DATE>.fasta;
	 
	```
	
2. Create a folder with the most up-to-date GISAID msa from gisaid.org (**this used to be done automatically using the following command, but there has been a change to the website. This website may have an answer to overcome this issue (https://askubuntu.com/questions/615142/wget-how-to-download-file-from-a-web-page-that-prompts-you-to-click-on-agree-c**) .

	```
	mkdir gisaid_msa;
	MSADATE=<Date (mmdd) of most recent msa>;
	wget --user <username> --password <password> https://www.epicov.org/epi3/entities/tmp/tmp_sd_2022_07_27_02_21_rfmbn1_2rih773b9433/msa_${MSADATE}.tar.xz;
	tar -xvf msa*;
	cd ../;
	```

3. Initialize variables based on run name, earliest and latest sample collection dates:

	```
	run=<NAME>;
	DATE1=<earliest sample collection date (yyyy-mm-dd)>;
	DATE2=<latest sample collection date (yyyy-mm-dd)>;
	```

## Step-by-step (automated)

1. Set variables based on last date of workflow run (already in dates.txt file) and today's date:

	```
	previous=$( tail -n 1 dates.txt );
	today=$( date +"%Y-%m-%d" )
	echo ${today} >> dates.txt;
	```

2. Preliminarily add samples to all previous trees using UShER in order to evaluate thir best-fitting cluster origin (creating new temporary pb files for today's date):
	
	```
	for i in $(ls | grep "${previous}.pb"); 
	do usher -i ${i} -v ${run}_${today}/${run}_florida_gisaid_${today}_masked.vcf -o -p ${i%_${previous}.pb}_${today}_prelim.pb; done	
	```

3. Compute parsimony score for new samples (vcf) assigned to each annotated tree (.pb files):
	
	```
	for i in $(ls | grep "${previous}.pb"); do 
	mkdir ${i%_${previous}.pb}_${today};
	usher -i ${i} -v {run}_${today}/${run}_florida_gisaid_${today}_masked.vcf -p -d ${i%_${previous}.pb}_${today};
	done
	```
	
4. Previous step will generate parsimony scores in tsv files, which we need to rename and copy to single folder for R analysis:
	
	```
	mkdir BPS_${today};
	for i in ./*_${today}; do mv ${i}/parsimony-scores.tsv BPS_${today}/${i}_parsimony-scores.tsv; done
	```
	
5. Use R script to evaluate parsimony scores from tsv files and output new fasta files for sequences needing to be placed on trees (currently not working owing to change in usher output):
	
	```
	cd BPS_${today}
	Rscript branch_support_eval.R
	```	


6. Add reference sequence again to each of the newly generated fasta files (from previous step):
	
	```
	for i in ./*.fasta; do cat ../cov_reference/cov_reference.fasta >> ${i}; done
	```

7. Place updated fasta in new folders:
	
	```
	for i in ./*updated.fasta; do mkdir .${i%.fasta}; cp ${i} .${i%.fasta}; done;
	cd ../
	``
	
8. Create a list of files that need to be updated for easier scripting:
	
	```
	find . -type d -name "*${today}_updated" > source.txt
	```

	
9. Convert updated fastas to vcf:
	
	 ```
	for i in $(cat source.txt); do faTovcf -ref=MN908947.3 ${i}/*.fasta ${i}/*.vcf; done	
	```

10. Officially add samples to the cluster phylogenies they belong (creating new, final pb files for today):
	
	```
	for i in $(cat source.txt); 
	do usher -i ${i%_${today}*}_${previous}.pb -v ${i}/*.vcf -u -d -o ${i%_updated}.pb; done
	```	
	
11. Not all trees (.pb) are going to be updated with sequences, so need to find all trees without today's date and make a copy of previous pb with today's date:
	
	```
	for i in ./*${previous}.pb; 
	do mv -vn ${i} ${i%_${previous}*}_${today}.pb; done
	```

12. Move new unannotated tree (.nh) files generated from the previous step into one folder (as well as original combined sample fasta) for characterization:
	
	```
	mkdir "trees_${today}";
	for i in $(cat source.txt); do cp ${i}/*final-tree.nh trees_${today}/${i}.tree; done;
	cp ?????? trees_${today}
	```
	
13. Run modified DYNAMITE to search for (and prune) clusters within background tree: 
	
	```
	cd ./background_${today}
	Rscript dynamite_background.R
	```	

14. Run R script to characterize added sequences (add new folder to save discard results):	*This needs to be modified so that when discard tree exists, it will create just a fasta that will need to be processed*
	
	```
	mkdir "./discard_${today}"
	cd trees_${today}
	Rscript ../fitness_calc.R
	```

Analysis of the discarded sequences currently not included but may be added later (see bottom of page).	

15. Place pruned background sequences onto previous annotated tree (*Brittany needs to re-evaluate this description*):
	
	```
	cd ../background_${today}_updated;
	rm -v !(*.fasta);
	cd ../;
	python Fasta2UShER.py -inpath ./background_${today}_updated -output ./background_${today}.vcf -reference ./cov_reference/cov_reference.fasta;
	usher -v background_${today}.vcf -i ./background_${previous}.pb -T 4 -c -u -d ./background_${today}_updated -o background_${today}.pb;
	```
	
16. Concatenate previous fasta with updated fasta for tree optimization :
	
	```
	for i in $(cat source.txt); do 
tail -n +2 ${i}/*.fasta >> ../full.fasta
	done
	```
	
17. FastTree is very picky about duplicate sequences, which can be generated from past regional sequences sharing best hits in FLACO BLAST with newly added regional sequences. UShER ignores new sequences if they are already in the tree, but they are still present in the FASTAs. Whereas I have seqkit to do this here, there are any number of tools that can do this easily:
	
	```
	for i in $(cat source.txt); cat ${i}/full.fasta | seqkit rmdup -n -o ${i}/clean.aln
	```

18. Force bifurcating tree in R (required for tree optimization):
	
	```
	for i in $(cat source.txt); do 
	Rscript ./bifurcate.R -t ${i}/uncondensed-final-tree.nh; done
	```

 
19. Re-optimize trees (only those that were updated) using FastTree: 
	
	```
	for i in $(cat source.txt); do FastTreeMP -nt -gamma -nni 0 -spr 2 -sprlength 1000 -boot 100 -log fasttree.log -intree ${i}/*.nh.tree ${i}/clean.aln > ${i}/reop.tree; done
	```
	

20. Re-annotate (UShER) updated fasttree trees, replacing old ones:
	
	```
	for i in ${cat source.txt); do
  	usher -v ${i%_updated}.vcf -t ${i}/reop.tree -T 4 -c -u -d ${i} -o ${i%_updated}.pb;
  	done
  	```

21. Test background and discard trees for clusters (separate R script using actual branchwise algorithm). R script will write new results for discard pile to a new folder, but will write results from background to same folder.

	If no clusters are found, results are not written, but downstream steps require that the new discard folder exist, so need to copy over.

	```
	mkdir ./discard_${today}_updated
	cd ./discard_${today}
	Rscript branchwise.R -t *.treefile -p ../discard_${today}_updated
	
	cd ../background_${today}_updated
	Rscript branchwise.R -t reop.tree -p ./
	cd ../
	```
		
22. Process fasta and tree for background (will require determining if new files added as a result of R script in step above)

	```
	cat ./cov_reference/cov_reference.fasta >> ./background_${today}_updated/*.fasta
	find . -type f -name "background_${today}_updated.tree" -empty -exec cp ./background_${today}/reop.tree ./background_${today}_updated.tree
	python3 Fasta2UShER.py -inpath background_${today}_updated -output background_${today}.vcf -reference ./cov_reference/cov_reference.fasta	
	usher -v background_${today}.vcf -t ./background_${today}_updated.tree -T 4 -c -u -d ./background_${today}_updated -o background_${today}.pb
	```
	
23. Add reference sequence to new cluster fastas and place in new directories:
	
	```
	for i in ./*.fasta; do  cat ./cov_reference/cov_reference.fasta >> ${i}; cp ${i} .${i%.fasta}; done
	```
	
24. Update source.txt file with new clusters and process fastas:
	
	```
	for i in ./*.fasta; do mkdir .${i%.fasta}; echo ${i%.fasta} > source.txt; done	
	ml python
	for i in $(cat source.txt); do python Fasta2UShER.py -inpath ${i} -output ${i}.vcf -reference ./cov_reference/cov_reference.fasta; done
	```
	
25. Perform subtree pre-processing for new clusters:
	
	```
	ml usher	
	for i in $(cat source.txt); do usher -v ${i}.vcf -t ${i}.tree -T 4 -c -u -d ./${i} -o ${i}.pb; done
	```

## Plotting

1. Plot summary of variants over time (needs modification)

```
meta=$(ls updated_metadata*)
Rscript variant_plotting.R -m ${meta} -v voc.tab
```

1. Plot new mutations for each cluster (relative to reference strain) *Needs modification to utilize vcf instead of MSA*. For mutations relative to corresponding lineage, see section below:

```
meta=$(ls updated_metadata*)
for i in *fasta; do
Rscript mut_analysis.R -q $i -m ${meta} -d ${today};
done
```


## For mapping new mutations (relative to corresponding lineage)

1. Clone git repository for pango reference sequences (not needed now)
	
	```
	mkdir ref_lineages
	cd ref_lineages
	git clone https://github.com/cov-lineages/pango-designation.git
	```

2. Extract sequence names only from pango reference file (not needed now)

	```
	awk -F"," `{print $1}`	pango-designation/lineages.csv | tail -n +2 > ref_lineages_${today}.txt
	```

3. Retrieve pango reference sequences from most up-to-date gisiad alignment (not needed now)

	```
	seq_cleaner.py -P ref_lineages_${today}.txt <path/to/gisaid/msa> > ref_lineages_${today}.fa
	```
	
4. Gap strip for alignment of pango reference sequences (not needed now)

	```
	seq_cleaner.py -g ref_lineages_${today}.fa > ref_lineages_${today}_gapstripped.fa
	```
	
5. Align pango reference sequences using viralmsa (not needed now)

	```
	file=$(ls *gapstripped.fa)
	ml viralmsa
	ViralMSA.py -s $file -t 4 -e brittany.rife@ufl.edu -o ${file%_gapstripped.fa} -r ../cov_reference/*fasta
	
	```

6. Mask uncertain sites in pango references (not needed now)

	```
	mask_aln_using_vcf.py -i ref_lineages_${today}/*.aln -o ref_lineages_${today}_masked.aln -v ../problematic_sites_sarsCov2.vcf
	```

7. Create R list of lineage-specific alignments for either whole genome (makeRefList.R) or Spike (not needed now)

	```
	R makeRefSpikeList.R -r ref_lineages_${today}_masked.aln	
	```		
	
	
## For processing of discarded sequences (i.e., tree reconstruction and cluster picking)

Below is disregarded for now, but will be added back in again later.

40. (Skip this step for now) If discard tree not present, create one. If so, process fasta from previous step and add to existing annotated tree, and also update source.txt file (sh scripts need to be included or the following code modified by Brittany):
	
	```
	if [ ! -e "discard_${previous}.vcf" ]; then
       cd 	../discard_${today}; sbatch ../iqtree.sh;
      else
    	cat ./cov_reference/cov_reference.fasta >> ./discard_${today}_updated/*.fasta
		ml python
		python Fasta2UShER.py -inpath discard_${today}_updated -output discard_${today}.vcf -reference ./cov_reference/cov_reference.fasta	
		ml usher
		usher -v discard_${today}.vcf -t ./discard_${today}_updated.tree -T 4 -c -u -d ./discard_${today}_updated -o discard_${today}.pb
		echo "discard_${today}_updated" >> source.txt
	fi 
	```

41. If first condition met above, then need to pre-process tree and fasta (but somehow need to make sure iqtree is finished:
	    
	
	```
	cat ./cov_reference/cov_reference.fasta >> ./discard_${today}_updated/*.fasta
	python Fasta2UShER.py -inpath discard_${today}_updated -output discard_${today}.vcf -reference ./cov_reference/cov_reference.fasta	usher -v discard_${today}.vcf -i ./discard_${previous}.pb -T 4 -c -u -d ./discard_${today}_updated -o discard_${today}.pb
	```
	






