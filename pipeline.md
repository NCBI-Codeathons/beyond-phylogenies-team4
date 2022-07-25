
# RETRIEVE REFERENCE LINEAGES USED BY PANGOLIN



1. Initiate FLACO environment (or simply add bin to PATH)

	```
	source bin/FLACO 
	```
	
2. Download most up-to-date gisaid alignment and place in "gisaid_msa" directory. Below no longer works for some reason, so have to do it manually.

	```
	mkdir gisaid_msa
	wget --user <username> --password <password> https://www.epicov.org/epi3/msa_dddd.tar.xz
	tar -xvf msa*
	cd ../
	```


3. Create a file with dates for later referencing and a variable ("today") from the dates file to help with naming:
	
	```
	echo $(date +"%Y-%m-%d") > dates.txt
	today=$( tail -n 1 dates.txt )
	```

4. Clone git repository for pango reference sequences (not needed now)
	
	```
	mkdir ref_lineages
	cd ref_lineages
	git clone https://github.com/cov-lineages/pango-designation.git
	```

5. Extract sequence names only from pango reference file (not needed now)

	```
	awk -F"," `{print $1}`	pango-designation/lineages.csv | tail -n +2 > ref_lineages_${today}.txt
	```

6. Retrieve pango reference sequences from most up-to-date gisiad alignment (not needed now)

	```
	seq_cleaner.py -P ref_lineages_${today}.txt <path/to/gisaid/msa> > ref_lineages_${today}.fa
	```
	
7. Gap strip for alignment of pango reference sequences (not needed now)

	```
	seq_cleaner.py -g ref_lineages_${today}.fa > ref_lineages_${today}_gapstripped.fa
	```
	
8. Align pango reference sequences using viralmsa (not needed now)

	```
	file=$(ls *gapstripped.fa)
	ml viralmsa
	ViralMSA.py -s $file -t 4 -e brittany.rife@ufl.edu -o ${file%_gapstripped.fa} -r ../cov_reference/*fasta
	
	```
	
9. Download most recent masked sites vcf (needs to be automated so that new file replaces old file instead of being renamed):

	```
	rm problematic_sites_sarsCov2.vcf
	wget - O ../problematic_sites_sarsCov2.vcf https://raw.githubusercontent.com/W-L/ProblematicSites_SARS-CoV2/master/problematic_sites_sarsCov2.vcf
	```

10. Mask uncertain sites in pango references (not needed now)

	```
	ml python
	python ../mask_aln_using_vcf.py -i ref_lineages_${today}/*.aln -o ref_lineages_${today}_masked.aln -v ../problematic_sites_sarsCov2.vcf
	```

11. Create R list of lineage-specific alignments for spike (not needed now)

	```
	sbatch makeRefSpikeList.sh -r ref_lineages_${today}_masked.aln	
	```
	 	
		
**Note:** If this already exists, just need to get sequences not already included in alignment above.


# ALIGNMENT OF INITAL DATA ##############################################################

1. Meanwhile, create directory for initial sequencing run and extract sequences from desired region from database that correspond to the relevant time frame:
	
	```
	RUN=<Sequencing run number/name>
	MSADATE=<Date (mmdd) of most recent gisaid msa deposit>
	DATE1=<earliest collection date (yyyy-mm-dd) of sampled sequences>
	DATE2=<latest collection date (yyyy-mm-dd) of sampled sequences>
	
	mkdir ${run}_${today}
	cd ${run}_${today}
	cp /path/to/consensus/sequences/FASTA ./${run}_${today}.fasta
	gisaidfilt.py ../gisaid_msa/msa_${MSADATE}/msa_${MSADATE}.fasta -n USA/FL -a ${DATE1}_00 -b ${DATE2}+00  -o florida_${today}.fasta 
	```

14. Combine in-house sequences with GISAID Floridian sequences:
	
	```
	cat ${run}_${today}.fasta florida_${today}.fasta > ${run}_florida_${today}.fasta
	```

15. Meanwhile, remove Floridian sequences from GISAID database (flaco_blast now in nextflow), which results in a new "clean.fa" gisaid database:	
16. 
	```
	mkdir ../flaco_blast;
	cd ../flaco_blast;
	flaco_blast.sh makedb ../gisaid_msa/msa_${MSADATE}/msa_${MSADATE}.fasta USA/FL
	```
	
16. FLACOBLAST combined sequences against new GISAID global database:
	
	```
	cd ../${run}_${today}
	seq_cleaner.py -g ${run}_florida_${today}.fasta > ${run}_florida_${today}_gapstripped.fa
	cd ../flaco_blast
	flaco_blast.sh run ../${run}_${today}/${run}_florida_${today}_gapstripped.fa msa_${MSADATE}.clean.fa
	``` 
	**Note**: Used to need to be run in the same folder as the folder containing the dbs from the previous step, but may no longer be the case (will ask)

18. Combine target sequences with original query sequences
	
	```
	cd ../${run}_${today}
	cat ${run}_florida_${today}_gapstripped.fasta ../flaco-blast/msa_${MSADATE}.clean.out.fa > ${run}_florida_gisaid_${today}_gapstripped.fa
	```
	
1. Now we can align sequences using viralmsa:

	```
	ViralMSA.py -s ${run}_florida_gisaid_${today}_gapstripped.fa -t 4 -e 	brittany.rife@ufl.edu -o ${run}_florida_gisaid_${today}_aln -r ../cov_reference/*fasta
	```
	 
2. Check for extensive gaps in sequences (viralmsa can sometimes result in extensive gaps if number of ambiguous sites is too high). This step does not necessarily have to be automated if earlier step filters out low-quality sequences.
	
	```
	grep -o "-" ${run}_florida_${today}_aln/${run}_florida_${today}.aln | wc -l | awk '$1>300{c++} END{print c+0}'
	```

3. Check for gaps in reference sequence (shouldn't be there). If found, need to either use mafft or report in metadata file:
	
	```
	head -n 2 ${run}_florida_${today}_aln/${run}_florida_${today}.aln | tail -n 1 | grep -o "-" | wc -l
	```
	
	(If mafft necessary):
	
	```
	mafft --thread -1 ${run}_florida_${today}.fasta > ${run}_florida_${today}.aln	mafft --retree 3 --maxiterate 10 --thread -1 --nomemsave --op 10 seqsCausingInsertionsInRef.fasta > seqsCausingInsertionsInRef_aligned.fasta
	mafft --thread 1 --quiet --keeplength --add sequencesNotCausingInsertionsInRef.fa seqsCausingInsertionsInRef_aligned.fasta > ${run}_florida_${today}.aln
	```

4. Mask uncertain sites

	```
	python mask_aln_using_vcf.py -i ${run}_florida_${today}_aln/${run}_florida_${today}.aln -o ${run}_florida_${today}_masked.aln -v ../bin/problematic_sites_sarsCov2.vcf
	```
	
5. Change first sequence header to just MN908947.3 (spaces in original name wreak havoc downstream):
	
	```
	REF=">MN908947.3"
	sed -i "1s/.*/$REF/" ${run}_florida_${today}_masked.aln
	```
	
		
# INITIAL TREE RECONSTRUCTION AND DYNAMITE CLUSTER IDENTIFICATION ####################################

1.  Create ML tree using IQ-TREE job script (which will automatically run on a fasta file):
	
	```
	iqtree2 -s ${run}_florida_${today}_masked.aln -bb 1000 -nt AUTO
	```
	
2. Create metadata.txt file containing minimal information used for DYNAMITE:
	
	```
	ml R
	Rscript metadata.R
	```
3. run DYNAMITE to identify clusters using tree and metadata (will output individual trees and fasta files for clusters and background):
	
	```
	sbatch dynamite.sh


# SUBTREE/ALN PROCESSING OF INITIAL TREE USING USHER #######################################################
	


8. Add reference sequence again to each of the newly generated fasta files:
	
	```
	for i in ./*.fasta; do cat ../cov_reference/cov_reference.fasta >> ${i}; done
	```
	
1. Copy FaToVcf to working directory (required for UShER):
		
	```
	rsync -aP rsync://hgdownload.soe.ucsc.edu/genome/admin/exe/linux.x86_64/faToVcf .
	chmod +x faToVcf
	```
				
2. Transform fasta sub-alignments to vcf
	
	```
	for i in ./*.fasta; do ./faToVcf -ref=MN908947.3 ${i} ${i%.fasta}.vcf; done
	```
	
3. Perform subtree pre-processing (need to modify to only select cluster and background trees and to exclude original, full tree) 
	
	```
	ml usher
	for i in ./*${today}.tree; do usher -v ${i%.tree}.vcf -t ${i} -T 4 -o ${i%.tree}.pb; done &
	```




# Wash, rinse, repeat
10. Update (assuming existing) dates.txt file with new date and assign current ("today") and previous ("previous") date variables:
	
	```
	echo $(date +"%Y-%m-%d") >> dates.txt
	today=$(tail -n 1 dates.txt)
	previous=$(sed 'x;$!d' <dates.txt)
	```

11. Make new directory for new round of samples:
	
	```
	mkdir ./<Run#>_${today}
	cd ./<Run#>_${today}
	```




19. Run pangolin on newly combined sequences (updated daily because known to change!)
	
	```
	ml pangolin 	
	pangolin <Run#>_florida_gisaid_${today}_gapstripped.fa --outfile <Run#>_florida_gisaid_${today}_lineages.csv
	```	
# ALIGNMENT OF NEW SEQUENCE DATA ########################################################

1. Align sequences using viralmsa:
	
	```
	sbatch ../viralmsa.sh
	```

If mafft necessary, see above.

2. Mask uncertain sites

	```
	ml python
	python ../mask_aln_using_vcf.py -i <Run#>_florida_gisaid_${today}/*.aln -o <Run#>_florida_gisaid_${today}_masked.aln -v ../problematic_sites_sarsCov2.vcf
	```
	
3. Change first sequence header to just "MN908947.3" (spaces in original name wreak havoc downstream):
	
	```
	var=">MN908947.3"
	sed -i "1s/.*/$var/" *.aln 
	```
	 
		
# SUBTREE/ALN PROCESSING USING USHER #######################################################
		
	
0.5 If full tree not in mutation-annotated form, then process using UsHER:

```
./faToVcf -ref=MN908947.3 updated_${today}_masked.aln updated_${today}_masked.vcf
ml usher
usher -v updated_${today}_masked.vcf -t updated_${previous}_reop.tree -T 4 -u -o updated_${today}.pb

```

1. Transform fasta to vcf
	
	```
	../faToVcf -ref=MN908947.3 ./<Run#>_florida_gisaid_${today}_masked.aln <Run#>_florida_gisaid_${today}_masked.vcf	&
	```		

2. Preliminarily add samples to all previous trees using usher (creating new pb files for today's date):
	
	```
	cd ../
	ml usher
	for i in $(ls | grep "${previous}.pb"); do usher -i ${i} -v <Run#>_florida_gisaid_${today}/<Run#>_florida_gisaid_${today}_masked.vcf -o -p ${i%_${previous}.pb}_${today}_prelim.pb; done	
	```

3. Compute parsimony score for new samples (vcf) assigned to each annotated tree (.pb files):
	
	```
	for i in $(ls | grep "${previous}.pb"); do 
	mkdir ${i%_${previous}.pb}_${today};
	usher -i ${i} -v <Run#>_florida_gisaid_${today}/<Run#>_florida_gisaid_${today}_masked.vcf -p -d ${i%_${previous}.pb}_${today};
	done
	```
	
4. Previous step will generate parsimony scores in tsv files, which we need to rename and copy to single folder for R analysis:
	
	```
	mkdir BPS_${today}
	for i in ./*_${today}; do mv ${i}/parsimony-scores.tsv BPS_${today}/${i}_parsimony-scores.tsv; done
	```
	
5. Use R script to evaluate parsimony scores from tsv files and output new fasta files for sequences needing to be placed on trees:
	
	```
	cd BPS_${today}
	ml R
	Rscript ../branch_support_eval.R
	```	


6. Add reference sequence again to each of the newly generated fasta files:
	
	```
	for i in ./*.fasta; do cat ../cov_reference/cov_reference.fasta >> ${i}; done
	```

7. Place updated fasta in new folders:
	
	```
	for i in ./*updated.fasta; do mkdir .${i%.fasta}; cp ${i} .${i%.fasta}; done;
	cd ../
	``
	
8. Replace old folders in source.txt file with updated folders:
	
	```
	find . -type d -name "*${today}_updated" > source.txt
	```

	
9. Convert to vcf as above
	
	 ```
	for i in $(cat source.txt); do faTovcf -ref=MN908947.3 ${i}/*.fasta ${i}/*.vcf; done	
	```

10. Officially add samples to all old trees using usher (creating new pb files for today):
	
	```
	ml usher
	for i in $(cat source.txt); do usher -i ${i%_${today}*}_${previous}.pb -v ${i}/*.vcf -u -d -o ${i%_updated}.pb; done
	```	
	
11. Not all trees (.pb) are going to be updated with sequences, so need to find all trees without today's date and make a copy of previous pb with today's date:
	
	```
	for i in ./*${previous}.pb; do mv -vn ${i} ${i%_${previous}*}_${today}.pb; done
	```

12. Move new unannotated tree (.nh) files generated from the previous step into one folder (as well as original combined sample fasta) for characterization:
	
	```
	mkdir "trees_${today}";
	for i in $(cat source.txt); do cp ${i}/*final-tree.nh trees_${today}/${i}.tree; done;
	cp ?????? trees_${today}
	```
	
13. Run modified DYNAMITE to search for (and prune) clusters within background tree: *Remind Brittany to modify dynamite to make cluster names in sequential order and if cluster already exists in other folder, add new digit place*
	
	```
	cd ./background_${today}
	sbatch ../dynamite.sh
	```	

14. Run R script to characterize added sequences (add new folder to save discard results):	*This needs to be modified so that when discard tree exists, it will create just a fasta that will need to be processed*
	
	```
	mkdir "./discard_${today}"
	cd trees_${today}
	ml R
	Rscript ../fitness_calc.R
	```
	
15. If discard tree not present, create one. If so, process fasta from previous step and add to existing annotated tree, and also update source.txt file
	
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

16. If first condition met above, then need to pre-process tree and fasta (but somehow need to make sure iqtree is finished:
	    
	
	```
	cat ./cov_reference/cov_reference.fasta >> ./discard_${today}_updated/*.fasta
	ml python
	python Fasta2UShER.py -inpath discard_${today}_updated -output discard_${today}.vcf -reference ./cov_reference/cov_reference.fasta	
	ml usher
	usher -v discard_${today}.vcf -i ./discard_${previous}.pb -T 4 -c -u -d ./discard_${today}_updated -o discard_${today}.pb
	```

17. While discard tree being reconstructed, place pruned background sequences onto previous annotated tree:
	
	```
	cd ../background_${today}_updated
	rm -v !(*.fasta)
	cd ../
	ml python
	python Fasta2UShER.py -inpath ./background_${today}_updated -output ./background_${today}.vcf -reference ./cov_reference/cov_reference.fasta
	ml usher
	usher -v background_${today}.vcf -i ./background_${previous}.pb -T 4 -c -u -d ./background_${today}_updated -o background_${today}.pb
	```
	

18. Concatenate previous fasta with updated fasta for tree optimization :
	
	```
	for i in $(cat source.txt); do 
tail -n +2 ${i}/*.fasta >> ../full.fasta
	done
	```
19. Remove duplicated sequences:
	
	```
	for i in $(cat source.txt); cat ${i}/full.fasta | seqkit rmdup -n -o ${i}/clean.aln
	```

20. Force bifurcating tree in R:
	
	```
	ml R
	for i in $(cat source.txt); do Rscript ./bifurcate.R -t ${i}/uncondensed-final-tree.nh; done
	```

 
21. Re-optimize trees (only those that were updated) using FastTree: 
	
	```
	module load fasttree/2.1.7
	for i in $(cat source.txt); do FastTreeMP -nt -gamma -nni 0 -spr 2 -sprlength 1000 -boot 100 -log fasttree.log -intree ${i}/*.nh.tree ${i}/clean.aln > ${i}/reop.tree; done
	```
	

22. Re-annotate updated fasttree trees (replace old ones):
	
	```
	ml usher
	for i in ${cat source.txt); do
  	usher -v ${i%_updated}.vcf -t ${i}/reop.tree -T 4 -c -u -d ${i} -o ${i%_updated}.pb;
  	done
  	```

23. Test background and discard trees for clusters (separate R script using actual branchwise algorithm). R script will write new results for discard pile to a new folder, but will write results from background to same folder.
	IF no clusters are found, results are not written, but downstream steps require that the new discard folder exist, so need to copy over.

	```
	mkdir ./discard_${today}_updated
	cd ./discard_${today}
	Rscript ../branchwise.R -t *.treefile -p ../discard_${today}_updated
	
	cd ../background_${today}_updated
	Rscript ../branchwise.R -t reop.tree -p ./
	cd ../
	```
		
24. Process fasta and tree for background (will require determining if new files added as a result of R script in step above)

	```
	cat ./cov_reference/cov_reference.fasta >> ./background_${today}_updated/*.fasta
	find . -type f -name "background_${today}_updated.tree" -empty -exec cp ./background_${today}/reop.tree ./background_${today}_updated.tree
	ml python
	python Fasta2UShER.py -inpath background_${today}_updated -output background_${today}.vcf -reference ./cov_reference/cov_reference.fasta	
	ml usher
	usher -v background_${today}.vcf -t ./background_${today}_updated.tree -T 4 -c -u -d ./background_${today}_updated -o background_${today}.pb
	```
	
25. Add reference sequence to new cluster fastas and place in new directories:
	
	```
	for i in ./*.fasta; do  cat ./cov_reference/cov_reference.fasta >> ${i}; cp ${i} .${i%.fasta}; done
	```
	
26. Update source.txt file with new clusters and process fastas:
	
	```
	for i in ./*.fasta; do mkdir .${i%.fasta}; echo ${i%.fasta} > source.txt; done	
	ml python
	for i in $(cat source.txt); do python Fasta2UShER.py -inpath ${i} -output ${i}.vcf -reference ./cov_reference/cov_reference.fasta; done
	```
	
27. Perform subtree pre-processing for new clusters:
	
	```
	ml usher	
	for i in $(cat source.txt); do usher -v ${i}.vcf -t ${i}.tree -T 4 -c -u -d ./${i} -o ${i}.pb; done
	```





			
			







