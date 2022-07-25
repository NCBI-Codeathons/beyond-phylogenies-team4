# Workflow

## Pre-requisite installs:
	* python3 (standard)
	* R (https://cloud.r-project.org/)
	* viralmsa (https://github.com/niemasd/ViralMSA)	* pangolin (https://cov-lineages.org/resources/pangolin/installation.html)
	* iqtree2 (https://apolo-docs.readthedocs.io/en/latest/software/applications/iqtree/2.1.2/index.html)
	* FastTreeMP (http://www.microbesonline.org/fasttree/#Install)
	* mafft (optional)	

## Initalize environment, alignment, and tree

1. Initiate FLACO environment (or simply add bin to PATH)

	```
	git clone https://github.com/salemilab/flaco.git;
	source /flaco/bin/FLACO 
	```
	
2. Download most up-to-date gisaid alignment and place in "gisaid_msa" directory. Below no longer works for some reason, so have to move in manually.

	```
	MSADATE=<Date (mmdd) of available GISAID msa>;
	mkdir gisaid_msa;
	wget --user <username> --password <password> https://www.epicov.org/epi3/msa_${MSADATE}.tar.xz;
	tar -xvf msa*;
	cd ../;
	```

3. Create a file with dates for later referencing and a variable ("today") from the dates file to help with naming:
	
	```
	echo $(date +"%Y-%m-%d") > dates.txt;
	today=$( tail -n 1 dates.txt );
	```

	
4. Download most recent masked sites vcf (needs to be automated so that new file replaces old file instead of being renamed):

	```
	rm problematic_sites_sarsCov2.vcf;
	wget - O ../bin/problematic_sites_sarsCov2.vcf https://raw.githubusercontent.com/W-L/ProblematicSites_SARS-CoV2/master/problematic_sites_sarsCov2.vcf &
	```


5. Meanwhile, create directory for initial sequencing run and extract sequences from desired region from database that correspond to the relevant time frame:
	
	```
	run=<Sequencing run number/name>;
	DATE1=<earliest collection date (yyyy-mm-dd) of sampled sequences>;
	DATE2=<latest collection date (yyyy-mm-dd) of sampled sequences>;
	
	mkdir ${run}_${today};
	cd ${run}_${today};
	cp /path/to/consensus/sequences/FASTA ./${run}_${today}.fasta;
	gisaidfilt.py ../gisaid_msa/msa_${MSADATE}/msa_${MSADATE}.fasta -n USA/FL -a ${DATE1}_00 -b ${DATE2}+00  -o florida_${today}.fasta; 
	```

6. Combine in-house sequences with GISAID Floridian sequences:
	
	```
	cat ${run}_${today}.fasta florida_${today}.fasta > ${run}_florida_${today}.fasta;
	```

7. Remove Floridian sequences from GISAID database, which results in a new "clean.fa" gisaid database:	
 
	```
	mkdir ../flaco_blast;
	cd ../flaco_blast;
	flaco_blast.sh makedb ../gisaid_msa/msa_${MSADATE}/msa_${MSADATE}.fasta USA/FL;
	```
	
8. (FLACO)BLAST combined sequences against new GISAID global database:
	
	```
	cd ../${run}_${today}
	seq_cleaner.py -g ${run}_florida_${today}.fasta > ${run}_florida_${today}_gapstripped.fa
	cd ../flaco_blast
	flaco_blast.sh run ../${run}_${today}/${run}_florida_${today}_gapstripped.fa msa_${MSADATE}.clean.fa
	```

9. Combine target sequences with original query sequences
	
	```
	cd ../${run}_${today}
	cat ${run}_florida_${today}_gapstripped.fasta ../flaco-blast/msa_${MSADATE}.clean.out.fa > ${run}_florida_gisaid_${today}_gapstripped.fa
	```
	
10. Run pangolin on newly combined sequences (updated daily because known to change!)
	
	```
	pangolin ${run}_florida_gisaid_${today}_gapstripped.fasta --outfile ${run}_florida_gisaid_${today}_lineages.csv	
	```
	
11. Multiple sequence alignment using viralmsa:

	```
	ViralMSA.py -s ${run}_florida_gisaid_${today}_gapstripped.fa -t 4 -e 	brittany.rife@ufl.edu -o ${run}_florida_gisaid_${today}_aln -r ../cov_reference/*fasta
	```
	
	**If no step prior to workflow exclude low-quality sequences (specifically high N content, proceed with Steps 12-14, else proceed to Step 15.**

12. Check for extensive gaps in sequences (viralmsa can sometimes result in extensive gaps if number of ambiguous sites is too high). This step does not necessarily have to be automated if earlier step filters out low-quality sequences.
	
	```
	grep -o "-" ${run}_florida_gisaid_${today}_aln/${run}_florida_gisaid_${today}.aln | wc -l | awk '$1>300{c++} END{print c+0}'
	```

13. Check for gaps in reference sequence (shouldn't be there). If found, need to either use mafft or report in metadata file:
	
	```
	head -n 2 ${run}_florida_gisaid_${today}_aln/${run}_florida_gisaid_${today}.aln | tail -n 1 | grep -o "-" | wc -l
	```
	
14. If mafft necessary, run the following:
	
	```
	mafft --thread -1 ${run}_florida_gisaid_${today}.fasta > ${run}_florida_gisaid_${today}_aln/${run}_florida_gisaid_${today}.aln	mafft --retree 3 --maxiterate 10 --thread -1 --nomemsave --op 10 seqsCausingInsertionsInRef.fasta > seqsCausingInsertionsInRef_aligned.fasta
	mafft --thread 1 --quiet --keeplength --add sequencesNotCausingInsertionsInRef.fa seqsCausingInsertionsInRef_aligned.fasta > ${run}_florida_gisaid_${today}.aln
	```

	```

15. Mask uncertain sites

	```
	python mask_aln_using_vcf.py -i ${run}_florida_gisaid_${today}_aln/${run}_florida_gisaid_${today}.aln -o ${run}_florida_gisaid_${today}_masked.aln -v ../bin/problematic_sites_sarsCov2.vcf
	```
	
16. Change first sequence header to just MN908947.3 (spaces in original name wreak havoc downstream):
	
	```
	REF=">MN908947.3"
	sed -i "1s/.*/$REF/" ${run}_florida_gisaid_${today}_masked.aln
	```
	
17.  Create maximum likelihood (ML) starting tree using IQ-TREE job script (which will automatically run on a fasta file):
	
	```
	iqtree2 -s ${run}_florida_gisaid_${today}_masked.aln -bb 1000 -nt AUTO &
	```
	
	
18. Create metadata.txt file containing minimal information used for DYNAMITE/PHYLOPART (included script written for full tree processing, rather than cluster processing so will need to be modified):
	
	```
	Rscript metadata.R --tree *treefile --metadata ${run}_metadata.tab --columnName SampleName --lineages  ${run}_florida_gisaid_${today}_lineages.csv
	```

19. Run DYNAMITE to identify clusters using tree and metadata (will output individual trees and fasta files for clusters and background):
	
	```
	sbatch dynamite.sh
	```


## SUBTREE/ALN PROCESSING OF INITIAL TREE(S) USING USHER #######################################################
	


20. Add reference sequence again to each of the newly generated fasta files:
	
	```
	for i in ./*.fasta; do cat ../cov_reference/cov_reference.fasta >> ${i}; done
	```
21. Copy FaToVcf to working directory (required for UShER):
		
	```
	rsync -aP rsync://hgdownload.soe.ucsc.edu/genome/admin/exe/linux.x86_64/faToVcf .
	chmod +x faToVcf
	```
				
22. Transform fasta sub-alignments to vcf
	
	```
	for i in ./*.fasta; do ./faToVcf -ref=MN908947.3 ${i} ${i%.fasta}.vcf; done
	```
	
23. Perform subtree pre-processing (need to modify to automatically only select cluster and background trees and to exclude original, full tree) 
	
	```
	for i in ./*${today}.tree; do usher -v ${i%.tree}.vcf -t ${i} -T 4 -o ${i%.tree}.pb; done &
	```


## Wash, rinse, repeat
24. Update dates.txt file with new date and assign current ("today") and previous ("previous") date variables and follow steps above:
	
	```
	echo $(date +"%Y-%m-%d") >> dates.txt
	today=$(tail -n 1 dates.txt)
	previous=$(sed 'x;$!d' <dates.txt)
	```

**May want to add additional details here since will be slightly different (e.g., no iqtree)**

## SUBTREE/ALN PROCESSING USING USHER #######################################################


25. Preliminarily add samples to all previous trees using usher (creating new pb files for today's date):
	
	```
	for i in $(ls | grep "${previous}.pb"); 
	do usher -i ${i} -v ${run}_${today}/${run}_florida_gisaid_${today}_masked.vcf -o -p ${i%_${previous}.pb}_${today}_prelim.pb; done	
	```

26. Compute parsimony score for new samples (vcf) assigned to each annotated tree (.pb files):
	
	```
	for i in $(ls | grep "${previous}.pb"); do 
	mkdir ${i%_${previous}.pb}_${today};
	usher -i ${i} -v {run}_${today}/${run}_florida_gisaid_${today}_masked.vcf -p -d ${i%_${previous}.pb}_${today};
	done
	```
	
27. Previous step will generate parsimony scores in tsv files, which we need to rename and copy to single folder for R analysis:
	
	```
	mkdir BPS_${today}
	for i in ./*_${today}; do mv ${i}/parsimony-scores.tsv BPS_${today}/${i}_parsimony-scores.tsv; done
	```
	
28. Use R script to evaluate parsimony scores from tsv files and output new fasta files for sequences needing to be placed on trees (currently not working owing to change in usher output):
	
	```
	cd BPS_${today}
	ml R
	Rscript ../branch_support_eval.R
	```	


29. Add reference sequence again to each of the newly generated fasta files (from previous step):
	
	```
	for i in ./*.fasta; do cat ../cov_reference/cov_reference.fasta >> ${i}; done
	```

30. Place updated fasta in new folders:
	
	```
	for i in ./*updated.fasta; do mkdir .${i%.fasta}; cp ${i} .${i%.fasta}; done;
	cd ../
	``
	
31. Replace old folders in source.txt file with updated folders (noticed old source.txt file generation step is missing...):
	
	```
	find . -type d -name "*${today}_updated" > source.txt
	```

	
32. Convert to vcf as above
	
	 ```
	for i in $(cat source.txt); do faTovcf -ref=MN908947.3 ${i}/*.fasta ${i}/*.vcf; done	
	```

33. Officially add samples to all old trees using usher (creating new pb files for today):
	
	```
	for i in $(cat source.txt); 
	do usher -i ${i%_${today}*}_${previous}.pb -v ${i}/*.vcf -u -d -o ${i%_updated}.pb; done
	```	
	
34. Not all trees (.pb) are going to be updated with sequences, so need to find all trees without today's date and make a copy of previous pb with today's date:
	
	```
	for i in ./*${previous}.pb; 
	do mv -vn ${i} ${i%_${previous}*}_${today}.pb; done
	```

35. Move new unannotated tree (.nh) files generated from the previous step into one folder (as well as original combined sample fasta) for characterization:
	
	```
	mkdir "trees_${today}";
	for i in $(cat source.txt); do cp ${i}/*final-tree.nh trees_${today}/${i}.tree; done;
	cp ?????? trees_${today}
	```
	
36. Run modified DYNAMITE to search for (and prune) clusters within background tree: *Remind Brittany to modify dynamite to make cluster names in sequential order and if cluster already exists in other folder, add new digit place*
	
	```
	cd ./background_${today}
	sbatch ../dynamite.sh
	```	

37. Run R script to characterize added sequences (add new folder to save discard results):	*This needs to be modified so that when discard tree exists, it will create just a fasta that will need to be processed*
	
	```
	mkdir "./discard_${today}"
	cd trees_${today}
	Rscript ../fitness_calc.R
	```
	
38. If discard tree not present, create one. If so, process fasta from previous step and add to existing annotated tree, and also update source.txt file (sh scripts need to be included or the following code modified by Brittany):
	
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

39. If first condition met above, then need to pre-process tree and fasta (but somehow need to make sure iqtree is finished:
	    
	
	```
	cat ./cov_reference/cov_reference.fasta >> ./discard_${today}_updated/*.fasta
	python Fasta2UShER.py -inpath discard_${today}_updated -output discard_${today}.vcf -reference ./cov_reference/cov_reference.fasta	usher -v discard_${today}.vcf -i ./discard_${previous}.pb -T 4 -c -u -d ./discard_${today}_updated -o discard_${today}.pb
	```

40. While discard tree being reconstructed, place pruned background sequences onto previous annotated tree:
	
	```
	cd ../background_${today}_updated
	rm -v !(*.fasta)
	cd ../
	python Fasta2UShER.py -inpath ./background_${today}_updated -output ./background_${today}.vcf -reference ./cov_reference/cov_reference.fasta
	ml usher
	usher -v background_${today}.vcf -i ./background_${previous}.pb -T 4 -c -u -d ./background_${today}_updated -o background_${today}.pb
	```
	

41. Concatenate previous fasta with updated fasta for tree optimization :
	
	```
	for i in $(cat source.txt); do 
tail -n +2 ${i}/*.fasta >> ../full.fasta
	done
	```
42. Raxml is very picky about duplicate sequences. Remove duplicated sequences:
	
	```
	for i in $(cat source.txt); cat ${i}/full.fasta | seqkit rmdup -n -o ${i}/clean.aln
	```

43. Force bifurcating tree in R:
	
	```
	for i in $(cat source.txt); do Rscript ./bifurcate.R -t ${i}/uncondensed-final-tree.nh; done
	```

 
44. Re-optimize trees (only those that were updated) using FastTree: 
	
	```
	for i in $(cat source.txt); do FastTreeMP -nt -gamma -nni 0 -spr 2 -sprlength 1000 -boot 100 -log fasttree.log -intree ${i}/*.nh.tree ${i}/clean.aln > ${i}/reop.tree; done
	```
	

45. Re-annotate updated fasttree trees (replace old ones):
	
	```
	for i in ${cat source.txt); do
  	usher -v ${i%_updated}.vcf -t ${i}/reop.tree -T 4 -c -u -d ${i} -o ${i%_updated}.pb;
  	done
  	```

46. Test background and discard trees for clusters (separate R script using actual branchwise algorithm). R script will write new results for discard pile to a new folder, but will write results from background to same folder.
	IF no clusters are found, results are not written, but downstream steps require that the new discard folder exist, so need to copy over.

	```
	mkdir ./discard_${today}_updated
	cd ./discard_${today}
	Rscript ../branchwise.R -t *.treefile -p ../discard_${today}_updated
	
	cd ../background_${today}_updated
	Rscript ../branchwise.R -t reop.tree -p ./
	cd ../
	```
		
47. Process fasta and tree for background (will require determining if new files added as a result of R script in step above)

	```
	cat ./cov_reference/cov_reference.fasta >> ./background_${today}_updated/*.fasta
	find . -type f -name "background_${today}_updated.tree" -empty -exec cp ./background_${today}/reop.tree ./background_${today}_updated.tree
	python Fasta2UShER.py -inpath background_${today}_updated -output background_${today}.vcf -reference ./cov_reference/cov_reference.fasta	
	usher -v background_${today}.vcf -t ./background_${today}_updated.tree -T 4 -c -u -d ./background_${today}_updated -o background_${today}.pb
	```
	
48. Add reference sequence to new cluster fastas and place in new directories:
	
	```
	for i in ./*.fasta; do  cat ./cov_reference/cov_reference.fasta >> ${i}; cp ${i} .${i%.fasta}; done
	```
	
49. Update source.txt file with new clusters and process fastas:
	
	```
	for i in ./*.fasta; do mkdir .${i%.fasta}; echo ${i%.fasta} > source.txt; done	
	ml python
	for i in $(cat source.txt); do python Fasta2UShER.py -inpath ${i} -output ${i}.vcf -reference ./cov_reference/cov_reference.fasta; done
	```
	
50. Perform subtree pre-processing for new clusters:
	
	```
	ml usher	
	for i in $(cat source.txt); do usher -v ${i}.vcf -t ${i}.tree -T 4 -c -u -d ./${i} -o ${i}.pb; done
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

7. Create R list of lineage-specific alignments for spike (not needed now)

	```
	makeRefSpikeList.sh -r ref_lineages_${today}_masked.aln	
	```		







