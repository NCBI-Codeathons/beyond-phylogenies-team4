# Workflow for initialization of mutation-annotated cluster phylogenies

## Dependencies
* [python3](https://www.python.org)
* [R v4.2.1](https://cloud.r-project.org/)
* [ViralmMSA](https://github.com/niemasd/ViralMSA)
* [mafft (optional)](https://mafft.cbrc.jp/alignment/software/)
* [Pangolin](https://cov-lineages.org/resources/pangolin/installation.html)
* [IQ-Tree2](https://apolo-docs.readthedocs.io/en/latest/software/applications/iqtree/2.1.2/index.html)
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
	
1. Initiate file with today's date, which will be populated later with dates for each worrkflow use for updating.

	```
	echo $(date +"%Y-%m-%d") > dates.txt;
	today=$( tail -n 1 dates.txt );
	```


2. Pull flaco dependencies from repo (required for sequence cleaning and FLACO BLAST):

	```
	git clone https://github.com/salemilab/flaco.git;
	PATH="$PATH:/flaco/bin";
	```
	
3. Download most recent masked sites vcf (needs to be automated so that new file replaces old file instead of being renamed):

	```
	rm problematic_sites_sarsCov2.vcf;
	wget - O ../bin/problematic_sites_sarsCov2.vcf https://raw.githubusercontent.com/W-L/ProblematicSites_SARS-CoV2/master/problematic_sites_sarsCov2.vcf &
	```


4. Extract sequences from desired region from database that correspond to the relevant time frame:

	```
	gisaidfilt.py ../gisaid_msa/msa_${MSADATE}/msa_${MSADATE}.fasta -n USA/FL -a ${DATE1}_00 -b ${DATE2}+00  -o florida_${today}.fasta;
	```

5. Combine in-house sequences with GISAID Floridian sequences:
	
	```
	cat ${run}_${today}.fasta florida_${today}.fasta > ${run}_florida_${today}.fasta;
	```

6. Remove Floridian sequences from GISAID database, which results in a new "clean.fa" gisaid database:	
 
	```
	mkdir ../flaco_blast;
	cd ../flaco_blast;
	mkdir tmp;
	flaco_blast.sh makedb ../gisaid_msa/msa_${MSADATE}/msa_${MSADATE}.fasta USA/FL;
	```
	
7. Strip sequences of gaps (improves BLAST performance):

	```
	cd ../${run}_${today}
	seq_cleaner.py -g ${run}_florida_${today}.fasta > ${run}_florida_${today}_gapstripped.fa
	```

10. (FLACO)BLAST combined sequences against new GISAID global database:
	
	```
	cd ../flaco_blast
	flaco_blast.sh run ../${run}_${today}/${run}_florida_${today}_gapstripped.fa msa_${MSADATE}.clean.fa
	```

11. Combine target sequences with original query sequences
	
	```
	cd ../${run}_${today}
	cat ${run}_florida_${today}_gapstripped.fa ${run}_florida_${today}_gapstripped.out.fa > ${run}_florida_gisaid_${today}_gapstripped.fa
	```
	
12. Run pangolin on newly combined sequences
	
	```
	pangolin ${run}_florida_gisaid_${today}_gapstripped.fasta --outfile ${run}_florida_gisaid_${today}_lineages.csv &	
	```
	
13. In the meantime, perform multiple sequence alignment using viralmsa:

	```
	ViralMSA.py -s ${run}_florida_gisaid_${today}_gapstripped.fa -t 4 -e 	brittany.rife@ufl.edu -o ${run}_florida_gisaid_${today}_aln -r ../cov_reference/*fasta
	```
	

14. **Steps 14-15 need to be automated so that if true, mafft is used; else, viralmsa alignment is kept**. Check for extensive gaps in sequences (viralmsa can sometimes result in extensive gaps if number of ambiguous sites is too high). This step does not necessarily have to be automated if earlier step filters out low-quality sequences.
	
	```
	grep -o "-" ${run}_florida_gisaid_${today}_aln/${run}_florida_gisaid_${today}.aln | wc -l | awk '$1>300{c++} END{print c+0}'
	```

15. Check for gaps in reference sequence (shouldn't be there). If found, need to either use mafft or report in metadata file:
	
	```
	head -n 2 ${run}_florida_gisaid_${today}_aln/${run}_florida_gisaid_${today}.aln | tail -n 1 | grep -o "-" | wc -l
	```
	
16. If mafft necessary, run the following:
	
	```
	mafft --thread -1 ${run}_florida_gisaid_${today}.fasta > ${run}_florida_gisaid_${today}_aln/${run}_florida_gisaid_${today}.aln	mafft --retree 3 --maxiterate 10 --thread -1 --nomemsave --op 10 seqsCausingInsertionsInRef.fasta > seqsCausingInsertionsInRef_aligned.fasta
	mafft --thread 1 --quiet --keeplength --add sequencesNotCausingInsertionsInRef.fa seqsCausingInsertionsInRef_aligned.fasta > ${run}_florida_gisaid_${today}.aln
	```

	```

17. Mask uncertain sites

	```
	python3 mask_aln_using_vcf.py -i ${run}_florida_gisaid_${today}_aln/${run}_florida_gisaid_${today}.aln -o ${run}_florida_gisaid_${today}_masked.aln -v problematic_sites_sarsCov2.vcf
	```
		
19.  Create maximum likelihood (ML) starting tree using IQ-TREE job script (which will automatically run on a fasta file):
	
	```
	iqtree2 -s ${run}_florida_gisaid_${today}_masked.aln -bb 1000 -nt AUTO &
	```
	
	
20. Merge metadata with pangolin file containing minimal information used for DYNAMITE/PHYLOPART (**May need to modify so that final column names are ID and DATE and first two columns**):
	
	```
	Rscript metadata.R  --metadata ${run}_metadata.tab --columnName SampleName --lineages  ${run}_florida_gisaid_${today}_lineages.csv
	```

21. Run DYNAMITE to identify clusters using tree and metadata (will output individual trees and fasta files for clusters and background):
	
	```
	tree=$(ls *treefile)
	Rscript dynamite.R -f ${run}_florida_gisaid_${today}_masked.aln -t ${tree} -m ../updated_metadata_${today}.tab -c c -q N 
	```
22. Add reference sequence again to each of the newly generated fasta files (original fasta is split according to cluster):
	
	```
	for i in ./*.fasta; do cat ../cov_reference/cov_reference.fasta >> ${i}; done
	```
	
23. Copy FaToVcf to working directory (required for UShER):
		
	```
	rsync -aP rsync://hgdownload.soe.ucsc.edu/genome/admin/exe/linux.x86_64/faToVcf .
	chmod +x faToVcf
	```
				
24. Transform fasta sub-alignments (now distinct for each cluster) to vcf format (for UShER processing of mutation-annotated trees):
	
	```
	for i in ./*.fasta; do ./faToVcf -ref=MN908947.3 ${i} ${i%.fasta}.vcf; done
	```
	
25. Perform subtree pre-processing with UShER: 
	
	```
	for i in ./*${today}.tree; do usher -v ${i%.tree}.vcf -t ${i} -T 4 -o ${i%.tree}.pb; done
	```

