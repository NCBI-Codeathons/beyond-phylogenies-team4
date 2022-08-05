# PhyloPRIME (Phylogeny-based Pipeline for RegIonal Molecular Epidemiology)


## Table of Contents
- [Introduction](#Introduction)
- [Methodology](#Methodology)
- [Results](#Results)
- [Future Work](#Future-Work)
- [Dependencies](#Dependencies)
- [Team](#Team)
- [Acknowledgment](#Acknowledgment)
- [References](#References)
- [License](#License)

## Introduction
The workflow provided in this repository is meant to provide a streamlined approach for regional molecular epidemiology analysis of continuously incoming SARS-CoV-2 sequencing data. Built-in functions for retrieving relevant sequence data from the global sequence repository [GISAID](https://gisaid.org) and pruning uninformative sequences help the user to draw meaningful conclusions in real-time from extremely large datasets. Morever, a [clustering approach][phylopart] that utilizes available metadata (including lineage information derived from [Pangolin][pangolin]) and a variety of built-in tree statistics allows for the user to focus on groups of sequences that may comprise specific risk factors that contribute to quantifiably increased or reduced transmission.

## Methodology

In-house sequences can be combined with additional sequences deposited into the region using an initial function as part of the [flaco_blast worflow][flacoblast], which searches an updated GISAID multiple sequence alignment for sequences belonging to the user-specified region based on name.

Each sequence can then be used to search the remaining GISAID database for epidemiologically relevant sequences - sequences from other regions collected within a certain time windown prior to, and following, collection date of the query sequence based on genetic similarity (best hit chosen using [flaco_blast worflow][flacoblast]).

These sequences are then concatenated, and quality control steps are taken to [mask problematic sites](https://virological.org/t/issues-with-sars-cov-2-sequencing-data/473) and [align sequence data][viralmsa] for tree reconstruction. 

Lineages are called using [Pangolin][pangolin], and the results are merged with the user-specified metadata for downstream interactive plotting an analysis in [R][cran].

[Initial tree reconstruction][iqtree] is performed using maximum likelihood. Following tree reconstruction, a depth-first search algorithm is used to cluster groups of sequences according to median branch length that falls within a user-specified threshold, as described in [Prosperi et al. (2011)][phylopart]. Lower thresholds (e.g., 0.01-0.05) correspond to a small portion of the total branch length distribution in the tree and thus to putative transmission chains, whereas larger thresholds will result in much larger clades within the tree (e.g., COVID-19 lineages).

Clusters/clades are pruned from the "background" tree/population, resulting in distinct subtrees on which the remaining workflow operates.


Additional sequence data are added to individual trees using [rapid phylogenetic placement][usher] so as to avoid time-consuming tree reconstruction for each round of added sequence data. Additional criteria are used, however, to determine if added sequences contribute to the corresponding cluster or are informative with regards to the background population:

- For sequence(s) added to a cluster, the median branch length for the cluster is recalculated, and the added sequence is reported as part of the cluster if this value is still below the user-specified threshold. If not, it is considered to have formed a new sub-cluster.
- For sequence(s) added to the background population, a fitness score is assigned to individual sequences based on their contribution to the 1) genetic diversity for that sequence's corresponding region and 2) sampling time distribution for that sequence's corresponding region, favoring a more uniform distribution. This process is described in more detail in [Marini et al. (2022)][tardis]. Sequences with low fitness scores are discarded (but not entirely deleted).

Trees are continuously [re-optimized][fastree] (branch lengths and placement) so as to prevent donstream propagation of placement errors.


## Results

The output from this workflow include the following:

1. updated tabular file with metadata, lineage, and cluster information
2. individual cluster and background trees
3. tabular file with mutational information for each cluster (new mutations relative to reference strain)
4. Plots for novel mutations in Spike, lineage proportions over time, and metadata column distributions over time for each cluster.

Plots can be viewed interactively in the [R Shiny][rshiny] application, including filtering of clusters based on transmission rate, calculated using the [Oster statistic][oster] and described in the context of respiratory virus transmission cluster dynamics by Sun et al. (2022)[deepdynatree].

## Future Work

Our future work for this pipeline includes:

* Expanding to other viruses without loss of automated, intuitive design
* Branch-site selection analysis
* Visualization of mutations on published protein structures
* Integration of phylogenetic regression for clinical data (e.g., viral load)


## Dependencies
* [python3](https://www.python.org)
* [R v4.2.1](https://cloud.r-project.org/)
* [ViralmMSA](https://github.com/niemasd/ViralMSA)
* [Pangolin](https://cov-lineages.org/resources/pangolin/installation.html)
* [IQ-Tree2](https://apolo-docs.readthedocs.io/en/latest/software/applications/iqtree/2.1.2/index.html)
* [FastTreeMP](http://www.microbesonline.org/fasttree/#Install)
* [mafft (optional)](https://mafft.cbrc.jp/alignment/software/)

**More on specifics for AWS and R libraries to come

## Team
- Aitor Serres Armero, National Human Genome Research Institute
- Anthony Fries, United States Air Force School of Aerospace Medicine
- Brittany Magalis, University of Florida
- Manoj M Wagle, Université Grenoble Alpes
- Marco Salemi, University of Florida
- Selva, MD Anderson Cancer Center
- Toby Koch, Affiliation

## Acknowledgment
We would like to thank the NIH National Library of Medicine/National Center for Biotechnology Information for providing all the required computational resources during the codeathon.

## References
[phylopart]: Prosperi MC, Ciccozzi M, Fanti I, Saladini F, Pecorari M, Borghi V, Di Giambenedetto S, Bruzzone B, Capetti A, Vivarelli A, Rusconi S, Re MC, Gismondo MR, Sighinolfi L, Gray RR, Salemi M, Zazzi M, De Luca A; ARCA collaborative group. A novel methodology for large-scale phylogeny partition. Nat Commun. 2011;2:321. doi: 10.1038/ncomms1325. Epub 2011 May 24. PMID: 21610724; PMCID: PMC6045912.

[pangolin]: O'Toole Á, Pybus OG, Abram ME, Kelly EJ, Rambaut A. Pango lineage designation and assignment using SARS-CoV-2 spike gene nucleotide sequences. BMC Genomics. 2022 Feb 11;23(1):121. doi: 10.1186/s12864-022-08358-2. PMID: 35148677; PMCID: PMC8832810.

[flacoblast]: Giovanetti M, Cella E, Benedetti F, Rife Magalis B, Fonseca V, Fabris S, Campisi G, Ciccozzi A, Angeletti S, Borsetti A, Tambone V, Sagnelli C, Pascarella S, Riva A, Ceccarelli G, Marcello A, Azarian T, Wilkinson E, de Oliveira T, Alcantara LCJ, Cauda R, Caruso A, Dean NE, Browne C, Lourenco J, Salemi M, Zella D, Ciccozzi M. SARS-CoV-2 shifting transmission dynamics and hidden reservoirs potentially limit efficacy of public health interventions in Italy. Commun Biol. 2021 Apr 21;4(1):489. doi: 10.1038/s42003-021-02025-0. PMID: 33883675; PMCID: PMC8060392.

[viralmsa]: Moshiri N. ViralMSA: massively scalable reference-guided multiple sequence alignment of viral genomes. Bioinformatics. 2021 May 5;37(5):714-716. doi: 10.1093/bioinformatics/btaa743. PMID: 32814953.

[iqtree]: Minh BQ, Schmidt HA, Chernomor O, Schrempf D, Woodhams MD, von Haeseler A, Lanfear R. IQ-TREE 2: New Models and Efficient Methods for Phylogenetic Inference in the Genomic Era. Mol Biol Evol. 2020 May 1;37(5):1530-1534. doi: 10.1093/molbev/msaa015. Erratum in: Mol Biol Evol. 2020 Aug 1;37(8):2461. PMID: 32011700; PMCID: PMC7182206.

[usher]: Turakhia Y, Thornlow B, Hinrichs AS, De Maio N, Gozashti L, Lanfear R, Haussler D, Corbett-Detig R. Ultrafast Sample placement on Existing tRees (UShER) enables real-time phylogenetics for the SARS-CoV-2 pandemic. Nat Genet. 2021 Jun;53(6):809-816. doi: 10.1038/s41588-021-00862-7. Epub 2021 May 10. PMID: 33972780; PMCID: PMC9248294.

[cran]: R Core Team. https://www.r-project.org, 2020.

[tardis]: Marini S, Mavian C, Riva A, Salemi M, Magalis BR. Optimizing viral genome subsampling by genetic diversity and temporal distribution (TARDiS) for phylogenetics. Bioinformatics. 2021 Oct 21;38(3):856–60. doi: 10.1093/bioinformatics/btab725. Epub ahead of print. PMID: 34672334; PMCID: PMC8756195.

[fasttree]: Zhou X, Shen XX, Hittinger CT, Rokas A. Evaluating Fast Maximum Likelihood-Based Phylogenetic Programs Using Empirical Phylogenomic Data Sets. Mol Biol Evol. 2018 Feb 1;35(2):486-503. doi: 10.1093/molbev/msx302. PMID: 29177474; PMCID: PMC5850867.

[rshiny]: R Core Team. https://CRAN.R-project.org/package=shiny, 2022.

[oster]: Oster AM, France AM, Panneer N, Bañez Ocfemia MC, Campbell E, Dasgupta S, Switzer WM, Wertheim JO, Hernandez AL. Identifying Clusters of Recent and Rapid HIV Transmission Through Analysis of Molecular Surveillance Data. J Acquir Immune Defic Syndr. 2018 Dec 15;79(5):543-550. doi: 10.1097/QAI.0000000000001856. PMID: 30222659; PMCID: PMC6231979.

[deepdynatree]: Sun C, Li Y, Marini S, Riva A, Wu DO, Salemi M, Rife Magalis B. Phylogenetic-informed graph deep learning to classify dynamic transmission clusters in infectious disease epidemics. bioRxiv 2022.04.10.487587; doi: https://doi.org/10.1101/2022.04.10.487587

## License
Licensed under MIT License - Copyright (c) 2022 NCBI-Codeathons (Refer to LICENSE file for more details)
