# flaco_blast

flaco_blast is a pipeline to perform BLAST analysis of a set of sequences against the GISAID database. It consists of a shell script that provides
six different commands. Normally, only the *makedb* and *run* commands will be needed.

## Prerequisites
This program requires [Nextflow](https://nextflow.io/) and BLAST. Please ensure they are properly installed and in path.

The program assumes that the input files are in FASTA format, and that sequence headers respect the GISAID format. Specifically, the header
should contain at least three fields separated by | characters: sequence name, sequence identifier, date in the format YYYY-MM-DD. For example:

```
>hCoV-19/Italy/LAZ-INMI1-isl/2020|EPI_ISL_410545|2020-01-29
```

## Execution

Basic usage is:

```
flaco_blast.sh command [options...]
```

where *command* is one of makedb, split, blast, parse, extract, or run. In normal situations, it should be enough to run the *makedb* command (that creates
the BLAST database) followed by the *run* command (which runs the other four commands in the correct order). For example, assuming that your target sequences
are in *target.fa* and the GISAID database is in *gisaid.fa*, the following two commands will run the whole pipeline:

```bash
$ flaco_blast.sh makedb gisaid.fa [pattern]
$ flaco_blast.sh run target.fa gisaid.clean.fa
```

Note that the *makedb* command only needs to be executed when the GISAID database changes. The *run* command can be called as many times as needed
with the same BLAST database.

The following sections describe each command in detail.

### makedb

This command creates a BLAST database from a FASTA file containing all GISAID sequences. Usage:

```
flaco_blast.sh makedb sequences.fa [remove...]
```

where *sequences.fa* is the file containing GISAID sequences. The command will:

- Clean each sequence by removing gaps ("-") in sequences, and replacing spaces in sequence headers with underscores (e.g., "Hong Kong" -> "Hong_Kong");
- Extract all sequences whose name contains one of the strings in `remove` and save them to file `sequences.excl.fa`;
- Save all other sequences to `sequences.clean.fa`;
- Split the sequences in `sequences.clean.fa` by month and create a BLAST database for each month.

The split sequences and the corresponding BLAST databases will be written to the `split-by-month` directory.

### run

Usage:

```
flaco_blast.sh run target.fa sequences.clean.fa
```

This command runs the split, blast, parse, and extract commands in succession (these commands are documented below, in case you would like
to run them individually).

It takes as input the FASTA file with target sequences, and the name of "clean" sequences file created by the *makedb* command. The output consists 
of the two files written by the *parse* command (main output table and matches table) and the FASTA file created by the *extract* command.

### split

Usage:

```
flaco_blast.sh split target.fa [moretargets.fa...]
```

This command extracts all sequences from the FASTA files to be analyzed and saves them to subdirectories of the form `split-by-month/month`, 
where *month* is extracted from the third field of the sequence header. Each sequence is saved to a separate file called *seqid*.fa, 
where *seqid* is the second field of the sequence header. 

Currently, the `split-by-month` directory name is hardcoded and cannot be changed; 
please rename this directory before launching another execution of *split* to avoid overwriting results from a previous run.

### blast

Usage:

```
flaco_blast.sh blast
```

This command runs BLAST on the sequences under `split-by-month` against the database created in the *makedb* step. It will submit a separate slurm job 
for each directory under `split-by-month`, and will wait until all jobs have terminated.

### parse

Usage:

```
flaco_blast.sh parse outfile matchesfile
```

This command parses the BLAST output files and creates two reports. *outfile* is the main output file. For each input sequence, it contains the best match found in the month before and in the month after the sequence's date. It has the following columns:

* SeqId - ID of the input sequence;
* Date - Input sequence date;
* Before - ID of best match in previous month;
* BefIdent - Identity score of best match in previous month;
* BefLength - Length of best match in previous month;
* After - ID of best match in following month;
* AftIdent - Identity score of best match in following month;
* AftLength - Length of best match in following month.

The file *matchesfile* contains information about all match sequences appearing in *outfile*, without repetitions. Its columns are:

* Sequence description (full header line from FASTA file);
* Sequence ID;
* Country;
* Sequence date;
* Number of occurrences in *outfile*.

### extract

Usage:

```
flaco_blast.sh extract sequences.fa matchesfile outfile.fa
```

This command extracts from the main GISAID sequence database (the one provided to the *makedb* command) all sequences listed in *matchesfile* (produced
by the *parse* command) writing them to *outfile.fa* in FASTA format.

