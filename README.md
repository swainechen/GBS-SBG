# GBS-SBG

## Group B Streptococcus Serotyping by Genome Sequencing

## Quick Start (assuming you have BLAST+ installed and available on your path):
```
git clone https://github.com/swainechen/GBS-SBG
cd GBS-SBG
./GBS-SBG.pl <assembly fasta file> -name <some name>
```

## Introduction
This repository contains a curated reference file which can be used for serotyping Streptococcus agalactiae (Group B Streptococcus, or GBS), in silico with whole genome sequencing data. The reference file (GBS-SBG.fasta) is designed to be usable for both short-read mapping and assembly-based strategies.

The fasta file is designed to be immediately usable with [SRST2](https://github.com/katholt/srst2), but this only means the fasta headers have been formatted specially. Otherwise it should still usable for other reference-based typing pipelines.

To do an assembly-based call, the GBS-SBG.pl script is included. This script uses only core Perl modules.
- It does require the [BLAST+](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download) suite of programs to be installed. The location can be specified with the -blast command line parameter.
- It also requires the GBS-SBG.fasta file. The location can be specified with the -ref command line parameter. However, if this file can't be found, the script will attempt to download a copy directly from this repository.

## Short Read serotype calling
We recommend using SRST2. This fasta file is already formatted for SRST2, so you can run (for example):
```
srst2 --input_pe strainA_1.fastq.gz strainA_2.fastq.gz --output strainA_test --log --gene_db GBS-SBG.fasta
```
For details on SRST2 installation, options, and other documentation, we refer you to the [SRST2 website](https://github.com/katholt/srst2).

## Assembly-based serotype calling
We have implemented a reasonably self-contained script that takes an assembly and produces a serotype call. You can use whatever assembly program you like. We do assume that you have an assembly for a single isolate (i.e. no mixtures or metagenomic assemblies).

To run this, clone this repository (or just download the GBS-SBG.pl script):
```
git clone https://github.com/swainechen/GBS-SBG
```
You should take a look at the built-in help for syntax:
```
./GBS-SBG.pl -help
```

You should see this output:
```
Usage: ./GBS-SBG.pl <assembly_fasta_file> [ -name <string> ] [ -best ] [ -blastn <path_to_blastn> ] [ -ref <GBS-SBG references> ] [ -debug ]

<assembly_fasta_file> should be a regular multi-fasta file with assembled contigs or a complete genome.
You should specify the -name parameter, all output will be prefixed by that string. Defaults to the input filename.
The -best option only prints out one call (with possible uncertainty information). Default behavior is to also print the next best call if any of the following are true for the "best" serotype call:
 - Overall coverage < 0.9
 - Overall percent ID < 90
 - Number of contigs > 1
 - Number of BLASTN hits > 1
 - BitScore for next best serotype is >0.9*BitScore for the best serotype

Requires BLAST+ version 2.7.0 or above. Will look on your path, or you can specify the blastn binary with -blast.
Requires GBS-SBG references (typically GBS-SBG.fasta), will try looking in a few places, otherwise will try to pull directly from https://github.com/swainechen/GBS-SBG.

More complete documentation available at https://github.com/swainechen/GBS-SBG.
```

The basic command is:
```
./GBS-SBG.pl assemblyfile.fasta -name assemblyname
```
If the program has problems finding blastn, you can tell it where to look with the `-blast /path/to/blastn` parameter.
The program will look in the same directory as the script itself for the GBS-SBG.fasta reference file. If it can't find it, it will download that file automatically from this repository into a temporary directory (which will get cleaned up automatically).
The program will also take care of making the blast databases, running blast, and then parsing the output to produce a call.
The default cutoffs correspond to those used in SRST2 - at least 90% identity over at least 90% of the total length of the reference sequence is required, and the highest identity and coverage determine the "best" call.

All output goes to STDOUT, so be sure to pipe it into a file if you want (especially if running on many assemblies).

If you are running this on many assemblies, you should make sure the GBS-SBG.fasta reference is already downloaded and findable by the program. You can specify the location with the `-ref /path/to/GBS-SBG.fasta` parameter. The blast databases will then only be made once and reused for each time the program runs.
You should also consider using the `-best` parameter, which will only output the "best" serotype call, along with any potential uncertainty flags.

## Uncertainty flags (for assembly-based calling)
While the basic call only requires 90% identity over 90% of the reference sequence length, GBS-SBG.pl will provide some uncertainty information in the following cases:
- Overall coverage < 99% of the reference sequence length
- Overall identity < 99% averaged over the aligned length to the reference sequence
- More than 1 blast hit (i.e. HSP) required to account for the total coverage of the reference sequence
- More than 1 contig (i.e. from the input assembly) required to account for the total coverage of the reference sequence
- The total BLASTN Bit Score for the next best serotype call is >95% of the total BLASTN Bit Score for the best hit

When these fall further below certain threshholds, default behavior is to also report the next best hit. This can be suppressed with the `-best` option (i.e. only print the best call). Alternatively, more information (including all the BLASTN data) can be obtained by using the `-debug` option - specifying this also will result in persistence of the temporary directory and all the intermediate files.

## Authors
- Suma Tiruvayipati
- Swaine Chen

## Reference
