========================
FaST
========================

This is the main command to run FaST, after generating a reference (with FaST-reference) and a map of 
the Illumina flowcell used to capture the RNA (using FaST-map).

======================
Usage
======================

Command syntax:

.. code-block:: bash

   $ FaST -1 <R1_fq.gz> -2 <R2.fq.gz> -s <species> -n <name of sample> -t <folder/with/tiles/map>  -c <file/with/tiles/offsets>  [-j n -T <float> -k <integer> -m <string> -f -S -R <float>]




===========   ===================
Option         Description
===========   ===================
**-1**        path to a gzipped fastq file containing the sequence of the barcodes. If the data is split in several pairs 
	      of R1 and R2 files, list all R1 files separated by commas.
**-2**	      path to a gzipped fastq file containing the RNA sequencing output.If the data is split in several pairs 
	      of R1 and R2 files, list all R2 files separated by commas. The order MUST match that of the R1 files.
**-s**	      either "human" or "mouse", instructs FaST to look for the appropriate annotation which is stored 
              in the FaST/reference/ folder, and shoud have been generated with FaST-reference ( FaST-reference -h for help)
**-n**	      name for the sample you want to analyse. The software will create a folder 
	      with this name (and path) all output will be written in that folder.
**-t**	      path to a folder containing barcodes with coordinates. This is generated from the 
	      first sequencing round of seq-scope, openST or novaseq protocols using FaST-map.
**-R**	      ratio of barcodes found to select a tile. Default is 0.1 (i.e. only retain tiles in which at least 10%
	      of the barcodes assigned to that tile are found in the sequencing)
**-j**	      max number of jobs to run in parallel during alignment (default=20, suggested value: num_threads -4).
**-k**        kernel density to use in cell segmentation. 11 > odd_integer > 1 (default=7 for human, k=5 for mouse).
**-T**        Threshold to use for EM_BP function during RNA segmentation. Values between 0.3 and 0.6 are generally
              appropriate. If omitted, spateo-release will evaluate one for each tile (but this can result into uneven 
              segmentation efficiency across tiles).
**-m**        Method to use to pick nuclear labels for segmentation. Either "apex" or "nuc". "apex" instructs
              FaST to run nuclear segmentation using Apex-seq set of nuclear genes. "nuc" innstead will first run
              a differential expression to identify genes co-expressed with introns and miR host genes and then use 
              those genes to guide nuclear segmentation. Default:"apex".
**-f**        No argument. Forces full re-analysis of the sample, including preprocessing alignment and DGE. If aligned
	      reads are found in the sample folder, default behaviour is to run only segmentation analysis (if -S is set).
**-H**	      No argument. Use if running FaST on an HPC cluster. This option instructs FaST to process each couple of 
	      fastq.gz files in parallel. Has no effect in case all reads are in a single pair (i.e. R1 and R2) of 
	      fastq files. Requires about 20 Gb of RAM for each pair of R1 and R2 fastq files.
**-S**	      No argument. If this is passed to FaST it will also perform RNA segmentation with spateo-release package.

===========   ===================


	
======================
Notes
======================

* FaST will perform complete analysis of High Resolution Spatial Transcriptomic data assuming that it follows the OpenST or seqscope format, and in particular: 

1) the reads in R1_fq.gz file(s) should match the barcode files in the folder provided as an argument to -t. To obtain details on how to generate the barcode files, please run FaST-map -h 

2) the reads in R2_fq.gz will contain a 9nt UMI at the 5'end followed by the "sense" mRNA sequence

FaST will take care of read trimming, there is no need to trim R1 and R2 before running FaST.

The output will be saved in the sample_name/ folder and will include:
-   bam files for alignment to reference genome of your species
-   several metrics (# of barcodes, # of segmented cells, automated analysis output including segmented cells area statistics, UMAP projection of clustered cells, spatial image of cluestered cells on tissue slice)
-   a full set of data on the RNA segmented cells including:
	- .h5ad anndata format files to load your data into Scanpy
	- images of the RNA segmentation output

Requirements:
FaST requires the bash shell, with support for command substitution, sed, awk, sort, uniq, cut, grep and xargs.
These features are already present in any recent linux distribution which has a bash shell.

The only required softwares are:
* bedtools v2.30.0 or later
* samtools v1.19.2 (but any version will likely work as the only samtool invoked by FaST is "samtools view")
* STAR v2.7.11a or later 
* spateo-release v1.0.2 (for RNA segmentation only)
* scanpy v1.10.1 or later (optional, for automated exploratory analysis)

To ensure appropriate segmentation by spateo-release v 1.0.2 the following versions of spateo-release dependencies should be installed
* scipy==1.12.0
* matplotlib==3.7.1
* louvain==0.8.0
* pyvista==0.42
      
All dependencies may be installed with conda, please refer to the FaST_env.yml file in the "data" folder to create a reproducible environment for FaST.

Hardware and time considerations:
FaST will take advantage of available threads, with the limit enforced by option -j in place for bowtie2 and STAR alignment. 
FaST will fit in 32Gb RAM, A typical 10 square mm tissue slice experiment will have ~ 40M "real" barcodes + ~ 60M "wrong barcodes 
resulting from sequencing errors, requiring ~20 Gb RAM during barcode selection. Low sequencing quality of R1 file may 
increase the number of wrong barcodes and thus result into an increased RAM requirement.
A typical sample (15 square mm slice of tissue, ~700M R2 reads) should take ~ 2 to 3 hours with 24 threads for complete analysis.
Moderately larger tissues slices or deeper sequencing should result in a linear increase of time requirements. 
If your data is split on several pairs of R1.fastq.gz/R2.fastq.gz files you may consider using option -H if running on a HPC or
other system with at least 64Gb RAM. This option will require significantly more RAM (about 20 Gb for each pair of R1/R2 fastq files).





