========================
FaST
========================

This is the main command to run FaST, after generating a reference (with FaST-reference) and a map of 
the Illumina flowcell used to capture the RNA (using FaST-map).


Usage
------------------------

Command syntax:

.. code-block:: bash

   $ FaST -1 <R1_fq.gz> -2 <R2.fq.gz> -s <species> -n <name of sample> -t <folder/with/tiles/map>   [-j <integer> -T <float> -K <integer> -f -S -R -P -I -M <float>]




===========   ===================
Option         Description
===========   ===================
**-1**        path to a gzipped fastq file containing the sequence of the barcodes. If the data is split in several pairs 
	      of R1 and R2 files, list all R1 files separated by commas (i.e. R1_seq1.fastq.gz,R1_seq2.fastq.gz,...). 
	      Do not leave empty spaces in between files paths and commas.
**-2**	      path to a gzipped fastq file containing the RNA sequencing output.If the data is split in several pairs 
	      of R1 and R2 files, list all R2 files separated by commas. The order MUST match that of the R1 files. 
	      Do not leave empty spaces in between files paths and commas.
**-f**        No argument. Forces full re-analysis of the sample, including preprocessing alignment and DGE. If aligned
	      reads are found in the sample folder, default behaviour is to run only segmentation analysis (if -S is set).
**-j**	      max number of jobs to run in parallel during alignment (default=20, suggested value: num_threads -4).
**-n**	      name for the sample you want to analyse. The software will create a folder with this name and all output
	      will be written in that folder.
**-p**	      ST protocol. Currently FaST supposrts either the "Illumina" (Seq-scope, Nova-scope, OpenST, Nova-ST) or "Stereo-seq".
**-r**	      ratio of barcodes found to select a tile. Default is 0.1 (i.e. only retain tiles in which at least 10%
	      of the barcodes assigned to that tile are found in the sequencing)
**-s**	      either "human" or "mouse", instructs FaST to look for the appropriate annotation which is stored 
              in the FaST/reference/ folder, and shoud have been generated with FaST-reference ( FaST-reference -h for help)
**-t**	      path to a folder containing barcodes with coordinates. This is generated from the first sequencing round
	      of seq-scope, OpenST, Stereo-seq or Nova-ST protocols using FaST-map.
**-B**	      Binsize to segment capture area during cell segmentation using Spateo (default: 20. Before changing this 
	      parameter, take the time to carefully read Spateo-release documentation.)
**-I**        If set, reads mapping to intergenic regions will be taken as evidence of a cell on that "puck" and will 
	      therefore contribute to cell segmentation. (Experimental, setting this option is not currently reccomended)
**-K**        kernel density to use in cell segmentation. 11 > odd_integer > 1 (default=3).
**-M**	      If set, reads mapping to multiple genes will be taken as evidence of a cell on that "puck" and will 
	      therefore contribute to cell segmentation. (Experimental, setting this option is not currently reccomended)
**-P**	      If set, gene classified as "pseudogenes" will be removed from the annotation before running FaST. About 10% 
	      of the reads are discarded as "mapping to multiple genes" just becase a (possibly very lowly or totally not
	      expressed) matching pseudogene exists in the genome. By setting this option, one can recover a significant
	      fraction of the reads. However, since many people may be legitimately interested in quantifying also pseudogenes,
	      default behaviour is to retain pseudogenes in the annotation. If you are not interested in tracking pseudogenes,
	      setting this option is reccommended.
**-S**	      No argument. If this is passed to FaST it will also perform RNA segmentation with spateo-release package.


===========   ===================


	

Notes
------------------------

* FaST will perform complete analysis of High Resolution Spatial Transcriptomic data assuming that it follows the OpenST/Seqscope/Nova-ST or Stereo-seq format, and in particular: 

1) the reads in R1_fq.gz file(s) should match the barcode files in the folder provided as an argument to -t. To obtain details on how to generate the barcode files, please run 

.. code-block:: bash
   $ path/to/FaST/scripts/FaST-map -h 

2) the reads in R2_fq.gz will contain a 9nt UMI at the 5'end followed by the "sense" mRNA sequence

FaST will take care of read trimming, there is no need to trim R1 and R2 before running FaST.

The output will be saved in the sample_name/ folder and will include:
-   bam files for alignment to reference genome of your species
-   several metrics (# of barcodes, # of segmented cells, automated analysis output including segmented cells area statistics, UMAP projection of clustered cells, spatial image of cluestered cells on tissue slice)
-   a full set of data on the RNA segmented cells including:
	- .h5ad anndata format files to load your data into Scanpy/Seurat
	- images of the RNA segmentation output

Requirements:
FaST requires the bash shell, with support for command substitution, sed, awk, sort, uniq, cut, grep and xargs.
These features are already present in any recent linux distribution which has a bash shell.

The only required softwares are:

* pigz

* perl-GD

* bedtools v2.30.0 or later

* samtools v1.19.2 (but any version will likely work as the only samtool invoked by FaST is "samtools view")

* STAR v2.7.11a or later 

* spateo-release v1.0.2 (for RNA segmentation only)

* scanpy v1.10.1 or later (optional, for automated exploratory analysis)
      
All dependencies may be installed with conda, please refer to the `FaST_env.yml <https://github.com/flcvlr/FaST/blob/main/data/FaST_env.yml>`_  file in the "data" folder to create a reproducible environment for FaST.

Hardware and time considerations:
FaST will take advantage of available threads, with the limit enforced by option -j in place for bowtie2 and STAR alignment. 
FaST will fit in 32Gb RAM. A typical sample (15 square mm slice of tissue, ~1 Billion reads) should take less than 2 hours 
with 24 threads for complete analysis. Moderately larger tissues slices or deeper sequencing should result in a linear increase of time
 requirements. If your data is split on several pairs of R1.fastq.gz/R2.fastq.gz files you may provide all of them at once, 
 check details of -1 and -2 options in the options descriptions.





