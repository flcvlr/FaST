===============================
Step by step guide to FaST pipeline
===============================


In this tutorial, I am assuming that you have downloaded and dumped the following data from SRA:

* `SRR27331459 <https://trace.ncbi.nlm.nih.gov/Traces/sra?run=SRR27331459>`_
* `SRR27331460 <https://trace.ncbi.nlm.nih.gov/Traces/sra?run=SRR27331460>`_
* `fc_1.fastq.gz.1 <https://sra-pub-src-2.s3.amazonaws.com/SRR27331427/fc_1.fastq.gz.1>`_

after dumping the SRA files (see `here <https://github.com/ncbi/sra-tools/wiki/HowTo:-fasterq-dump>`_ for details) you are supposed to have:

fc_1.fastq.gz.1			This is the file containing the full sequencing of the Illumina flowcell
SRR27331459_1.fastq.gz		This is the file containing read 1 (barcode sequences) of replicate 1
SRR27331459_2.fastq.gz		This is the file containing read 2 (mRNA sequences) of replicate 1
SRR27331460_1.fastq.gz		This is the file containing read 1 (barcode sequences) of replicate 2
SRR27331460_2.fastq.gz		This is the file containing read 2 (mRNA sequences) of replicate 2

You should also download FaST from the github repository (you have several options, use git or download the zipped folder) 
and that the content of the repository have been placed in a directory in the current directory, called FaST.
Inside the FaST/ folder you should have the following folders:

scripts/
docs/
data/

===================
Step 1: FaST-map
===================

You need:

* barcode files in either fastq.gz ("Illumina" protocols: OpenST, Seqscope, Novascope, Nova-ST) or .h5 ("Stereo-seq" protocol, BGI)
* a Hard Drive with sufficient free space, about 200 Gb for "Illumina" protocols, 50 to 300 Gb for BGI (depending on the capture area used)

.. code-block:: bash

   $ FaST/scripts/FaST-map fc_1.fastq.gz.1 tile_map

This command will take a few hours to read the whole 250 Gb fastq file and generate the single tile files in the tile_map/ folder.
After processing each tile, it will issue a line informing you about the number of duplicated barcodes removed from that tile. There are 
about 3500 tiles in an Illumina flowcell.


===================
Step 2: FaST-reference
===================

You need:

* a Hard Drive with sufficient free space, about 30 Gb for the mouse genome.
* a working internet connection to download data from the Gencode repository.

To generate the mouse reference using the latest release available on Gencode you should simply run:

.. code-block:: bash

   $ FaST/scripts/FaST-reference -s mouse

This command will take about 20 minutes to download genomic data and build a STAR index.
Beware, even if you already have a STAR index built on the mouse genome, this step is still required, 
as FaST-reference will also prepare a set of annotation files that should match the precise gtf and fasta you are using here.

===================
Step 3: FaST
===================

You need:

* a PC equipped with a multicore processor (recommended, at least 24 threads)
* 32 Gb of RAM
* A reasonable amount of swap memory (usually on Linux, this is set equal to the RAM).
* about 50 GB of free space on the HD

Now you can run FaST by simply typing the following command:

.. code-block:: bash

   $ FaST/scripts/FaST -1 SRR27331459_1.fastq.gz,SRR27331460_1.fastq.gz -2 SRR27331460_1.fastq.gz,SRR27331460_2.fastq.gz fc_1.fastq.gz.1 -n mouse_head -s mouse -t tile_map -P -S 
   
   
This command will run the FaST pipeline, including:
* collection of the ~ 1 billion barcodes (ETA: ~ 2 minutes)
* collapsing identical barcodes and removing barcodes that are duplicated in the selected tiles (ETA: ~ 1 minute)
* identification of the tiles of the Illumina flowcell used for the libraries (ETA: ~ 1 minute)
* mapping of the reads with STAR (ETA: ~ 30 minutes)
* Building Digital expression matrices (ETA: ~ 5 minutes)
* plotting read density (ETA: ~ 1 minute)
* Segmenting spatial data into cells (ETA: ~ 10 minutes)
* Running a rudimental automated analysis with scanpy to yield a basic clustering, and a UMAP (ETA: ~ 5 minutes)

The whole analysis should take no more than ~ 1 hour on a recent multicore processor equipped with 16 cores/24 threads.

You will find output in the following directories:

* mouse_head/seg_k_3_binsize_20/:	UMAP, Segmented cells in spatial coordinates (with clusters highlighted), a short log reporting the number of cells, histograms reporting cell area, counts and genes
* mouse_head/Aligned.bam:      	Bam file containing the alignments of R2 reads, with barcodes and coordinates as bam tags.
* mouse_head/logs/run.log		A short log reporting the command line options, for your future reference. This log is created when you first run FaST and align reads. If you run again FaST with different options, FaST will append info to this file, so that you will be able to know what you have done. This log is erased and overwritten if you use option -f to run again the pre-processing and alignment of the reads.
* mouse_head/logs/		Several logs about single tiles statistics, including counts and UMIs
     
     
In case you want to re-run FaST with different parameters (usually, to repeat segmentation) just run:

.. code-block:: bash

   $ FaST/scripts/FaST -n mouse_head -s mouse -P -S -K 5 -B 15
   
This will look for aligned reads in the mouse_head/ folder and, if those are found, re-run segmentation, saving output to mouse_head/seg_k_5_binsize_15/. 
FaST will also append a new line to mouse_head/logs/run.log to let you keep track of the different runs you did on the same sample.


   
   






