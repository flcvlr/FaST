===============================================
FaST-map
===============================================

The FaST-map script will read a single gzipped fastq file and populate an output folder with the tile-specific
HDMI barcode collection files and the corresponding indexes.
Upon completion of the run the output folder will contain a few thousand files named after Illumina flowcell tiles
(e.g. 2_2456), two files for each tile, (2_2456.txt.gz and 2_2456.idx). While FaST-map is running you can monitor 
its work by simply opening a new terminal window and listing the contents of the output directory. FaST-map is 
expected to process (and therefore generate outpu files) several tiles per minute. Generation of .idx files will
start with a small lag, approximately after the first 20 tiles have been processed.
 
Command syntax:

.. code-block:: bash

   $ FaST-map <path/to/first/seq.fastq.gz> </target/folder>

The first argument should be the path to the gzipped fastq file containin the first sequencing results. 
If the results were (rigorously, i.e. the fastq file was split in between records before compression of each chunk)
split into several different files, those can be simply concatenated by:

.. code-block:: bash
    $ cat file1.fastq.gz file2.fastq.gz file3.fastq.gz > big.file.fastq.gz

or, in case you do not care of saving on your filesystem the whole file (which will just occupy space for duplicated
data)

.. code-block:: bash

   $ FaST-map <(cat file1.fastq.gz file2.fastq.gz file3.fastq.gz) </target/folder>
   
should also work.
    
The second argument is the path to a folder (FaST-map will create it if that does not exist yet) to save the barcodes
per tile map. You will have to provide this path to FaST as an argument to the ``-t`` option in the analysis of the 
2nd sequencing.

================
Notes on fastq.gz input files
================
    
*Please bear in mind that order matters in this case as the reads are expected to occur in order by tile as per 
Illumina default policy. If you rejoin fastq.gz chunks in the wrong order data of the tiles across two files will
be corrupted.
*If you run only the first chunk of a fastq.gz file you have trunkated (because you know that the tiles you are 
interested in are in that chunk and the file was too large for your repository or transfer tool, or just to save disk 
space) FaST-map will throw an error upon reaching the trunkated end of the gzipped file. In that case a few corrupted 
tiles may be generated. Those will correspond to files in the output directory with size = 0. Remove all (both .txt.gz 
and .idx) output for all tiles in which either of the files has size = 0 bytes or later FaST will fail.
*Please note that if a large fastq.gz has been split into multiple chunks arbitrarily (i.e. without decompressing 
first, then splitting in between two fastq records and finally recompressing the resulting chunks) it should be 
appropriately rejoined before processing. If you are unsure about that, ask those who split the file how to 
restore it.


