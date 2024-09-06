========================
FaST-reference
========================

This command should be run before processing a library, in order to build the appropriate reference.
FaST-reference requires the reference genome fasta file and a annotation file in gtf format.
For human and mouse, if the fasta and gtf files are not provided, FaST-refernce will instead download 
these information from Gencode website, defaulting to the latest release. However, you can specify a
given release of the Gencode data.
Currently the only supported species are human and mouse, but support for custom genomes will be
added soon


======================
Usage
======================

Command syntax:

.. code-block:: bash

   $ FaST-reference -s <species>  [ -v gencode_version | -f genome.fa -g annotation.gtf ] 


Options:
-s	Mandatory. Currently either 'human' or 'mouse'.
-v	Optional, gencode annotation version (i.e. 'v46' for human version 46, 'vM25' for mouse version 25).
	If omitted and no fasta or gtf are selected, the latest Gencode release for the selected species will be 
	downloaded.
-f	Optional. Path to a fasta file corresponding to the species of interest. 
	If omitted, the latest Gencode data for the selected species will be downloaded.
-g	Optional, Path to gtf format file containing transcriptome annotation for the selected species.
	If omitted, the latest Gencode data for the selected species will be downloaded.
	
======================
Notes
======================
	
* If you do not have a local fasta and gtf file, internet connection will be required to download gencode data.
* If specifying a local fasta and/or gtf file please make sure that they correspond to the same genome assembly.
* If specifying only a fasta (or gtf) file, please make sure that they match the same genome assembly of
the Gencode gtf (or fasta) that FaST-reference will download (the latest available if -v is omitted)




