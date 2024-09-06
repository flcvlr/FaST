========================
FaST-reference
========================

This command should be run before processing a library, in order to build the appropriate reference.
FaST-reference requires the reference genome fasta file and a annotation file in gtf format.
For human and mouse, if the fasta and gtf files are not provided, FaST-refernce will instead download 
these information from Gencode website, defaulting to the latest release. However, you can specify a
given release of the Gencode data.
If your species is not human nor mouse, you have only the option of providing a fasta and gtf file.
FaST-reference assumes that the fasta is a plain non gzipped fasta and that the gtf file is also a 
plain, non gzipped gtf file.
The reference will be saved to a folder with the species name (overwriting its content if it already exists).
However, by providing different names, you can keep several different references for the same species (i.e. 
``hg38``, ``hg39``).


