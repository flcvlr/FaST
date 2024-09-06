=============================
Input data
=============================

FaST has been specifically designed to analyse Spatial Transcriptomic datasets 
obtained by repurposing the Illumina Flowcells as a RNA capture device.
Those protocols include:
* SeqScope
* Open-ST
* Nova-ST
* NovaScope

Briefly, all these protocols rely on two sequencing steps:
#. HDMI barcode oligos are spot on the Illumina flowcell and sequenced (1st sequencing)
   the output of the first sequencing is tipycally a huge (hundreds of Gb) set of fastq.gz files,
   containing all the sequences of the HDMI barcodes in the Illumina flowcell
   Occasionally, these data may be provided in small chunks corresponding to the tiles of a (set of) experiment(s).
#. After recovering the Illumina flowcell (with the HDMI barcodes attached) and using the polyT end of the 
   HDMI barcodes to capture RNA from a tissue slice which is lying on a small area of the flowcell, 
   a new library is prepared and sequenced (2nd sequencing).
  
*FaST-map* will process the (usually huge) fastq.gz files obtained in the 1st sequencing.
*FaST-reference* should be fed with a (non-gzipped) genome fasta file and a corresponding (non-gzipped) gtf file.
*FaST* will instead process the paired fastq.gz files obtained in the 2nd sequencing.

