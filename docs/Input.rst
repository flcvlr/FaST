=============================
Input data
=============================

FaST has been specifically designed to analyse Spatial Transcriptomic datasets 
obtained by repurposing the NGS devices (Illumina or BGI) as an RNA capture device.
Those protocols include:

* SeqScope
* OpenST
* Nova-ST
* NovaScope
* Stereo-seq

Briefly, all these protocols rely on two sequencing steps:

#. HDMI barcode oligos are spot (by soehow complex protocols!) on the Illumina/BGI device 
   and sequenced (1st sequencing). The output of the first sequencing is tipycally a huge 
   (hundreds of Gb) set of fastq.gz file(s), containing all the sequences of the HDMI barcodes 
   in the Illumina flowcell or BGI device. Occasionally, these data may be provided in small chunks 
   corresponding to the tiles of a (set of) experiment(s). BGI provides these data in HDF5 format.
#. After recovering the Illumina flowcell (with the HDMI barcodes attached) and using the polyT end of the 
   HDMI barcodes to capture RNA from a tissue slice which is lying on a small area of the flowcell, 
   a new library is prepared and sequenced (2nd sequencing).
  
**FaST-map** will process the (usually huge) fastq.gz (or .h5, for Stereo-seq/BGI) files obtained in the 1st sequencing.

**FaST-reference** should be fed with a (non-gzipped) genome fasta file and a corresponding (non-gzipped) gtf file.

**FaST** will instead process the paired fastq.gz files obtained in the 2nd sequencing.

