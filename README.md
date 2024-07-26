# FaST
Fast analysis of Spatial Transcriptomics

FaST is designed to allow quick analysis of large spatial transcriptomics datasets based on the 
repurposing of the Illumina flow-cell as an RNA capture device. Published protocols supported by FaST
include OpenST, seqscope and novascope.

Requirements:

bash shell environment
bedtools v2.30.0 or later
samtools v1.19.2 (but any version will likely work as the only samtool invoked by FaST is "samtools view")
STAR v2.7.11a or later 
bowtie2 v2.5.4 or later
spateo-release v1.0.2 (for RNA segmentation only)
scanpy v1.9.8 or later (optional, for automated exploratory analysis)


