=======================================
System requirements
=======================================

FaST requires at least 32 Gb of RAM and a multicore processor with at least 8 cores/16 threads.
More threads may marginally speed up the analysis, especially if using a SSD drive.
All the currently available public datasets compatible with FaST can be analysed with 32 Gb of RAM.

=======================================
Software requirements
=======================================

FaST has been developed to run in a bash shell on Linux OS.
FaST currently requires:
* bedtools v2.30.0 or later
* samtools v1.19.2 (but any version will likely work as the only samtool invoked by FaST is "samtools view")
* STAR v2.7.11a or later 
* spateo-release v1.0.2 (for RNA segmentation only)
* scanpy v1.10.1 or later (optional, for automated exploratory analysis)

It is advisable to take advantage of the requirements.yaml file in the data folder to create a dedicated
environment. Upon successful creation (and activation) of the environment, spateo-release should be installed with pip

.. code-block:: bash
    
    $ pip install spateo-release==1.0.2

=======================================
Installing FaST
=======================================

Fast comes as a set of bash, Perl and Python scripts. To install it you should simply clone the github repository
to a suitable location of your filesystem. Optionally, you can add the /path/to/FaST/scripts/ directory 
to your PATH evironmental variable. This will allow you to run FaST by just typing FaST

.. code-block:: bash
  
    $ FaST -h

instead of 

.. code-block:: bash
  
    $ /path/to/FaST/scripts/FaST -h




