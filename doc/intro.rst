.. _intro:

###############
Getting Started
###############

What is SHEAR?
==============

SHEAR (Simply Handler for Error and Adapter Removal) is a short-read trimmer for high-throughput sequencing fastq files.
SHEAR first scans the fastq file(s) and automatically detects likely adapter and primer contaminants.  Then it calls 
Scythe (https://github.com/vsbuffalo/scythe), which removes reads using a Bayesian error-tolerant approach that 
effectively removes adapters even if the adapter itself contains a sequence variation due to sequencing error.
SHEAR then trims and filters reads based on minimum quality and content cutoffs.  Additionally, SHEAR is designed to 
automate the simultaneous trimming and concatenation of multiple paired read files into a single trimmed file ready
for mapping or assembly.  

Requirements
============

* Python 2.7.x or 3.x (3.x recommended)

Optional
--------

* Scythe: https://github.com/vsbuffalo/scythe (**Strongly recommended**)

Installation
============

No installation is necessary, simply clone the repository from GitHub. Scythe should be installed according to its instructions.::

  git clone https://www.github.com/jbpease/shear


Preparing your data
===================
Standard fastq files with four lines per entry (header, sequence, gap, quality) should be used.  ABI solid colorspace reads are not currently supported.  When using paired-end mode reads only need to be sorted and contain the same number of reads if the ``-U/--filter-unpaired`` mode is used, since this will remove both reads from a pair when either of them is filtered out.


Usage
=====

Paired-end
----------

 ::

  python shear.py --fq1 FASTQ.SAMPLE1.p1.fastq FASTQ.SAMPLE2.p1.fastq ... --fq2 FASTQ.SAMPLE1.p2.fastq FASTQ.SAMPLE2.p2.fastq ...  --out1 FASTQ.sheared.p1.fq --out2 FASTQ.sheared.p2.fq

Single-read
----------

 ::

  python shear.py --fq1 FASTQ.SAMPLE1.p1.fastq FASTQ.SAMPLE2.p1.fastq ... --fq2 FASTQ.SAMPLE1.p2.fastq FASTQ.SAMPLE2.p2.fastq ...  --out1 FASTQ.sheared.p1.fq --out2 FASTQ.sheared.p2.fq


Config file alternative
-----------------------
Alternatively to a full set of command line arguments you can enter a single positional argument that points in a text file.  You can then specify command line arguments more neatly over several lines as in the example:

 ::

  --fq1 
  FASTQ.SAMPLE1.p1.fastq 
  FASTQ.SAMPLE2.p1.fastq 
  --fq2 
  FASTQ.SAMPLE1.p2.fastq 
  FASTQ.SAMPLE2.p2.fastq 
  --out1 FASTQ.sheared.p1.fq 
  --out2 FASTQ.sheared.p2.fq 

Version History
===============

v. 2017-06-20
-------------
Major upgrade, fixed issues with compatibility with Scythe, 
added automatic adapter finding with adapt.py, 
included new manual and documentation with Sphinx.  
WARNING: program options have changed significantly to
be more conventional and consistent in format. Please be 
advised that old commands will need to be updated in terms
of flag names.

v. 2015-09-13 
--------------
Fixes for Python3 compatbility , add gzip capability

v. 0.007
--------
Fixes to default parameters and option processing

v. 0.006
---------
Minor fixes to default parameters

v. 0.005 
--------
Major update, fixed filtered output, clean-up, removed GC content filter

v. 0.004
--------
Added support for combining multiple pairs of input files

v. 0.003
--------
Alpha Release
