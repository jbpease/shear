.. adapt:

adapt
=====

Description
-----------

This progam (as part of SHEAR) searches adapter sequences
and generates an adapter file for use with SHEAR/Scythe.


Parameters
----------

``-h/--help``
^^^^^^^^^^^^^

**Description:** show this help message and exit

**Type:** boolean flag



``--fq1`` (required)
^^^^^^^^^^^^^^^^^^^^

**Description:** one or more fastq file paths, separated by spaces

**Type:** file path; **Default:** None



``-o/--out`` (required)
^^^^^^^^^^^^^^^^^^^^^^^

**Description:** output FASTA of adapters detected

**Type:** file path; **Default:** None



``--fq2``
^^^^^^^^^

**Description:** one or more fastq file paths separated by spaces, only use this for paired-end fastq files and enter these files in the same order as their counterparts in --fq1

**Type:** file path; **Default:** None



``-k/--end-klength``
^^^^^^^^^^^^^^^^^^^^

**Description:** Length of end kmer to tabulate for possible adapter matches.

**Type:** integer; **Default:** 16



``-m/--mode``
^^^^^^^^^^^^^

**Description:** known=only use list of known adapters;endmer=search for common 3'end sequences;both=both known and endmers

**Type:** None; **Default:** known

**Choices:** ('known', 'endmer', 'both')


``--quiet``
^^^^^^^^^^^

**Description:** Suppress progress messages

**Type:** boolean flag



``-M/--min-match``
^^^^^^^^^^^^^^^^^^

**Description:** Minimum proportion of read match required to report the kmer as a possible match.

**Type:** float; **Default:** 0.0001



``-N/--number-of-reads``
^^^^^^^^^^^^^^^^^^^^^^^^

**Description:** Number of reads to search in each fastq

**Type:** integer; **Default:** 200000


