# SHEAR: Simple Handler for Error and Adapter Removal

This program provides an integrated system to trim and filter Illumina short reads. 
Barcode-customized adapter tags are prepared for removal by Scythe (https://github.com/vsbuffalo/scythe),
then reads are filtered by length, quality scores, and in the case of highly repetition of sequence.
This provides a quick, integrated, and installation-free method to reduce error in your Illumina short reads.
You can also skip Scythe and use your own adapter trimmer, but still use SHEAR for all other trimming.

## Key features:
* Simultaneous pair trimming and discarding of orphaned reads (i.e. the other pair of a filtered read). 
* Poly-A and Poly-T tail trimming (good for RNA-Seq)
* Removes dinucleotide repeat reads.
* Automatically builds your adapter sequences with custom barcodes from your samples.

##AUTHORS: James B. Pease

##Contents:
shear.py  (main script)

##Requirements:
Python 2.5+, or 3.x

**Optional, but strongly recommended: **
Scythe (any version, https://github.com/vsbuffalo/scythe)

##Installation:
+ Install Scythe ([https://github.com/vsbuffalo/scythe](https://github.com/vsbuffalo/scythe), *optional but strongly recommdended*)
+ Download and run script.

##Usage:
Full usage options can be found by running: 
`> shear.py --help`
