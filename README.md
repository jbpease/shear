# SHEAR: Simple Handler for Error and Adapter Removal

## What is SHEAR?

This program removes adapters, trims, and quality-filters Illumina short reads. 
The program first will automatically find adapters sequences in your reads and remove them using Scythe (https://github.com/vsbuffalo/scythe), then reads are filtered by length, quality scores, and in the case of highly repetition of sequence. This provides a quick, integrated, and installation-free method to reduce error in your Illumina short reads. You can also skip Scythe and use your own adapter trimmer, but still use SHEAR for all other trimming.

## Key features

* Simultaneous pair trimming and discarding of orphaned reads (i.e. the other pair of a filtered read). 
* Poly-A and Poly-T tail trimming (good for RNA-Seq)
* Removes dinucleotide repeat reads.
* Automatically finds adapter sequences from your samples.

## Manual
See the SHEAR manual here: (https://github.com/jbpease/shear/blob/master/shear.pdf)

## Authors

James B. Pease
https://www.peaselab.github.io

## Requirements

Python 3.x (recommend) or 2.7.x

**Optional, but strongly recommended:**

Scythe (any version, https://github.com/vsbuffalo/scythe)

## License
This file is part of shear.

shear is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as publihed by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

shear is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. You should have received a copy of the GNU General Public License along with Foobar. If not, see (http://www.gnu.org/licenses/).
