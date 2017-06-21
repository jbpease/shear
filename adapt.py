#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
This progam (as part of SHEAR) searches adapter sequences
and generates an adapter file for use with SHEAR/Scythe.
"""

import sys
import os
import argparse
import re
import gzip

_LICENSE = """

SHEAR: Simple Handler for Error and Adapter Removal
James B. Pease
http://www.github.com/jbpease/shear

This file is part of SHEAR.

SHEAR is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

SHEAR is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with SHEAR.  If not, see <http://www.gnu.org/licenses/>.
"""

ADAPTERS = [
    # GATCGGAAGAGCACACGTCTGAACTCCAGTCAC[6 bases]  ATCTCGTATGCCGTCTTCTGCTTG
    # GATCGGAAGAGCACACGTCTGAACTCCAGTCAC[6 bases]CAATCTCGTATGCCGTCTTCTGCTTG
    # GATCGGAAGAGCACACGTCTGAACTCCAGTCAC[6 bases]GTATCTCGTATGCCGTCTTCTGCTTG
    # GATCGGAAGAGCACACGTCTGAACTCCAGTCAC[6 bases]GAATCTCGTATGCCGTCTTCTGCTTG
    # GATCGGAAGAGCACACGTCTGAACTCCAGTCAC[6 bases]CGATCTCGTATGCCGTCTTCTGCTTG
    # GATCGGAAGAGCACACGTCTGAACTCCAGTCAC[6 bases]ACATCTCGTATGCCGTCTTCTGCTTG
    # GATCGGAAGAGCACACGTCTGAACTCCAGTCAC[6 bases]TTATCTCGTATGCCGTCTTCTGCTTG
    # GATCGGAAGAGCACACGTCTGAACTCCAGTCAC[6 bases]TAATCTCGTATGCCGTCTTCTGCTTG
    # GATCGGAAGAGCACACGTCTGAACTCCAGTCAC[6 bases]ATATCTCGTATGCCGTCTTCTGCTTG
    # GATCGGAAGAGCACACGTCTGAACTCCAGTCAC[7 bases]  ATCTCGTATGCCGTCTTCTGCTTG
    # GATCGGAAGAGCACACGTCTGAACTCCAGTCAC[6 bases]  ATCTCGTATGCCGTCTTCTGCTTG
    "A?GATCGGAAGAGCACACGTCTGAACTCCAGTCAC([ATGC]{6,8})"
    "ATCTCGTATGCCGTCTTCTGCTTG",
    "AATGATACGGCGACCACCGAGATCTACAC([ATGC]{6,8})"
    "ACACTCTTTCCCTACACGACGCTCTTCCGATCT",
    # AATGATACGGCGACCACCGAGATCTACAC[5 bases]ACACTCTTTCCCTACACGACGCTCTTCCGATCT
    # AATGATACGGCGACCACCGAGATCTACAC[5 bases]TCGTCGGCAGCGTC
    # AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT
    # AATGATACGGCGACCACCGAGATCTACACGTTCAGAGTTCTACAGTCCGA
    # AATGATACGGCGACCACCGACAGGTTCAGAGTTCTACAGTCCGA
    "AATGATACGGCGACCACCGAGATCTACAC([ACGT]{5,7})"
    "ACACTCTTTCCCTACACGACGCTCTTCCGATCT",
    "AATGATACGGCGACCACCGAGATCTACAC([ACGT]{5,7})TCGTCGGCAGCGTC",
    "AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT",
    "AATGATACGGCGACCACCGAGATCTACACGTTCAGAGTTCTACAGTCCGA",
    "AATGATACGGCGACCACCGACAGGTTCAGAGTTCTACAGTCCGA",
    # CAAGCAGAAGACGGCATACGAGAT[6 bases]GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT
    # CAAGCAGAAGACGGCATACGAGAT[6 bases]GTGACTGGAGTTCCTTGGCACCCGAGAATTCCA
    # CAAGCAGAAGACGGCATACGAGATCGGTCTCGGCATTCCTGCTGAACCGCTCTTCCGATCT
    # CAAGCAGAAGACGGCATACGAGAT[6 bases]GTGACTGGAGTTC
    # CAAGCAGAAGACGGCATACGAGAT[7 bases]GTCTCGTGGGCTCGG
    # CAAGCAGAAGACGGCATACGAGCTCTTCCGATCT
    "CAAGCAGAAGACGGCATACGAGAT([ATGC]{5,7})GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT",
    "CAAGCAGAAGACGGCATACGAGAT([ATGC]{5,7})GTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
    "CAAGCAGAAGACGGCATACGAGATCGGTCTCGGCATTCCTGCTGAACCGCTCTTCCGATCT",
    "CAAGCAGAAGACGGCATACGAGAT([ATGC]{5,7})GTGACTGGAGTTC",
    "CAAGCAGAAGACGGCATACGAGAT([ATGC]{5,7})GTCTCGTGGGCTCGG",
    "CGGTTCTTCCCTGCCGAACCCTATCTTCGTCGGCAGCGTCAGATGTGTATAAGAGACAGTACGCTTGCAT",
    "TTTTTAATGATACGGCGACCACCGAGATCTACACACACTCTTTCCCTACACGACGCTCTTCCGATCT",
    "ATGATACGGCGACCACCGAGATCTACACGTTCAGAGTTCTACAGTCCGACG",
    # ACAGGTTCAGAGTTCTACAGTCCGAC
    # ACAGGTTCAGAGTTCTACAGTCCGACATG
    # CCGACAGGTTCAGAGTTCTACAGTCCGACATG
    # CGACAGGTTCAGAGTTCTACAGTCCGACGATC
    # GTTCAGAGTTCTACAGTCCGACGATC
    "C?C?G?A?C?A?G?GTTCAGAGTTCTACAGTCCGACA?T?G?A?T?C?",
    #
    # GATCGGAAGAGCACACGTCTGAACTCCAGTCAC
    # AGATCGGAAGAGCACACGTCT
    # GATCGGAAGAGCACACGTCTGAACTCCAGTCAC
    # AGATCGGAAGAGCACACGTCT
    # GATCGGAAGAGCACACGTCT
    "A?GATCGGAAGAGCACACGTCTGAACTCCAGTCAC",
    "A?GATCGGAAGAGCACACGTCT",
    "AGATCGGAAGAGCGTCGGTGTAGGGAAAG",
    "AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAG",
    "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA",
    "ATCTCGTATGCCGTCTTCTGCTTG",
    "CAAGCAGAAGACGGCATACGA",
    "CAAGCAGAAGACGGCATACGAGCTCTTCCGATCT",
    "ACACTCTTTCCCTACACGACGCTCTTCCGATCT",
    "CGGTCTCGGCATTCCTGCTGAACCGCTCTTCCGATCT",
    "CTGTCTCTTATACACATCTCCGAGCCCACGAGAC",
    "CTGTCTCTTATACACATCTGACGCTGCCGACGA",
    "GAAUUCCACCACGUUCCCGUGG",
    "GATCGGAAGAGCGGTTCAGCAGGAATGCCGAG",
    "GATCGGAAGAGCTCGTATGCCGTCTTCTGCTTG",
    "GATCGTCGGACTGTAGAACTCTGAAC",
    "GCCTTGGCACCCGAGAATTCCA",
    "GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG",
    "GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT",
    "GUUCAGAGUUCUACAGUCCGACGAUC",
    "TCGGACTGTAGAACTCTGAAC",
    "TCGTATGCCGTCTTCTGCTTG",
    "TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG",
    "TGGAATTCTCGGGTGCCAAGG",
    "AGATCGGAAGAG",
    "CTGTCTCTTATA",
    "CGCCTTGGCCGT"
    "ATCGTCGGACT",
    "GGAATTCTCGG",
    ]


def gopen(path, mode):
    """Automatically invokes gzip file opening when path ends with .gz"""
    return path.endswith('.gz') and gzip.open(path, mode) or open(path, mode)


def revcom(adapter):
    revstr = []
    i = 0
    while i < len(adapter):
        if adapter[i] in 'ATGCU':
            if adapter[i+1:i+2] in '?+*':
                revstr.append(adapter[i:i+2])
                i += 2
            else:
                revstr.append(adapter[i])
                i += 1
        elif adapter[i] == '(':
            j = adapter.find(")", i)
            revstr.append(adapter[i:j+1])
            i = j + 2
    revstr = [complement(x) if ('(' not in x and "?" not in x)
              else x for x in revstr[::-1]]
    revcom = ''.join(revstr)
    return revcom


def complement(seq):
    """Make complement of DNA/RNA sequence"""
    return seq[::-1].lower().replace("a", "U" if "u" in seq else "T").replace(
        "c", "G").replace("g", "C").replace("t", "A").replace(
            "u", "A")


class AdaptFind(object):

    def __init__(self):
        self.known_adapters = []
        for adapter in ADAPTERS:
            self.known_adapters.append(re.compile(adapter))
            self.known_adapters.append(re.compile(revcom(adapter)))
        self.found = {}
        self.barcodes = set([])
        # print([x.pattern for x in self.known_adapters])

    def find_known_adapter(self, read):
        for pattern in self.known_adapters:
            match = re.match(pattern, read)
            if match is not None:
                self.found[match.group(0)] = self.found.get(
                    match.group(0), 0) + 1
                if len(match.groups()) > 0:
                    if match.groups()[0] is not None:
                        if len(match.groups()[0]) > 0:
                            self.barcodes.update([match.groups()[0]])
        return ''


class Adaptamers(object):

    def __init__(self):
        self.endmers = {}

    def add_read_end(self, read, k=16):
        endmer = read[len(read) - k:]
        if "N" in endmer:
            return ''
        self.endmers[endmer] = self.endmers.get(endmer, 0) + 1
        return ''


def generate_argparser():
    parser = argparse.ArgumentParser(
        prog="adapt.py",
        description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        epilog=_LICENSE)
    parser.add_argument("--fq1", required=True, nargs='*',
                        type=os.path.abspath,
                        help=("one or more fastq file paths, "
                              "separated by spaces"))
    parser.add_argument("--fq2", nargs='*', type=os.path.abspath,
                        help=("one or more fastq file paths separated "
                              "by spaces, only use this for "
                              "paired-end fastq files and enter "
                              "these files in the same order "
                              "as their counterparts in --fq1"))
    parser.add_argument("-m", "--mode", choices=("known", "endmer", "both"),
                        default="known",
                        help=("known=only use list of known adapters;"
                              "endmer=search for common 3'end sequences;"
                              "both=both known and endmers"))
    parser.add_argument("-o", "--out", type=os.path.abspath, required=True,
                        help=("output FASTA of adapters detected"))
    parser.add_argument("-N", "--number-of-reads", type=int, default=200000,
                        help="Number of reads to search in each fastq")
    parser.add_argument("-k", "--end-klength", type=int, default=16,
                        help=("Length of end kmer to tabulate for possible "
                              "adapter matches."))
    parser.add_argument("-E", "--end-min-match", type=float, default=0.0001,
                        help=("Minimum proportion of read match required to "
                              "report the endmer as a possible match."))
    parser.add_argument("-M", "--known-min-match", type=float, default=-1,
                        help=("Minimum proportion of read match required to "
                              "report the endmer as a possible match."
                              "Set to -1 (default) to accept all matches"))
    parser.add_argument("--quiet", action="store_true",
                        help="Suppress progress messages")
    parser.add_argument("--version", action="version", version="2017-06-21",
                        help="Display software version")
    return parser


def main(arguments=None):
    """Main method"""
    arguments = sys.argv[1:] if arguments is None else arguments
    # time0 = time()
    parser = generate_argparser()
    args = parser.parse_args(args=arguments)
    args.end_min_match = int(args.end_min_match * args.number_of_reads)
    args.known_min_match = int(args.known_min_match * args.number_of_reads)
    # ===== BEGIN ITERATION =====
    ndex = 0
    paired_end = False
    knownfinder1 = AdaptFind()
    endmerfinder1 = Adaptamers()
    knownfinder2 = AdaptFind()
    endmerfinder2 = Adaptamers()
    if args.fq2 is not None:
        paired_end = True
        input_fq = zip(args.fq1, args.fq2)
    else:
        input_fq = zip(args.fq1, ['']*len(args.fq1))
    for (fq1, fq2) in input_fq:
        infq1 = gopen(fq1, 'rt')
        if paired_end is True:
            infq2 = gopen(fq2, 'rt')
        while 1:
            # ===== Read1 Filtering =====
            infq1.readline()
            line1_seq = infq1.readline().rstrip()
            infq1.readline()
            infq1.readline()
            if not line1_seq:
                break
            if args.mode in ("known", "both"):
                knownfinder1.find_known_adapter(line1_seq)
            if args.mode in ("endmer", "both"):
                endmerfinder1.add_read_end(
                    line1_seq, k=args.end_klength)
            # ===== Read2 Filtering =====
            if paired_end is True:
                infq2.readline()
                line2_seq = infq2.readline()
                infq2.readline()
                infq2.readline()
                if not line2_seq:
                    break
                if args.mode in ("known", "both"):
                    knownfinder2.find_known_adapter(line2_seq)
                if args.mode in ("endmer", "both"):
                    endmerfinder2.add_read_end(
                        line2_seq, k=args.end_klength)
            ndex += 1
            if ndex % 10000 == 0:
                print(ndex, 'reads read')
            if ndex == args.number_of_reads:
                break
        infq1.close()
        if paired_end:
            infq2.close()
    adapter_entries = []
    for i, kadapt, ematch in (
            (1, knownfinder1, endmerfinder1),
            (2, knownfinder2, endmerfinder2)):
        if len(kadapt.found) > 0:
            print("=== Known adapters found in FQ{} ===".format(i))
            for entry, val in kadapt.found.items():
                print("{} found {} times.".format(
                    entry, val))
                if args.known_min_match == -1 or val >= args.known_min_match:
                    print("Min matches met, added to output.")
                    adapter_entries.append(entry)
            for bcode in kadapt.barcodes:
                print("Possible barcode found: {}".format(bcode))
        if len(ematch.endmers) > 0:
            print("=== Possible adapters inferred from FQ{} ===".format(i))
            for entry, val in ematch.endmers.items():
                print("{} found {} times.".format(entry, val))
                if val >= args.end_min_match:
                    print("Min matches met, added to output.")
                    adapter_entries.append(entry)
    final_entries = []
    for entry in adapter_entries:
        if not any(entry in x for x in adapter_entries if x != entry):
            final_entries.append(entry)
    with open(args.out, 'w') as outfile:
        for i, entry in enumerate(final_entries):
            outfile.write(">ADAPTER{}\n{}\n".format(i + 1, entry))
            outfile.write(">ADAPTER{}R\n{}\n".format(i + 1, complement(entry)))
    return ''


if __name__ == "__main__":
    main()
