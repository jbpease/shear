#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
SHEAR is a read trimmer that coordinates the automatic
finding of adapter sequences, removes adapters using Scythe,
implements various trimming and filtering options
for high-throughput short read sequences, and allows coordinated
removal of paired end sequence files.
"""

import sys
import os
import argparse
import subprocess
import gzip
from time import time
from datetime import datetime
from math import log
from adapt import main as makeadapt

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

_VERSION = "2017-06-21"

STAT_STRING = ("""
original reads:\t{}
original bases:\t{}
reads with adapters removed:\t{}
reads partially trimmed for quality:\t{}\t{}
bases trimmed for quality:\t{}\t{}
reads filtered for length:\t{}\t{}
bases filtered for length:\t{}\t{}
reads filtered for low avg quality:\t{}\t{}
bases filtered for low avg quality:\t{}\t{}
reads filtered for ambiguity:\t{}\t{}
bases filtered for ambiguity:\t{}\t{}
reads filtered for low info:\t{}\t{}
bases filtered for low info:\t{}\t{}
reads filtered because unpaired:\t{}\t{}
bases filtered because unpaired:\t{}\t{}
reads removed:\t{}\t{}
reads trimmed:\t{}\t{}
bases removed:\t{}\t{}""")


class TrimStat(object):
    """Trimming Statistics Handler"""

    def __init__(self):
        # Quartets represent [FQ1reads, FQ1bases, FQ2reads, FQ2bases]
        self.stats = dict(tuple([x, [0, 0, 0, 0]]) for x in [
            'scythe', 'total', 'ambig', 'totalqual',
            'trim', 'lowinfo', 'unpaired', 'tooshort'])

    def p1add(self, statname, value):
        """Add values to p1 stats"""
        if statname is False:
            return ''
        self.stats[statname][0] += 1
        self.stats[statname][1] += value
        return ''

    def p2add(self, statname, value):
        """Add values tp p2 stats"""
        if statname is False:
            return ''
        self.stats[statname][2] += 1
        self.stats[statname][3] += value
        return ''

    def p1scythe(self, value):
        """Add to scythe stat for p1"""
        self.stats['scythe'][0] += value

    def p2scythe(self, value):
        """Add to scythe stat for p2"""
        self.stats['scythe'][2] += value

    def stat_list(self, pair1=False, pair2=False, ndigits=3):
        """Create statistics output strings"""
        reads_rem = 0
        bases_rem = 0
        reads_trim = 0
        (i, j) = pair1 and (0, 1) or pair2 and (2, 3)
        xlist = [self.stats['total'][i],
                 self.stats['total'][j],
                 self.stats['scythe'][i]]

        for statname in ['trim', 'tooshort', 'totalqual', 'ambig',
                         'lowinfo', 'unpaired']:
            xlist.extend([self.stats[statname][i],
                          zdiv(self.stats[statname][i],
                               self.stats['total'][i], ndigits),
                          self.stats[statname][j],
                          zdiv(self.stats[statname][j],
                               self.stats['total'][j], ndigits)])
            if statname in ['trim', 'scythe']:
                reads_trim += self.stats[statname][i]
            else:
                reads_rem += self.stats[statname][i]
            bases_rem += self.stats[statname][j]
        xlist.extend([reads_rem, zdiv(reads_rem,
                                      self.stats['total'][i], ndigits),
                      reads_trim, zdiv(reads_trim,
                                       self.stats['total'][i], ndigits),
                      bases_rem, zdiv(bases_rem,
                                      self.stats['total'][j], ndigits)])
        return STAT_STRING.format(*xlist)


def timestamper():
    """Returns dash-separated timestamp string"""
    return datetime.now().strftime('%Y-%m-%d-%H-%M-%S-%f')


def read_config(configfilepath):
    """Reads the config file and returns command line arguments
    """
    args = []
    with open(configfilepath, 'r') as cfile:
        for line in [l.strip() for l in cfile.readlines()]:
            if len(line) < 1 or line[0] == '#':
                continue
            if '=' in line:
                elems = [d.strip() for d in line.split('=')]
            else:
                elems = [d.strip() for d in line.split()]
            args.extend(elems)
    return args


def gopen(path, mode):
    """Automatically invokes gzip file opening when path ends with .gz"""
    return path.endswith('.gz') and gzip.open(path, mode) or open(path, mode)


def zdiv(num, den, ndigits=-1):
    """Calculate a fraction, but avoid zero denominator errors"""
    if den == 0:
        return 0.0
    else:
        num = float(num) / float(den)
        if ndigits != -1:
            num = round(num, ndigits)
        return num


def mutual_info(seq):
    """Calculate a simple dinucleotide mutual information score"""
    if len(seq) < 3:
        return 0.0
    dinucs = {}
    monucs = {}
    for j in range(len(seq) - 1):
        dinucs[seq[j:j+2]] = dinucs.get(seq[j:j+2], 0) + 1
        monucs[seq[j]] = monucs.get(seq[j], 0) + 1
    monucs[seq[-1]] = monucs.get(seq[-1], 0) + 1
    total_monucs = float(sum(monucs.values()))
    total_dinucs = float(sum(dinucs.values()))
    mutual_information = 0
    if total_dinucs == 0:
        return 0.0
    for dinuc in dinucs:
        p_dinuc = float(dinucs[dinuc]) / total_dinucs
        mutual_information += (p_dinuc * log(p_dinuc / (
            (float(monucs[dinuc[0]]) / total_monucs) *
            (float(monucs[dinuc[1]]) / total_monucs))))
    return mutual_information


def run_scythe(fastqpath, params, execpath='scythe', logfilepath='scythe.log',
               quiet=False):
    """Runs the Scythe Bayesian Trimmer"""
    logx = open(logfilepath, 'w')
    cmd = "{} {} {}".format(execpath,
                            ' '.join(["{} {}".format(k, v)
                                      for k, v in params.items()]),
                            fastqpath)
    print(cmd)
    if quiet is False:
        print("executing subprocess \"{}\"".format(cmd))
    process = subprocess.Popen(
        cmd, stdout=logx, stderr=subprocess.STDOUT, shell=True)
    process.communicate()
    logx.close()
    contaminated = 0
    scythelog = gopen(logfilepath, 'rt')
    line = scythelog.readline()
    while line:
        if line.startswith('contaminated'):
            contaminated = int(line[line.find(': ') + 2:line.find(', ')])
            break
        line = scythelog.readline()
    scythelog.close()
    return contaminated


def detect_fastq_format(filepath, ntest=200000):
    """Determine FASTQ Format"""
    i = 0
    with gopen(filepath, 'rt') as fastqfile:
        code = ''
        while code in ['', 'phred33', 'phred64'] and i < ntest:
            i += 1
            fastqfile.readline()  # Skip Header
            bases = fastqfile.readline()
            fastqfile.readline()  # Skip spacer
            quals = fastqfile.readline()
            bases = bases.rstrip()
            if code == 'phred33':
                code = "sanger"
            elif code == 'phred64':
                if any(x in quals for x in ';<=>?'):
                    code = "solexa"
                elif any(x in quals for x in '@"'):
                    code = 'illumina'
            else:
                if any(x in quals for x in '!"#$%&\'()*+,-./0123456789:'):
                    code = "phred33"
                elif any(x in quals for x in 'KLMNOPQRSTUVWXYZ[]^_`abcdefgh'):
                    code = 'phred64'
                    if any(x in quals for x in ';<=>?'):
                        code = "solexa"
                    elif any(x in quals for x in '@"'):
                        code = 'illumina'
        while i < ntest:
            fastqfile.readline()
            bases = fastqfile.readline()
            fastqfile.readline()
            fastqfile.readline()
            if bases == '':
                break
            bases = bases.rstrip()
    if code == 'phred64':
        code = "illumina"
    elif code == 'phred33':
        code = "sanger"
    return code


def colon_sep_check(arg, name):
    try:
        return tuple([int(x) for x in arg.split(':')][0:2])
    except:
        msg = "{} must be in form FRONT:END".format(name)
        raise SyntaxError(msg)


def generate_argparser():
    parser = argparse.ArgumentParser(
        prog="shear.py",
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
    parser.add_argument("--out1", required=True, type=os.path.abspath,
                        help=("Output fastq file path. Note this is a single "
                              "output file that concatenates the processed "
                              "outputs from all files in --fq1"))
    parser.add_argument("--out2", type=os.path.abspath,
                        help=("output fastq file path. Note this is a single "
                              "output file that concatenates the processed "
                              "outputs from all files in --fq2"))
    parser.add_argument("--filt1", type=os.path.abspath,
                        help=("Output fastq file path for sequences that were"
                              "filtered outs. "
                              "Note this is a single "
                              "output file that concatenates the rejected"
                              "outputs from all files in --fq1."
                              "Default is [--out1]_filtered_1.fastq"))
    parser.add_argument("--filt2", type=os.path.abspath,
                        help=("Output fastq file path for sequences that were "
                              "filtered out. Note this is a single "
                              "output file that concatenates the rejected "
                              "outputs from all files in --fq2. "
                              "Default is [--out2]_filtered_2.fastq"))
    parser.add_argument("-a", "--adapters", type=os.path.abspath, nargs='*',
                        help=("Skip adapter finding an use these "
                              "adapter files. Either enter (1) a single "
                              " file to use for all fastq files, (2) "
                              " one file per single end file or pair "
                              " of paired-end files. "
                              "Adapter file(s) should be in FASTA format."))
    parser.add_argument("-t", "--platform",
                        choices=('TruSeq', 'TruSeqDualIndex'),
                        default='TruSeq', help="Sequencing Platform")
    parser.add_argument("-f", "--trim-fixed", default="0:0",
                        help=("Trim a fixed number of bases "
                              "from the FRONT:END of the sequence "
                              "(NOT recommended)."))
    parser.add_argument("-q", "--trim-qual", default="20:20",
                        help=("Trim bases below this quality score from "
                              "the FRONT:END of each read."))
    parser.add_argument("--trim-qual-pad", default="0:0",
                        help=("Trim additional bases next to low-quality "
                              "bases specified by --trim-qual "
                              "from the FRONT:END of each read."))
    parser.add_argument("-p", "--trim-poly", type=int, default=12,
                        help=("Trim poly-A or poly-T repeats of "
                              "at least this length from the front or end."))
    parser.add_argument("-n", "--retain-ambig", action="store_true",
                        help=("By default ambiguous nucleotides (N) are "
                              "removed from both ends of each read. If this "
                              "flag is specified, N's are retained."))
    parser.add_argument("-L", "--filter-length", type=int, default=30,
                        help=("Filter out reads that contain fewer than this "
                              "many characters after trimming."))
    parser.add_argument("-Q", "--filter-quality", type=int, default=3,
                        help=("Filter out reads with a mean quality score "
                              "below this value (before trimming)."))
    parser.add_argument("-U", "--filter-unpaired", action="store_true",
                        help=("If either read in a read pair is filtered out, "
                              "the counterpart reads is also filtered out "
                              "regardless of quality."))
    parser.add_argument("-I", "--filter-low-info", type=float, default=0.,
                        help=("Filter out reads with mutual information "
                              "scores exceeding this value "
                              "(ADVANCED, removes highly repetitive reads)."))
    parser.add_argument("-A", "--filter-ambig", type=int, default=5,
                        help=("Filter reads with more than this number"
                              "of ambiguous nucleotides "
                              "(N's; set as 0 to skip"))
    parser.add_argument("-z", "--trim-pattern-3", nargs=1,
                        help="Comma-separated list of specific sequences to "
                             "trim from the 3' end. Can be used for extra "
                             "stringent adapter trimming.")
#    parser.add_argument("-y", "--trim-pattern-5", nargs=1,
#                        help="Comma-separated list of specific sequences to "
#                             "trim from the 5' end. (Not recommmended).")
    parser.add_argument("--scythe-skip", action="store_true",
                        help="Skip scythe 3' adapter removal.")
    parser.add_argument("-X", "--scythe-executable", default="scythe",
                        help="Set the path of the scythe executable manually.")
    parser.add_argument("-s", "--scythe-prior", type=float, default=0.1,
                        help=("Bayesian prior for proportion of adapters "
                              "expected to be sampled in Scythe."))
    parser.add_argument("--scythe-match", type=int, default=5,
                        help=("Minimum number of bases required for a "
                              "match in Scythe."))
    parser.add_argument("-m", "--adapter-mode",
                        choices=("known", "endmer", "both"),
                        default="known",
                        help=("(Adapter finding) "
                              "known=only use list of known adapters;"
                              "endmer=search for common 3'end sequences;"
                              "both=both known and endmers"))
    parser.add_argument("-N", "--adapter-number-of-reads",
                        type=int, default=200000,
                        help=("(Adapter finding) "
                              "Number of reads to search in each fastq"))
    parser.add_argument("-k", "--adapter-end-klength", type=int, default=16,
                        help=("(Adapter finding) Length of 3'-end kmer "
                              "to tabulate for possible adapter matches."))
    parser.add_argument("-E", "--adapter-end-min-match", type=float,
                        default=0.0001,
                        help=("(Adapter finding) "
                              "Minimum proportion of read matches required to "
                              "report the 3'-end-mer as a possible match."))
    parser.add_argument("-M", "--adapter-known-min-match", type=float,
                        default=-1,
                        help=("(Adapter finding) "
                              "Minimum proportion of read matches required to "
                              "report a known contaminant as a "
                              "possible match. "
                              "Use -1 (default) to accept all matches."))
    parser.add_argument("--log-path", type=os.path.abspath,
                        help=("Manually specify log file path, "
                              "default is 'shear_TIMESTAMP'"))
    parser.add_argument("--retain-temp",
                        choices=["none", "tempfastq",
                                 "exceptadapters", "all"],
                        default="none",
                        help=("Retain temporary files "
                              "(none=remove all; "
                              "tempfastq=remove temporary fastq files "
                              "from scythe; "
                              "exceptadapters=remove temp fastq and log "
                              "files from Scythe but keep adapter file; "
                              "all=keep all temporary file)"))
    parser.add_argument("--temp-dir", default='.', type=os.path.abspath,
                        help="directory to use for temporary files")
    parser.add_argument("--clean-header", action="store_true",
                        help=("removes any additional terms from the"
                              "header line (useful after STAR)"))
    parser.add_argument("--quality-scale", choices=(
                            'sanger', 'illumina', 'phred', 'solexa'),
                        help=("Quality scale is usually automatically "
                              "determined, but use this to set manually."))
    parser.add_argument("--quiet", action="store_true",
                        help="Suppress progress messages")
    parser.add_argument("--version", action="version", version=_VERSION,
                        help="Display software version")
    return parser


def main(arguments=None):
    """Main method"""
    arguments = sys.argv[1:] if arguments is None else arguments
    if (len(arguments) == 1 and arguments[0] not in
            ('-h', '--help', '--version')):
        arguments = read_config(arguments[0])
        print("Config file used.")
        print("Executing with arguments: ", " ".join(arguments))
    time0 = time()
    parser = generate_argparser()
    args = parser.parse_args(args=arguments)

    # ===== Establish Log and Temporary Directory =====
    timestamp = timestamper()
    if not os.path.exists(args.temp_dir):
        os.mkdir(args.temp_dir)
    if args.log_path is None:
        args.log_path = os.path.abspath("shear_" + timestamp + ".log")
    mainlog = gopen(args.log_path, 'wt')
    mainlog.write("shear version={}\n{}\n".format(_VERSION,
                                                  " ".join(sys.argv)))
    # ===== Additional Argument Checks =====
    trimfixed = colon_sep_check(args.trim_fixed, "--trim-fixed")
    args.trim_pattern_3 = (None if args.trim_pattern_3 is None
                           else args.trim_pattern_3[0].split(','))
    trimqual = colon_sep_check(args.trim_qual, "--trim-qual")
    trimqualpad = colon_sep_check(args.trim_qual_pad, "--trim-qual-pad")
    if args.fq2 is not None:
        if len(args.fq1) != len(args.fq2):
            raise RuntimeError(
                "Enter an equal number of files in --fq1 and --fq2")
    polynucs = []
    if not args.retain_ambig:
        polynucs.append(('N', 1))
    if args.trim_poly:
        polynucs.extend([('A', args.trim_poly), ('T', args.trim_poly)])
    if not args.scythe_skip and any(
            x.endswith('.gz') for x in args.fq1 + args.fq2):
        msg = "If Scythe is used inputs cannot be GZIPed"
        mainlog.write(msg + "\n")
        raise SyntaxError(msg)
    # Adapter file setup
    if args.adapters is None:
        args.adapters = tuple([None] * len(args.fq1))
    elif len(args.adapters) == 1:
        args.adapters = tuple([args.adapters[0]] * len(args.fq1))
    elif len(args.adapters) != len(args.fq1):
        msg = ("--adapters must either contain "
               "one file, or a number of files the same as "
               "--fq1")
        mainlog.write(msg + "\n")
        raise SyntaxError(msg)
    # ===== Establish file paths and paired/single mode =====
    args.fq1 = [os.path.abspath(x) for x in args.fq1]
    args.out1 = os.path.abspath(args.out1)
    outfq1 = gopen(args.out1, 'wt')
    filterfq1 = gopen(os.path.abspath("{}_filtered_1.fastq".format(
        args.filt1 if args.filt1 is not None else args.out1)), 'wt')
    paired_end = False
    if args.fq2:
        paired_end = True
        if not args.out2:
            msg = "--fq2 specified without --out2"
            mainlog.write(msg + "\n")
            raise RuntimeError(msg)
        args.out2 = os.path.abspath(args.out2)
        filterfq2 = gopen(os.path.abspath("{}_filtered_2.fastq".format(
            args.filt2 if args.filt2 is not None else args.out2)), 'wt')
        outfq2 = gopen(args.out2, 'wt')
        input_fq = zip(args.fq1, args.fq2, args.adapters)
    else:
        input_fq = zip(args.fq1, ['']*len(args.fq1),
                       args.adapters[0])
    # ===== Establish stats container =====
    stats = TrimStat()
    # ===== Scythe Adapter Preparation =====
    if args.scythe_skip:
        mainlog.write("Skipping Scythe...\n")
    else:
        mainlog.write("Using temporary working directory: {}\n".format(
            args.temp_dir))

    # ===== BEGIN ITERATION =====
    ndex = 0
    adapterfiles = []
    for (fq1, fq2, adapterfp) in input_fq:
        if args.scythe_skip is False:
            if adapterfp is None:
                adapterfp = os.path.join(
                    args.temp_dir,
                    "scythe_adapters_{}.{}.fasta".format(
                        timestamp, ndex))
                makeadaptargs = [str(x) for x in [
                    "--fq1", fq1, "--fq2", fq2,
                    "--out", adapterfp,
                    "-m", args.adapter_mode,
                    "-N", args.adapter_number_of_reads,
                    "-E", args.adapter_end_min_match,
                    "-M", args.adapter_known_min_match,
                    "-k", args.adapter_end_klength]]
                if paired_end:
                    makeadaptargs.extend(["--fq2", fq2])
                makeadapt(arguments=makeadaptargs)
                mainlog.write("Adapters written to: {}".format(adapterfp))
            adapterfiles.append(adapterfp)
        code_p1 = (detect_fastq_format(fq1,)
                   if args.quality_scale is None else
                   args.quality_scale)
        if not code_p1:
            msg = "fq1 quality score scale cannot be determined"
            mainlog.write(msg + "\n")
            raise RuntimeError(msg)
        else:
            mainlog.write("fq1 quality scale detected: {}\n".format(code_p1))
        if paired_end:
            code_p2 = (detect_fastq_format(fq1,)
                       if args.quality_scale is None else
                       args.quality_scale)
            if not code_p2:
                msg = "fq2 quality score scale cannot be determined"
                mainlog.write(msg + "\n")
                raise RuntimeError(msg)
            else:
                mainlog.write("fq2 quality scale detected: {}\n".format(
                    code_p2))
            if code_p1 != code_p2:
                msg = ("fq1 and fq2 quality score scales do not match "
                       "({}, {})").format(code_p1, code_p2)
                mainlog.write(msg + "\n")
                raise RuntimeError(msg)
        # ===== Scythe Adapter Trimming on Pair 1 =====
        if not args.scythe_skip:
            tempfile1 = os.path.abspath(os.path.join(
                args.temp_dir, "scythetemp_p1_{}_{}.fastq".format(
                    timestamp, ndex)))
            logfile1 = os.path.abspath(os.path.join(
                args.temp_dir, "scythetemp_p1_{}_{}.log".format(
                    timestamp, ndex)))
            mainlog.write("Running Scythe on fq1...\n")
            scythe_params = {'-a': adapterfp,
                             '-q': code_p1,
                             '-p': args.scythe_prior,
                             '-n': args.scythe_match,
                             '-o': tempfile1}
            stats.p1scythe(run_scythe(fq1, scythe_params, logfilepath=logfile1,
                                      execpath=args.scythe_executable,
                                      quiet=args.quiet))
            mainlog.write("Scythe on fq1 [{}] complete at {} seconds. "
                          "{} reads had adapters removed. ".format(
                            fq1, round(time() - time0, 2),
                            stats.stats['scythe'][0]))
        # ===== Scythe Adapter Trimming on Pair 2 =====
        if paired_end and not args.scythe_skip:
            tempfile2 = os.path.abspath(os.path.join(
                args.temp_dir, "scythetemp_p2_{}_{}.fastq".format(
                    timestamp, ndex)))
            logfile2 = os.path.abspath(os.path.join(
                args.temp_dir, "scythetemp_p2_{}_{}.log".format(
                    timestamp, ndex)))
            mainlog.write("Running Scythe on fq2...\n")
            scythe_params = {'-a': adapterfp,
                             '-q': code_p2,
                             '-p': args.scythe_prior,
                             '-n': args.scythe_match,
                             '-o': tempfile2}
            stats.p2scythe(run_scythe(fq2, scythe_params, logfilepath=logfile2,
                                      execpath=args.scythe_executable,
                                      quiet=args.quiet))
            mainlog.write("Scythe on fq2 [{}] complete at {} seconds. "
                          "{} reads had adapters removed. ".format(
                              fq2, round(time() - time0, 2),
                              stats.stats['scythe'][2]))
        # ===== Quality Trimming Setup =====
        lowq_front = set([])
        lowq_end = set([])
        qual_offset = 32
        if code_p1 in ['solexa', 'illumina_1.3', 'illumina_1.5']:
            qual_offset = 64
        if trimqual[0]:
            lowq_front = set([chr(x) for x in range(
                qual_offset, qual_offset + trimqual[0])])
        if trimqual[1]:
            lowq_end = set([chr(x) for x in range(qual_offset,
                                                  qual_offset + trimqual[1])])
        # ===== Open Files and Establish Trimming Counters =====
        if args.scythe_skip is True:
            infq1 = gopen(fq1, 'r')
        else:
            infq1 = gopen(tempfile1, 'r')
        if paired_end is True:
            if args.scythe_skip is True:
                infq2 = gopen(fq2, 'r')
            else:
                infq2 = gopen(tempfile2, 'r')
        while 1:
            read1_removed = False
            read2_removed = False
            # ===== Read1 Filtering =====
            line1_header = infq1.readline().rstrip()
            line1_seq = infq1.readline().rstrip()
            infq1.readline()
            line1_qual = infq1.readline()
            if not line1_qual:
                break
            if args.clean_header:
                line1_header = line1_header.split()[0]
            line1_qual = line1_qual.rstrip()
            len_read1 = len(line1_qual)
            stats.p1add('total', len_read1)
            if len_read1 < args.filter_length:
                read1_removed = 'tooshort'
            if read1_removed is False and args.filter_ambig > 0:
                if line1_seq.count('N') > args.filter_ambig:
                    read1_removed = 'ambig'
            if read1_removed is False and args.filter_quality > 0:
                if sum([ord(x) - qual_offset for x in line1_qual])/(
                        len_read1) < args.filter_quality:
                    read1_removed = 'totalqual'
            if read1_removed is False and args.filter_low_info > 0:
                mutinfo = mutual_info(line1_seq)
                if mutinfo >= args.filter_low_info:
                    read1_removed = 'lowinfo'
            # ===== Read2 Filtering =====
            if paired_end:
                line2_header = infq2.readline().rstrip()
                line2_seq = infq2.readline().rstrip()
                infq2.readline()
                line2_qual = infq2.readline()
                if not line2_qual:
                    break
                if args.clean_header:
                    line2_header = line2_header.split()[0]
                line2_qual = line2_qual.rstrip()
                len_read2 = len(line2_qual)
                stats.p2add('total', len_read2)
                if (read1_removed is not False and
                        args.filter_unpaired is True):
                    read2_removed = 'unpaired'
                elif len_read2 < args.filter_length:
                    read2_removed = 'tooshort'
                if read2_removed is False and args.filter_ambig > 0:
                    if line2_seq.count('N') > args.filter_ambig:
                        read2_removed = 'ambig'
                if read2_removed is False and args.filter_quality > 0:
                    if sum([ord(x) - qual_offset for x in line2_qual])/(
                            len_read2) < args.filter_quality:
                        read2_removed = 'totalqual'
                if read2_removed is False and args.filter_low_info > 0:
                    mutinfo = mutual_info(line2_seq)
                    if mutinfo >= args.filter_low_info:
                        read2_removed = 'lowinfo'
            if (read2_removed is not False and
                    args.filter_unpaired is True and paired_end is True):
                if read1_removed is False:
                    read1_removed = 'unpaired'
            # ===== Read1 trimming from 3' end =====
            if read1_removed is False:
                coord_end1 = len_read1 - trimfixed[1] - 1
                if args.trim_pattern_3 is not None:
                    for elem in args.trim_pattern_3:
                        if line1_seq.endswith(elem):
                            coord_end1 = len_read1 - (
                                max(len(elem), trimfixed[1]) - 1)
                            break
                for nuc, num in polynucs:
                    if line1_seq.endswith(nuc * num):
                        j = len(line1_seq) - num
                        while line1_seq[j] == nuc and j:
                            j -= 1
                        coord_end1 = j + 1
                while read1_removed is False:
                    if coord_end1 <= args.filter_length + trimfixed[0]:
                        read1_removed = 'totalqual'
                        break
                    if line1_qual[coord_end1] not in lowq_end:
                        if not (set(line1_qual[coord_end1 - trimqualpad[1]:
                                               coord_end1:-1]) & lowq_end):
                            break
                    coord_end1 -= 1
                # ===== Read1 trimming for 5' end =====
                coord_start1 = trimfixed[0] + 0
                for nuc, num in polynucs:
                    if line1_seq.startswith(nuc * num):
                        j = num + 0
                        while line1_seq[j] == nuc:
                            j += 1
                            if j >= len(line1_seq):
                                break
                        coord_start1 = j
                while not read1_removed:
                    if coord_start1 >= coord_end1:
                        read1_removed = 'totalqual'
                        break
                    if line1_qual[coord_start1] not in lowq_front:
                        if not (lowq_front & set(line1_qual[
                                coord_start1:coord_start1 + trimqualpad[0]])):
                            break
                    coord_start1 += 1
            # ===== Read2 trimming of 3' end =====
            if read2_removed is False and paired_end is True:
                coord_end2 = len_read2 - trimfixed[1] - 1
                if args.trim_pattern_3:
                    for elem in args.trim_pattern_3:
                        if line2_seq.endswith(elem):
                            coord_end2 = len_read2 - (
                                max(len(elem), trimfixed[1]) - 1)
                            break
                for nuc, num in polynucs:
                    if line2_seq.endswith(nuc * num):
                        j = len(line2_seq) - num
                        while line2_seq[j] == nuc:
                            j -= 1
                            if j == 0:
                                break
                        coord_end2 = j + 1
                while read2_removed is False:
                    if coord_end2 <= args.filter_length + trimfixed[0]:
                        read2_removed = 'totalqual'
                        break
                    if line2_qual[coord_end2] not in lowq_end:
                        if not (lowq_end &
                                set(line2_qual[coord_end2 - trimqualpad[1]:
                                               coord_end2])):
                            break
                    coord_end2 -= 1
                # ===== Read2 trimming of 5' end =====
                coord_start2 = trimfixed[0] + 0
                for nuc, num in polynucs:
                    if line2_seq.startswith(nuc * num):
                        j = num + 0
                        while line2_seq[j] == nuc:
                            j += 1
                            if j >= len(line2_seq):
                                break
                        coord_start2 = j
                while read2_removed is False:
                    if coord_start2 >= len_read2 or coord_start2 >= coord_end2:
                        read2_removed = 'totalqual'
                        break
                    if line2_qual[coord_start2] not in lowq_front:
                        if not (lowq_front & set(line2_qual[
                                coord_start2: coord_start2 + trimqualpad[0]])):
                            break
                    coord_start2 += 1
                if read2_removed is False and (
                        coord_end2 - coord_start2 < args.filter_length):
                    read2_removed = 'tooshort'
            if read1_removed is not False and (
                    args.filter_unpaired is True and paired_end is True):
                if read2_removed is False:
                    read2_removed = 'unpaired'
            if read2_removed and args.filter_unpaired is True and (
                    paired_end is True):
                if read1_removed is False:
                    read1_removed = 'unpaired'
            if read1_removed is False:
                outfq1.write("{}\n{}\n+\n{}\n".format(
                    line1_header,
                    line1_seq[coord_start1:coord_end1 + 1],
                    line1_qual[coord_start1:coord_end1 + 1]))
                ntrim = len_read1 - (coord_end1 - coord_start1)
                if ntrim:
                    stats.p1add('trim', ntrim)
            else:
                stats.p1add(read1_removed, len_read1)
                filterfq1.write("{}_{}\n{}\n+\n{}\n".format(
                    line1_header, read1_removed, line1_seq, line1_qual))
            if paired_end is True:
                if read2_removed is not False:
                    stats.p2add(read2_removed, len_read2)
                    filterfq2.write("{}_{}\n{}\n+\n{}\n".format(
                        line2_header, read2_removed,
                        line2_seq, line2_qual))
                else:
                    outfq2.write("{}\n{}\n+\n{}\n".format(
                        line2_header,
                        line2_seq[coord_start2:coord_end2 + 1],
                        line2_qual[coord_start2:coord_end2 + 1]))
                    ntrim = len_read2 - (coord_end2 - coord_start2)
                    if ntrim:
                        stats.p2add('trim', ntrim)

        infq1.close()
        if paired_end:
            infq2.close()
        ndex += 1
    # FINISH ITERATION
    outfq1.close()
    filterfq1.close()
    if paired_end:
        outfq2.close()
        filterfq2.close()
    # ===== Write Statistics =====
    mainlog.write("\nfq1_stats\n{}\n\nfq2 stats\n{}\n".format(
        stats.stat_list(pair1=True), stats.stat_list(pair2=True)))
    # ===== Optional File Cleanup =====
    if args.retain_temp != "all":
        os.remove(tempfile1)
        if paired_end:
            os.remove(tempfile2)
        if args.retain_temp != "tempfastq":
            os.remove(logfile1)
            if paired_end:
                os.remove(logfile2)
            if args.retain_temp != "exceptadapters":
                if args.scythe_skip is False:
                    os.remove(adapterfp)
    # ===== Finish log and close =====
    mainlog.write("\nTotal running time: {} seconds.\n".format(
        round(time() - time0, 2)))
    mainlog.close()
    return ''


if __name__ == "__main__":
    main()
