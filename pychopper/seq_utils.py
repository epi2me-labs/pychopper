# -*- coding: utf-8 -*-

from six.moves import reduce
from numpy.random import random
from pysam import FastxFile
from math import log
import sys
from pychopper.common_structures import Seq

""" Utilities manipulating biological sequences and formats. Extensions to biopython functionality.
"""

# Reverse complements of bases, taken from dragonet:
comp = {
    'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'X': 'X', 'N': 'N',
    'a': 't', 't': 'a', 'c': 'g', 'g': 'c', 'x': 'x', 'n': 'n',
    '-': '-'
}


def base_complement(k):
    """ Return complement of base.

    Performs the subsitutions: A<=>T, C<=>G, X=>X for both upper and lower
    case. The return value is identical to the argument for all other values.

    :param k: A base.
    :returns: Complement of base.
    :rtype: str

    """
    try:
        return comp[k]
    except KeyError:
        sys.stderr.write(
            "WARNING: No reverse complement for {} found, returning argument.".format(k))
        return k


def reverse_complement(seq):
    """ Return reverse complement of a string (base) sequence.

    :param seq: Input sequence.
    :returns: Reverse complement of input sequence.
    :rtype: str

    """
    if len(seq) == 0:
        return seq
    return reduce(lambda x, y: x + y, map(base_complement, seq[::-1]))


def readfq(fastq, sample=None, min_qual=0, rfq_sup={}):  # this is a generator function
    """Read fastx files.

    This is a generator function that yields sequtils.Seq objects.
    Optionally filter by a minimum mean quality (min_qual).
    Optionally subsample the fastx file using sample (0.0 - 1.0)
    """
    sup = ("out_fq" in rfq_sup) and (rfq_sup["out_fq"] is not None)
    tsup = "total" in rfq_sup
    if sup:
        fh = open(rfq_sup["out_fq"], "w")

    with FastxFile(fastq) as fqin:
        for fx in fqin:
            if sample is None or (random() < sample):
                if tsup:
                    rfq_sup["total"] += 1
                qual_array = fx.get_quality_array()
                if qual_array is None or \
                        not (sum(qual_array) / len(qual_array)) < min_qual:
                    if tsup:
                        rfq_sup["pass"] += 1
                    yield Seq(
                        Id=fx.name,
                        Name=f"{fx.name} {fx.comment}" if fx.comment else fx.name,
                        Seq=fx.sequence, Qual=fx.quality, Umi=None)
                else:
                    if sup:
                        fh.write(f"{fx}\n")
    if sup:
        fh.flush()
        fh.close()


def writefq(r, fh):
    "Write read to fastq file"
    q = r.Qual
    if q is None:
        q = "!" * len(r.Seq)
    fh.write("@{}\n{}\n+\n{}\n".format(r.Name, r.Seq, q))


def revcomp_seq(seq):
    """ Reverse complement sequence record """
    qual = seq.Qual
    if qual is not None:
        qual = qual[::-1]
    return Seq(seq.Id, seq.Name, reverse_complement(seq.Seq), qual, seq.Umi)


def get_runid(desc):
    """ Parse out runid from sequence description. """
    tmp = [t for t in desc.split(" ") if t.startswith("runid")]
    if len(tmp) != 1:
        return "NA"
    return tmp[0].rsplit("=", 1)[1]


def record_size(read, in_format='fastq'):
    """ Calculate record size. """
    if read.Qual is None:
        in_format = 'fasta'
    dl = len(read.Name)
    sl = len(read.Seq)
    if in_format == 'fastq':
        bl = dl + 2 * sl + 6
    elif in_format == 'fasta':
        bl = dl + sl + 3
    else:
        raise Exception("Unkonwn format!")
    return bl


def get_primers(primers):
    "Load primers from fasta file"
    all_primers = {}
    for primer in readfq(primers):
        all_primers[primer.Name] = primer.Seq
        all_primers['-' + primer.Name] = reverse_complement(primer.Seq)
    return all_primers


def errs_tab(n):
    """Generate list of error rates for qualities less than equal than n."""
    return [10**(q / -10) for q in range(n + 1)]


def mean_qual(quals, qround=False, tab=errs_tab(128)):
    """Calculate average basecall quality of a read.
    Receive the ascii quality scores of a read and return the average quality for that read
    First convert Phred scores to probabilities,
    calculate average error probability
    convert average back to Phred scale
    """
    if quals:
        mq = -10 * log(sum([tab[ord(q) - 33] for q in quals]) / len(quals), 10)
        if qround:
            return round(mq)
        else:
            return mq
    else:
        return 0.0
