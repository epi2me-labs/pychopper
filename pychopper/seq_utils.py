# -*- coding: utf-8 -*-

from six.moves import reduce
from numpy.random import random
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


def readfq(fp, sample=None, min_qual=None, rfq_sup={}):  # this is a generator function
    """
    Below function taken from https://github.com/lh3/readfq/blob/master/readfq.py
    Much faster parsing of large files compared to Biopyhton.
    """
    sup = ("out_fq" in rfq_sup) and (rfq_sup["out_fq"] is not None)
    tsup = "total" in rfq_sup
    if sup:
        fh = open(rfq_sup["out_fq"], "w")
    last = None  # this is a buffer keeping the last unprocessed line
    while True:  # mimic closure; is it a bad idea?
        if not last:  # the first record or a record following a fastq
            for l in fp:  # search for the start of the next record
                if l[0] in '>@':  # fasta/q header line
                    last = l[:-1]  # save this line
                    break
        if not last:
            break
        name, seqs, last = last[1:], [], None
        for l in fp:  # read the sequence
            if l[0] in '@+>':
                last = l[:-1]
                break
            seqs.append(l[:-1])
        if not last or last[0] != '+':  # this is a fasta record
            if sample is None or (random() < sample):
                if tsup:
                    rfq_sup["total"] += 1
                yield Seq(name.split(" ", 1)[0], name, ''.join(seqs), None)  # yield a fasta record
            if not last:
                break
        else:  # this is a fastq record
            seq, leng, seqs = ''.join(seqs), 0, []
            for l in fp:  # read the quality
                seqs.append(l[:-1])
                leng += len(l) - 1
                if leng >= len(seq):  # have read enough quality
                    last = None
                    if sample is None or (random() < sample):
                        quals = "".join(seqs)
                        oseq = Seq(Id=name.split(" ", 1)[0], Name=name, Seq=seq, Qual=quals)
                        if tsup:
                            rfq_sup["total"] += 1
                        if not (min_qual is not None and min_qual > 0 and mean_qual(quals) < min_qual):
                            if tsup:
                                rfq_sup["pass"] += 1
                            yield oseq
                        else:
                            if sup:
                                writefq(oseq, fh)
                    break
            if last:  # reach EOF before reading enough quality
                if sample is None or (random() < sample):
                    yield Seq(name.split(" ", 1)[0], name, seq, None)  # yield a fasta record instead
                break

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
    return Seq(seq.Id, seq.Name, reverse_complement(seq.Seq), qual)


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
    for primer in readfq(open(primers, 'r')):
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
