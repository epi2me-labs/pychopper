# -*- coding: utf-8 -*-

import sys
import numpy as np
from pychopper import seq_utils as seu
from pychopper import hmmer_backend, edlib_backend
from pychopper.common_structures import Segment, Seq
from pychopper.alignment_hits import process_hits


def _build_segments(hits, config):
    "Build tuple of segments"
    segments = []
    if len(hits) == 0:
        return tuple(segments)
    for s in zip(hits, hits[1:]):
        strand = None
        seg_len = 0
        c = (s[0].Query, s[1].Query)
        if c in config:
            strand = config[c]
            seg_len = s[1].RefStart - s[0].RefEnd
        segments.append(Segment(s[0].RefStart, s[0].RefEnd, s[1].RefStart, s[1].RefEnd, strand, seg_len))
    return tuple(segments)


def analyse_hits(hits, config):
    """ Segment reads based on alignment hits using dynamic programming.
        The algorithm is based on the rule that each primer alignment hit can be used only once.
        Hence if a segment is included, the next one has to be excluded.
    """
    segments = _build_segments(hits, config)
    nr = len(segments)

    # No hits, return no segments:
    if nr == 0:
        return (), (), 0

    # Initialize DP matrix and traceback dictionary:
    M = np.zeros((2, nr), dtype=int)
    B = dict()

    # Initialize entry for the first segment:
    M[1, 0] = segments[0].Len

    # Fill in DP matrix:
    for j in range(1, nr):
        for i in range(2):
            if i == 0:
                # First row holds excluded segments.
                # The can transition from eiter exluded or included segments:
                M[i, j] = M[0, j - 1]
                B[i, j] = (0, j - 1)
                if M[1, j - 1] > M[0, j - 1]:
                    M[i, j] = M[1, j - 1]
                    B[i, j] = (1, j - 1)
            elif i == 1:
                # Included segments can only transition from previosuly excluded segments:
                M[i, j] = M[0, j - 1] + segments[j].Len
                B[i, j] = (0, j - 1)

    tlen = np.argmax(M[:, nr - 1])
    valid_segments = []

    # Traceback and build solution:
    p = (tlen, nr - 1)
    while True:
        if p[0] == 1:
            s = segments[p[1]]
            # Filter out invalid segments which can be part of the
            # solution:
            if s.Len > 0:
                valid_segments.append(s)
        if p[1] == 0:
            break
        p = B[p]

    return tuple(valid_segments), hits, tlen


def segments_to_reads(read, segments, keep_primers, bam_tags, detect_umis):
    "Convert segments to output reads with annotation"
    for s in segments:
        Start = s.Start
        End = s.End
        if keep_primers:
            Start = s.Left
            End = s.Right

        # Format FASTQ name and comment
        sr_id = "{}:{}|".format(Start, End)
        if bam_tags:
            try:
                name, comment = read.Name.split(" ", 1)
            except ValueError:
                name = read.Name
                comment = ""
            sr_name = sr_id + name + " CO:Z:" + comment + "\tTS:A:{}".format(s.Strand)
        else:
            sr_name = sr_id + read.Name + " strand=" + s.Strand

        if len(segments) > 1:
            sr_name += " rescue=1"

        umi = None
        if detect_umis:
            max_umi_ed = 3
            # Detect UMIs
            # Get adapters for UMI search
            padding = 40
            p1_from = max(0, s.Left - padding)
            p1_to = min(len(read.Seq), s.Start + padding)
            p_1 = read.Seq[p1_from:p1_to]
            p2_from = max(0, s.End - padding)
            p2_to = min(len(read.Seq), s.Right + padding)
            p_2 = read.Seq[p2_from:p2_to]

            umi, _ = edlib_backend.find_umi_single(
                [p_1 + 'NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN' + p_2, max_umi_ed])
            if bam_tags:
                sr_name += "\tRX:Z:{}".format(umi)
            else:
                sr_name += " umi={}".format(umi)

        sr_seq = read.Seq[Start:End]
        sr = Seq(sr_id + read.Id, sr_name, sr_seq, read.Qual[Start:End] if read.Qual is not None else None, umi)
        if s.Strand == '-':
            sr = seu.revcomp_seq(sr)
        yield sr


def chopper_phmm(reads, phmm_file, config, cutoff, threads, pool, min_batch):
    "Segment using the profile HMM backend"
    batch_hits = hmmer_backend.find_locations(reads, phmm_file, E=cutoff, pool=pool, min_batch=min_batch)
    for i, hits in enumerate(batch_hits):
        hits = process_hits(hits, cutoff)
        yield reads[i], analyse_hits(hits, config)


def chopper_edlib(reads, primers, config, max_ed, cutoff, pool, min_batch):
    "Segment using the edlib/parasail backend"
    batch_hits = edlib_backend.find_locations(reads, primers, max_ed=max_ed, pool=pool, min_batch=min_batch)
    for i, hits in enumerate(batch_hits):
        hits = process_hits(hits, cutoff)
        yield reads[i], analyse_hits(hits, config)
