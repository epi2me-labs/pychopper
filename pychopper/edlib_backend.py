# -*- coding: utf-8 -*-

import edlib

from pychopper import seq_utils as seu
from pychopper.common_structures import Hit
from pychopper import utils
from pychopper.parasail_backend import refine_locations


def find_locations(reads, all_primers, max_ed, pool, min_batch):
    "Find alignment hits of all primers in all reads using the edlib/parasail backend"
    for batch in utils.batch(reads, min_batch):
        for res in pool.map(_find_locations_single, zip(batch, [(all_primers, max_ed)] * len(batch))):
            try:
                yield res
            except StopIteration:
                return


def find_umi_single(params):
    "Find UMI in a single reads using the edlib/parasail backend"
    read = params[0]
    max_ed = params[1]

    pattern_list = [
        (
            "TTTVVVVTTVVVVTTVVVVTTVVVVTTT",
            "V",
            [("V", "A"), ("V", "G"), ("V", "C")],
            True
        ),
        (
            "AAABBBBAABBBBAABBBBAABBBBAAA",
            "B",
            [("B", "T"), ("B", "G"), ("B", "C")],
            False
        )
    ]

    best_ed = None
    best_pattern = None
    best_result = None
    for pattern, wildcard, equalities, forward in pattern_list:
        result = edlib.align(pattern,
                             read,
                             task="path",
                             mode="HW",
                             k=max_ed,
                             additionalEqualities=equalities)
        if result['editDistance'] == -1:
            continue
        if not best_ed or result['editDistance'] < best_ed:
            best_ed = result['editDistance']
            best_pattern = (pattern, wildcard, equalities, forward)
            best_result = result

    if best_ed is None:
        return None, None

    # Extract UMI
    pattern, wildcard, equalities, forward = best_pattern
    ed = best_result["editDistance"]

    locs = best_result["locations"][0]
    umi = read[locs[0]:locs[1]+1]
    if not forward:
        umi = seu.reverse_complement(umi)
    # Do not assign UMIs where UMI probe has aligned to the 'N' spacer.
    if 'N' in umi:
        return None, None
    return umi, ed


def _find_locations_single(params):
    "Find alignment hits of all primers in a single reads using the edlib/parasail backend"
    read = params[0]
    all_primers, max_ed = params[1]
    all_locations = []
    for primer_acc, primer_seq in all_primers.items():
        primer_max_ed = int(max_ed * len(primer_seq))
        result = edlib.align(primer_seq, read.Seq,
                             mode="HW", task="locations", k=primer_max_ed)
        ed = result["editDistance"]
        locations = result["locations"]
        if locations:
            # all_locations[primer_acc] = []
            for refstart, refend in locations:
                refend += 1
                # ('Hit', 'Ref RefStart RefEnd Query QueryStart QueryEnd Score')
                hit = Hit(read.Name, refstart, refend, primer_acc, 0,
                          len(primer_seq),  ed / len(primer_seq))
                all_locations.append(hit)
    refined_locations = refine_locations(read, all_primers, all_locations)
    return refined_locations
