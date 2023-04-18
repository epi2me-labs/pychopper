# -*- coding: utf-8 -*-

import sys
import os
import re
import subprocess as sp
import itertools
from collections import defaultdict
from pychopper.common_structures import Hit
from pychopper import utils


def _parse_hmmscan_tab(lines, reads):
    "Prase nhmmscan tabular output"
    res = []
    for line in lines:
        line = line.strip()
        if len(line) == 0:
            continue
        if line.startswith('#'):
            continue
        tmp = re.split(r'\s+', line)
        hit = Hit(tmp[2], int(tmp[6]), int(tmp[7]) + 1, tmp[0], int(tmp[4]), int(tmp[5]) + 1, float(tmp[12]))
        res.append(hit)

    buff = defaultdict(list)
    for x, y in itertools.groupby(res, lambda x: x.Ref):
        buff[x] = list(y)
    for r in reads:
        yield buff[r.Id.split()[0]]


def find_locations(reads, phmm_file, E, pool, min_batch):
    "Find alignment hits of all primers in all reads using the pHMM/nhmmscan backend"
    batches = list(utils.batch(reads, min_batch))
    for res in pool.map(_find_locations_single, zip(batches, [(phmm_file, E, 1)] * len(batches))):
        for h in res:
            yield list(h)


def _find_locations_single(params):
    "Find alignment hits of all primers in a single reads using the pHMM/nhmmscan backend"
    reads = params[0]
    phmm_file, E, threads = params[1]
    if not os.path.isfile(phmm_file):
        raise Exception("Profile HMM file is invalid: " + phmm_file)
    cmd = "nhmmscan --notextw --max -E {} --cpu {} --watson -o /dev/null --tblout /dev/stdout {} -"
    cmd = cmd.format(E, threads, phmm_file)
    with sp.Popen(cmd, shell=True, stdin=sp.PIPE, stdout=sp.PIPE) as proc:
        in_data = ""
        for read in reads:
            in_data += ">{}\n{}\n".format(read.Id, read.Seq)
        dout, derr = proc.communicate(in_data.encode())

    if derr is not None:
        print("Failed to run hmmscan with model {} and read {}:".format(phmm_file, read.Name), file=sys.stderr)
        sys.exit(1)

    return list(_parse_hmmscan_tab(dout.decode().split("\n"), reads))
