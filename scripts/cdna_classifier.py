#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import os
import io
import sys
import numpy as np
import pandas as pd
from collections import OrderedDict, defaultdict
import concurrent.futures
import tqdm
import gzip

from pychopper import seq_utils as seu
from pychopper import utils
from pychopper import chopper, report
import pychopper.phmm_data as phmm_data
import pychopper.primer_data as primer_data


"""
Parse command line arguments.
"""
parser = argparse.ArgumentParser(
    description='Tool to identify, orient and rescue full-length cDNA reads.')
parser.add_argument(
    '-b', metavar='primers', type=str, default=None, help="Primers fasta.", required=False)
parser.add_argument(
    '-g', metavar='phmm_file', type=str, default=None, help="File with custom profile HMMs (None).", required=False)
parser.add_argument(
    '-c', metavar='config_file', type=str, default=None, help="File to specify primer configurations for each direction (None).")
parser.add_argument(
    '-k', metavar='kit', type=str, default="PCS109", help="Use primer sequences from this kit (PCS109).")
parser.add_argument(
    '-q', metavar='cutoff', type=str, default=None, help="Cutoff parameter (autotuned).")
parser.add_argument(
    '-Q', metavar='min_qual', type=float, default=7.0, help="Minimum mean base quality (7.0).")
parser.add_argument(
    '-z', metavar='min_len', type=int, default=50, help="Minimum segment length (50).")
parser.add_argument(
    '-r', metavar='report_pdf', type=str, default="cdna_classifier_report.pdf", help="Report PDF (cdna_classifier_report.pdf).")
parser.add_argument(
    '-u', metavar='unclass_output', type=str, default=None, help="Write unclassified reads to this file.")
parser.add_argument(
    '-l', metavar='len_fail_output', type=str, default=None, help="Write fragments failing the length filter in this file.")
parser.add_argument(
    '-w', metavar='rescue_output', type=str, default=None, help="Write rescued reads to this file.")
parser.add_argument(
    '-S', metavar='stats_output', type=str, default="cdna_classifier_report.tsv", help="Write statistics to this file.")
parser.add_argument(
    '-K', metavar='qc_fail_output', type=str, default=None, help="Write reads failing mean quality filter to this file.")
parser.add_argument(
    '-Y', metavar='autotune_nr', type=float, default=10000, help="Approximate number of reads used for tuning the cutoff parameter (10000).")
parser.add_argument(
    '-L', metavar='autotune_samples', type=int, default=30, help="Number of samples taken when tuning cutoff parameter (30).")
parser.add_argument(
    '-A', metavar='scores_output', type=str, default=None, help="Write alignment scores to this BED file.")
parser.add_argument(
    '-m', metavar='method', type=str, default="phmm", help="Detection method: phmm or edlib (phmm).")
parser.add_argument(
    '-x', metavar='rescue', type=str, default=None, help="Protocol-specific read rescue: DCS109 (None).")
parser.add_argument(
    '-p', action='store_true', default=False, help="Keep primers, but trim the rest.")
parser.add_argument(
    '-t', metavar='threads', type=int, default=8, help="Number of threads to use (8).")
parser.add_argument(
    '-B', metavar='batch_size', type=int, default=1000, help="Maximum number of reads processed in each batch (1000000).")
parser.add_argument(
    '-D', metavar='read stats', type=str, default=None, help="Tab separated file with per-read stats (None).")
parser.add_argument(
    '-y', action='store_true', default=False, help="Output FASTQ comment as BAM tags. Use with minimap2 -y to pass UMI and additional info into BAM file.")
parser.add_argument(
    '-U', action='store_true', default=False, help="Detect UMIs")
parser.add_argument('input_fastx', metavar='input_fastx', type=str, help="Input file.")
parser.add_argument('output_fastx', metavar='output_fastx', nargs="?", type=str, default="-", help="Output file.")


def _new_stats():
    "Initialize a new statistic dictionary"
    st = OrderedDict()
    st["Classification"] = OrderedDict([('Primers_found', 0), ('Rescue', 0), ('Unusable', 0)])
    st["Strand"] = OrderedDict([('+', 0), ('-', 0)])
    st["RescueStrand"] = OrderedDict([('+', 0), ('-', 0)])
    st["RescueSegmentNr"] = defaultdict(int)
    st["RescueHitNr"] = defaultdict(int)
    st["UnclassHitNr"] = defaultdict(int)
    st["Unusable"] = defaultdict(int)
    st["Hits"] = defaultdict(int)
    st["LenFail"] = 0
    st["QcFail"] = 0
    st["PassReads"] = 0
    return st


def _update_stats(st, d_fh, segments, hits, usable_len, read):
    "Update stats dictionary with properties of a read"
    st["PassReads"] += 1
    if len(hits) > 0:
        h = ",".join([x.Query for x in hits])
        st["Hits"][h] += 1
    if len(segments) == 0:
        st["Classification"]["Unusable"] += 1
        st["UnclassHitNr"][len(hits)] += 1
        if d_fh is not None:
            d_fh.write("{}\t{}\t0\t{}\t{}\t{}\n".format(read.Name, len(read.Seq), -1, -1, "."))
    elif len(segments) == 1:
        st["Classification"]["Primers_found"] += 1
        st["Strand"][segments[0].Strand] += 1
        st["Unusable"][int(segments[0].Len / len(read.Seq) * 100)] += 1
        if d_fh is not None:
            rs = segments[0]
            d_fh.write("{}\t{}\t1\t{}\t{}\t{}\n".format(read.Name, len(read.Seq), rs.Start, rs.End, rs.Strand))
    else:
        for rs in segments:
            st["Classification"]["Rescue"] += 1
            st["RescueStrand"][rs.Strand] += 1
            st["RescueHitNr"][len(hits)] += 1
            if d_fh is not None:
                d_fh.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(read.Name, len(read.Seq), len(segments), rs.Start, rs.End, rs.Strand))
        st["Unusable"][len(read.Seq) - int(sum([s.Len for s in segments]))] += 1
        st["RescueSegmentNr"][len(segments)] += 1


def _process_stats(st):
    "Convert stats dictionary into a data frame"
    res = OrderedDict([("Category", []), ("Name", []), ("Value", [])])
    res["Category"] += ["ReadStats"]
    res["Name"] += ["PassReads"]
    res["Value"] += [st["PassReads"]]
    res["Category"] += ["ReadStats"]
    res["Name"] += ["LenFail"]
    res["Value"] += [st["LenFail"]]
    res["Category"] += ["ReadStats"]
    res["Name"] += ["QcFail"]
    res["Value"] += [st["QcFail"]]
    for k, v in st['Classification'].items():
        res["Category"] += ["Classification"]
        res["Name"] += [k]
        res["Value"] += [v]
    for k, v in st['Strand'].items():
        res["Category"] += ["Strand"]
        res["Name"] += [k]
        res["Value"] += [v]
    for k, v in st['RescueStrand'].items():
        res["Category"] += ["RescueStrand"]
        res["Name"] += [k]
        res["Value"] += [v]
    for k, v in sorted(st['UnclassHitNr'].items(), key=lambda x: x[0]):
        res["Category"] += ["UnclassHitNr"]
        res["Name"] += [k]
        res["Value"] += [v]
    for k, v in sorted(st['RescueHitNr'].items(), key=lambda x: x[0]):
        res["Category"] += ["RescueHitNr"]
        res["Name"] += [k]
        res["Value"] += [v]
    for k, v in sorted(st['RescueSegmentNr'].items(), key=lambda x: x[0]):
        res["Category"] += ["RescueSegmentNr"]
        res["Name"] += [k]
        res["Value"] += [v]
    for k, v in sorted(st['Unusable'].items(), key=lambda x: x[0]):
        res["Category"] += ["Unusable"]
        res["Name"] += [k]
        res["Value"] += [v]
    for k, v in sorted(st['Hits'].items(), key=lambda x: x[1], reverse=True):
        res["Category"] += ["Hits"]
        res["Name"] += [k]
        res["Value"] += [v]
    res = pd.DataFrame(res)
    return res


def _detect_anomalies(st, config):
    raw_anom = []
    total = st["PassReads"]
    for k, v in sorted(st['Hits'].items(), key=lambda x: x[1], reverse=True):
        x = tuple(k.split(","))
        if len(x) != 2:
            raw_anom.append((k, v))
        elif x not in config:
            raw_anom.append((k, v))
    anom = []
    for tmp in raw_anom:
        pc = tmp[1] * 100.0 / total
        if pc >= 3.0:
            anom.append((tmp[0], tmp[1], pc))
    if len(anom) > 0:
        sys.stderr.write("Detected {} potential artefactual primer configurations:\n".format(len(anom)))
        ml = max([len(t[0]) for t in anom])
        mn = max([len(str(t[1])) for t in anom])
        sys.stderr.write(("{:" + str(ml) + "}\t{:" + str(mn) + "}\t{}\n").format("Configuration", "NrReads", "PercentReads"))
        for x in anom:
            fs = "{:" + str(ml) + "}\t{:" + str(mn) + "}\t{:.2f}%\n"
            sys.stderr.write(fs.format(x[0], str(x[1]), x[2]))


def _plot_pd_bars(df, title, report, alpha=0.7, xrot=0, ann=False):
    "Generate bar plot from data frame"
    df = df.drop("Category", axis=1)
    df = df.set_index("Name")
    report.plt.clf()
    if sum(df.Value) > 0:
        ax = df.plot(kind='bar', align='center', alpha=alpha, rot=xrot)
        report.plt.title(title)
        if ann:
            y_offset = 1
            for rect in ax.patches:
                height = rect.get_height()
                ax.text(rect.get_x() + rect.get_width() / 2., height + y_offset, "{:.0f}".format(height), ha='center', va='bottom', rotation=0)
                ax.text(rect.get_x() + rect.get_width() / 2., 0.5 * height - 2 * y_offset, "{:.1f}%".format(100 * height / sum(df.Value)), ha='center', va='bottom', rotation=0)
        report.save_close()


def _plot_pd_line(df, title, report, alpha=0.7, xrot=0, vline=None):
    "Generate line plot from a data frame"
    df = df.drop("Category", axis=1)
    df = df.set_index("Name")
    report.plt.clf()
    if sum(df.Value) > 0:
        df.plot(kind='line', rot=xrot)
        if vline is not None:
            report.plt.axvline(x=vline, color="red")
        report.plt.title(title)
        report.save_close()


def _plot_stats(st, pdf):
    "Generate plots and save to report PDF"
    R = report.Report(pdf)
    rs = st.loc[st.Category == "Classification", ]
    _plot_pd_bars(rs.copy(), "Classification of output reads", R, ann=True)
    found, rescue, unusable = float(rs.loc[rs.Name == "Primers_found", ].Value), float(rs.loc[rs.Name == "Rescue", ].Value), float(rs.loc[rs.Name == "Unusable", ].Value)
    total = found + rescue + unusable
    found = found / total * 100
    rescue = rescue / total * 100
    unusable = unusable / total * 100
    sys.stderr.write("-----------------------------------\n")
    sys.stderr.write("Reads with two primers:\t{:.2f}%\nRescued reads:\t\t{:.2f}%\nUnusable reads:\t\t{:.2f}%\n".format(found, rescue, unusable))
    sys.stderr.write("-----------------------------------\n")
    _plot_pd_bars(st.loc[st.Category == "Strand", ].copy(), "Strand of oriented reads", R, ann=True)
    _plot_pd_bars(st.loc[st.Category == "RescueStrand", ].copy(), "Strand of rescued reads", R, ann=True)
    _plot_pd_bars(st.loc[st.Category == "UnclassHitNr", ].copy(), "Number of hits in unclassified reads", R)
    _plot_pd_bars(st.loc[st.Category == "RescueHitNr", ].copy(), "Number of hits in rescued reads", R)
    _plot_pd_bars(st.loc[st.Category == "RescueSegmentNr", ].copy(), "Number of usable segments per rescued read", R)
    if q_bak is None:
        _plot_pd_line(st.loc[st.Category == "AutotuneSample", ].copy(), "Usable bases as function of cutoff(q). Best q={:.4g}".format(args.q), R, vline=args.q)
    udf = st.loc[st.Category == "Unusable", ].copy()
    udf.Name = np.log10(1.0 + np.array(udf.Name, dtype=float))
    _plot_pd_line(udf, "Log10 length distribution of trimmed away sequences.", R)
    R.close()


def _opener(filename, mode, encoding='utf8'):
    if filename == '-':
        sys.stderr.write("Reading from stdin\n")
        #return open(sys.stdin.buff, mode, encoding=encoding)
        return io.TextIOWrapper(sys.stdin.buffer, encoding=encoding)
    elif filename.endswith('.gz'):
        return gzip.open(filename, mode + 't', encoding=encoding)
    else:
        return open(filename, mode, encoding=encoding)

if __name__ == '__main__':
    args = parser.parse_args()

    if args.m == "phmm":
        utils.check_command("nhmmscan -h > /dev/null")
        utils.check_min_hmmer_version(3, 2)

    CONFIG = "+:SSP,-VNP|-:VNP,-SSP"
    if args.c is not None:
        CONFIG = open(args.c, "r").readline().strip()

    kits = {"PCS109": {"HMM": os.path.join(os.path.dirname(phmm_data.__file__), "cDNA_SSP_VNP.hmm"), "FAS": os.path.join(os.path.dirname(primer_data.__file__), "cDNA_SSP_VNP.fas"), }, "PCS110": {
                                           "HMM": os.path.join(os.path.dirname(phmm_data.__file__), "PCS110_primers.hmm"), "FAS": os.path.join(os.path.dirname(primer_data.__file__), "PCS110_primers.fas")}}

    if args.g is None:
        args.g = kits[args.k]["HMM"]

    if args.b is None:
        args.b = kits[args.k]["FAS"]

    if args.x is not None and args.x in ('DCS109'):
        if args.x == "DCS109":
            CONFIG = "-:VNP,-VNP"

    config = utils.parse_config_string(CONFIG)
    sys.stderr.write("Using kit: {}\n".format(args.k))
    sys.stderr.write("Configurations to consider: \"{}\"\n".format(CONFIG))

    #in_fh = sys.stdin
    #if args.input_fastx != '-':
    in_fh = _opener(args.input_fastx, "r")

    out_fh = sys.stdout
    if args.output_fastx != '-':
        out_fh = open(args.output_fastx, "w")

    u_fh = None
    if args.u is not None:
        u_fh = open(args.u, "w")

    l_fh = None
    if args.l is not None:
        l_fh = open(args.l, "w")

    w_fh = None
    if args.w is not None:
        w_fh = open(args.w, "w")

    a_fh = None
    if args.A is not None:
        a_fh = open(args.A, "w")

    d_fh = None
    if args.D is not None:
        d_fh = open(args.D, "w")
        d_fh.write("Read\tLength\tStatus\tStart\tEnd\tStrand\n")

    st = _new_stats()
    #input_size = None
    #if args.input_fastx != "-":
    #    input_size = os.stat(args.input_fastx).st_size

    if args.q is None and args.Y <= 0:
        sys.stderr.write("Please specifiy either -q or -Y!")

    if args.m == "edlib":
        all_primers = seu.get_primers(args.b)

    if args.m == "phmm":
        if args.q is None:
            args.q = "0.00001"
        def backend(x, pool, q=None, mb=None):
            return chopper.chopper_phmm(x, args.g, config, q, args.t, pool, mb)
    elif args.m == "edlib":
        if args.q is None:
            args.q = "0.5"
        def backend(x, pool, q=None, mb=None):
            return chopper.chopper_edlib(x, all_primers, config, q * 1.2, q, pool, mb)
    else:
        raise Exception("Invalid backend!")

    # Pick the -q maximizing the number of classified reads using grid search:
    #nr_records = None
    tune_df = None
    q_bak = args.q
    if args.q == 'auto':
        if args.input_fastx == '-':
            sys.stderr.write("Can't auto tune -q when reading from stdin.\n")
            sys.exit(1)
        nr_cutoffs = args.L
        cutoffs = np.linspace(0.0, 1.0, num=nr_cutoffs)
        cutoffs = cutoffs / cutoffs[-1]
        if args.m == "phmm":
            cutoffs = np.linspace(10**-5, 5.0, num=nr_cutoffs)
        class_reads = []
        class_readLens = []
        target_prop = 1.0
        sys.stderr.write("Counting fastq records in input file: {}\n".format(args.input_fastx))
        #sys.stderr.write("Total fastq records in input file: {}\n".format(nr_records))
        read_sample = list(seu.readfq(_opener(args.input_fastx, "r"), sample=target_prop, min_qual=args.Q))
        sys.stderr.write("Tuning the cutoff parameter (q) on {} sampled reads ({:.1f}%) passing quality filters (Q >= {}).\n".format(len(read_sample), target_prop * 100.0, args.Q))
        sys.stderr.write("Optimizing over {} cutoff values.\n".format(args.L))
        for qv in tqdm.tqdm(cutoffs):
            clsLen = 0
            cls = 0
            with concurrent.futures.ProcessPoolExecutor(max_workers=args.t) as executor:
                for batch in utils.batch(read_sample, int((len(read_sample)))):
                    for read, (segments, hits, usable_len) in backend(batch, executor, qv, max(1000, int((len(read_sample)) / args.t))):
                        flt = list([x.Len for x in segments if x.Len > 0])
                        if len(flt) == 1:
                            clsLen += sum(flt)
                            cls += 1
            class_reads.append(cls)
            class_readLens.append(clsLen)
        best_qi = np.argmax(class_readLens)
        args.q = cutoffs[best_qi]
        tune_df = OrderedDict([("Category", []), ("Name", []), ("Value", [])])
        for i, c in enumerate(cutoffs):
            tune_df["Category"] += ["AutotuneSample"]
            tune_df["Name"] += [c]
            tune_df["Value"] += [class_reads[i]]

        tune_df["Category"] += ["Parameter"]
        tune_df["Name"] += ["Cutoff(q)"]
        tune_df["Value"] += [args.q]
        if best_qi == (len(class_reads) - 1):
            sys.stderr.write("Best cuttoff value is at the edge of the search interval! Using tuned value is not safe! Please pick a q value manually and QC your data!\n")
        sys.stderr.write("Best cutoff (q) value is {:.4g} with {:.0f}% of the reads classified.\n".format(args.q, class_reads[best_qi] * 100 / len(read_sample)))
    else:
        args.q = float(args.q)
    #if nr_records is not None:
    #    input_size = nr_records
    #    if args.B > nr_records:
    #        args.B = nr_records
    #    if args.B == 0:
    #        args.B = 1

    sys.stderr.write("Processing the whole dataset using a batch size of {}:\n".format(args.B))
    pbar = tqdm.tqdm(unit=" reads")
    min_batch_size = max(int(args.B / args.t), 1)
    rfq_sup = {"out_fq": args.K, "pass": 0, "total": 0}
    with concurrent.futures.ProcessPoolExecutor(max_workers=args.t) as executor:
        for batch in utils.batch(seu.readfq(in_fh, min_qual=args.Q, rfq_sup=rfq_sup), args.B):
            for read, (segments, hits, usable_len) in backend(batch, executor, q=args.q, mb=min_batch_size):
                if args.A is not None:
                    for h in hits:
                        a_fh.write(utils.hit2bed(h, read) + "\n")
                _update_stats(st, d_fh, segments, hits, usable_len, read)
                if args.u is not None and len(segments) == 0:
                    seu.writefq(read, u_fh)
                for trim_read in chopper.segments_to_reads(read, segments, args.p, args.y, args.U):
                    if args.l is not None and len(trim_read.Seq) < args.z:
                        st["LenFail"] += 1
                        seu.writefq(trim_read, l_fh)
                        continue
                    if len(segments) == 1:
                        seu.writefq(trim_read, out_fh)
                    if args.w is not None and len(segments) > 1:
                        seu.writefq(trim_read, w_fh)
                #if nr_records is None:
                #    pbar.update(seu.record_size(read, 'fastq'))
                #else:
                pbar.update(1)
    pbar.close()
    sys.stderr.write("Finished processing file: {}\n".format(args.input_fastx))
    fail_nr = rfq_sup["total"] - rfq_sup["pass"]
    fail_pc = (fail_nr * 100 / rfq_sup["total"])
    st["QcFail"] = fail_nr
    sys.stderr.write("Input reads failing mean quality filter (Q < {}): {} ({:.2f}%)\n".format(args.Q, fail_nr, fail_pc))
    sys.stderr.write("Output fragments failing length filter (length < {}): {}\n".format(args.z, st["LenFail"]))

    # Save stats as TSV:
    stdf = None
    if args.S is not None or args.r is not None:
        stdf = _process_stats(st)
        if tune_df is not None:
            stdf = stdf.append(pd.DataFrame(tune_df))

    _detect_anomalies(st, config)

    if args.S is not None:
        stdf.to_csv(args.S, sep="\t", index=False)

    for fh in (in_fh, out_fh, u_fh, l_fh, w_fh, a_fh, d_fh):
        if fh is None:
            continue
        fh.flush()
        fh.close()

    if args.r is not None:
        _plot_stats(stdf, args.r)
