![ONT_logo](/ONT_logo.png)

-----------------------------

Pychopper
=========

Pychopper v2 is a tool to identify, orient and trim full-length Nanopore cDNA reads. The tool is also able to rescue fused reads.

Background
----------
The general approach of Pychopper v2 is the following:

- Pychopper first identifies alignment hits of the primers across the length of the sequence. The default method for doing this is using `nhmmscan` with the pre-trained strand specific profile HMMs, included with the package. Alternatively, one can use the `edlib` backend, which uses a combination of global and local alignment to identify the primers within the read.
- After identifying the primer hits by either of the backends, the reads are divided into segments defined by two consecutive primer hits. The score of a segment is its length if the configuration of the flanking primer hits is valid (such as `SPP,-VNP` for forward reads) or zero otherwise.
- The segments are assigned to rescued reads using a dynamic programming algorithm maximizing the sum of used segment scores (hence the amount of rescued bases). A crucial observation about the algorithm is that if a segment is included as a rescued read, then the next segment must be excluded as one of the primer hits defining it was "used up" by the previous segment. This put constraints on the dynamic programming graph, as illustrated in the figure below. The arrows in read define the optimal path for rescuing two fused reads with the a total score of `l1 + l3`.

![dp_segmentation](/dp_segmentation.png)

- A crucial parameter of Pychopper v2 is `-q`, which determines the stringency of primer alignment (E-value in the case of the pHMM backend). This can be explicitly specified by the user, however by default it is optimized on a random sample of input reads to produce the maximum number of classified reads.

Getting Started
================

## Installation

Install using conda :

```bash
conda install -c epi2melabs -c bioconda -c conda-forge “pychopper>=2.7.0”
```

## Usage

```
usage: pychopper [-h] [-b primers] [-g phmm_file] [-c config_file]
                          [-k kit] [-q cutoff] [-Q min_qual] [-z min_len]
                          [-r report_pdf] [-u unclass_output]
                          [-l len_fail_output] [-w rescue_output]
                          [-S stats_output] [-K qc_fail_output]
                          [-Y autotune_nr] [-L autotune_samples]
                          [-A scores_output] [-m method] [-x rescue] [-p]
                          [-t threads] [-B batch_size] [-D read stats]
                          input_fastx output_fastx

Tool to identify, orient and rescue full-length cDNA reads.

positional arguments:
  input_fastx          Input file.
  output_fastx         Output file.

optional arguments:
  -h, --help           show this help message and exit
  -b primers           Primers fasta.
  -g phmm_file         File with custom profile HMMs (None).
  -c config_file       File to specify primer configurations for each
                       direction (None).
  -k kit               Use primer sequences from this kit (PCS109).
  -q cutoff            Cutoff parameter (autotuned).
  -Q min_qual          Minimum mean base quality (7.0).
  -z min_len           Minimum segment length (50).
  -r report_pdf        Report PDF (pychopper_report.pdf).
  -u unclass_output    Write unclassified reads to this file.
  -l len_fail_output   Write fragments failing the length filter in this file.
  -w rescue_output     Write rescued reads to this file.
  -S stats_output      Write statistics to this file.
  -K qc_fail_output    Write reads failing mean quality filter to this file.
  -Y autotune_nr       Approximate number of reads used for tuning the cutoff
                       parameter (10000).
  -L autotune_samples  Number of samples taken when tuning cutoff parameter
                       (30).
  -A scores_output     Write alignment scores to this BED file.
  -m method            Detection method: phmm or edlib (phmm).
  -x rescue            Protocol-specific read rescue: DCS109 (None).
  -p                   Keep primers, but trim the rest.
  -t threads           Number of threads to use (8).
  -B batch_size        Maximum number of reads processed in each batch
                       (1000000).
```

*WARNING: Do not turn on trimming during basecalling as it will remove the primers needed for classifying the reads!*

### Basic usage

Example usage with default PCS109/DCS109 primers using the default pHMM backend:

```bash
pychopper -r report.pdf -u unclassified.fq -w rescued.fq input.fq full_length_output.fq
```

Example usage with default PCS109/DCS109 primers using the edlib/parasail backend:

```bash
pychopper -m edlib -r report.pdf -u unclassified.fq -w rescued.fq input.fq full_length_output.fq
```
Example usage with default PCS109/DCS109 primers using the default pHMM backend:

```bash
pychopper -r report.pdf -A aln_hits.bed -S statistics.tsv -u unclassified.fq -w rescued.fq input.fq full_length_output.fq
```

### Advanced usage

The fasta files with custom primers used by the `edlib/parasail` backend can be specified through `-b`, while the valid primer configurations are specified through `-c`:

```bash
pychopper -m edlib -b custom_pimers.fas -c primer_config.txt input.fq full_length_output.fq
```
Where the contents of `primer_config.txt` looks like `+:MySSP,-MyVNP|-:MyVNP,-MySSP`.

The `pHMM` alignment backend takes a "compressed" profile HMM trained from a multiple sequence alignment using the [hmmer](http://hmmer.org/) package. Custom profile HMMs can be trained from a fastq of reads and a fasta file with the primer sequences using the [hammerpede](https://github.com/nanoporetech/hammerpede) package. The path to the custom profile HMM can be specified using `-g`:

```bash
pychopper -m phmm -g MySSP_MyVNP.hmm -c primer_config.txt input.fq full_length_output.fq
```

Help
====

## Licence and Copyright

(c) 2019 Oxford Nanopore Technologies Ltd.

This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this
file, You can obtain one at http://mozilla.org/MPL/2.0/.

## FAQs and tips

## References and Supporting Information

See the post announcing the tools at the Oxford Nanopore Technologies Community [here](https://community.nanoporetech.com/posts/new-transcriptomics-analys).

## Research Release

Research releases are provided as technology demonstrators to provide early access to features or stimulate Community development of tools. Support for this software will be minimal and is only provided directly by the developers. Feature requests, improvements, and discussions are welcome and can be implemented by forking and pull requests. However much as we would like to rectify every issue and piece of feedback users may have, the developers may have limited resource for support of this software. Research releases may be unstable and subject to rapid iteration by Oxford Nanopore Technologies.


