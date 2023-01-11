# Comparison of Multiple alignment Programs for Alternative Splice Site discovery (COMPASS)
COMPASS identifies splice junctions in RNA-seq data with high precision and sensitivity. It is specifically designed to detect alternative splicing events, especially those involving unannotated and non-canonical splice sites.

## Why COMPASS?
COMPASS was inspired by the observation that different aligners have complementary strengths and weaknesses with mapping junctions. In many cases, it is straightforward to determine which aligner performed the "best" and therefore more likely aligned a read or read pair correctly. We therefore decided to develop a systematic computational approach to integrate the best alignments from a set of aligners with diverse mapping strategies.

## What does COMPASS do?
At a high-level, COMPASS simply passes reads through different aligners and selects the best alignment based on fewest mismatches with the reference. Ties (alignments with the same score but differing on junction location) are broken by specific criteria. Under the hood, COMPASS handles many implementation details prior to scoring each alignment. This includes trimming and quality filtering of the raw reads, disabling soft-clipping to ensure that aligners map the reads in their entirety, and reformatting alignment representations to ensure consistency for comparisons.

## 1) COMPASS_install_required_programs.sh
  * Creates the compass environment in conda.
  * Installs samfixcigar from jvarkit.

## 2) COMPASS.sh (run for each sample):
This is the core program. COMPASS.sh calls the following scripts.

### a) process_reads_and_align.sh:
  * Read 3′ trimming for polyA tails and base calling quality scores (cutadapt).
  * Assigns consecutive numbers to each read in the fastq file (awk).
  * Alignment with multiple aligners: BBMap, STAR (both default and noncanonical splicing modes), HISAT2 (both default and noncanonical splicing modes), Magic-BLAST, and GSNAP.
  * Prior to or as part of the alignment program calls, genome indices are built for each aligner only if necessary (i.e. when running the first time).
  * Sam files are converted to .bam format (samtools view), sorted by mapped reference coordinates (samtools sort), cigars reformatted to SAM format 1.4 (samfixcigar), and then sorted by read numbers (samtools sort).

### b) compare_splice_junctions_from_multiple_aligners.py:
  * Comparison of Multiple Programs for Accurate Splice Site discovery (COMPASS) core program
  * Adjustment of ambiguous junctions to most likely splice sites based on species-specific splice signals

### c) add_exonic_intronic_sequence.py:
  * Adjustment of ambiguous junctions to most likely splice sites based on species-specific splice signals

### d) create_splice_site_bed.py
  * Writes a bed file from the obtained splice sites

### e) add_unspliced_read_counts_to_junctions.py

## 3) COMPASS_combine_junction_tables_from_multiple_samples.py
  * This script performs alternative splice junction calling and junction quality filtering.
