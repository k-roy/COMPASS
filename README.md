# Comparison Of Multiple alignment Programs for Alternative Splice Site discovery (COMPASS)
COMPASS identifies splice junctions in RNA-seq data with high precision and sensitivity. It is designed to detect alternative splicing events, especially those involving unannotated and non-canonical splice sites. In addition, COMPASS can be used to obtain an optimal set of alignments for reads in general.

## Why COMPASS?
COMPASS was inspired by the observation that different aligners have complementary strengths and weaknesses when mapping spliced reads across introns. In many cases, it is straightforward to determine which aligner performed the "best" and therefore more likely aligned a read or read pair correctly. To improve read mapping accuracy, we developed a systematic computational approach to select the best alignment for each read from a panel of aligners with diverse mapping strategies.

## What does COMPASS do?
COMPASS first passes RNA-seq reads through different aligners and then selects the best alignment based on fewest mismatches with the reference. Ties (alignments with the same score but differing on junction location) are broken by specific criteria, first favoring ungapped alignments over gapped ones, then annotated introns over unannotated ones, and finally shorter introns over longer ones. Under the hood, COMPASS handles numerous implementation details prior to scoring each alignment, including trimming and quality filtering of raw reads, disabling soft-clipping to ensure that aligners map the reads in their entirety, and reformatting alignment representations to ensure consistency for comparisons.

## 1) COMPASS_install_required_programs.sh
  * Creates the compass environment in the conda package manager.
  * Installs samfixcigar from jvarkit.

## 2) COMPASS.sh:
This is the core program. COMPASS.sh calls the following scripts.

### a) process_reads_and_align.sh:
  * Read 3′ trimming for polyA tails and base calling quality scores (cutadapt).
  * Assigns consecutive numbers to each read in the fastq file (awk).
  * Alignment with multiple aligners: BBMap, STAR (both default and noncanonical splicing modes), HISAT2 (both default and noncanonical splicing modes), Magic-BLAST, and GSNAP.
  * Prior to or as part of the alignment program calls, genome indices are built for each aligner only if necessary (i.e. when running the first time).
  * Sam files are converted to bam format (samtools view), sorted by mapped reference coordinates (samtools sort), cigars reformatted to SAM format 1.4 (samfixcigar), and then sorted by read numbers (samtools sort).

### b) compare_splice_junctions_from_multiple_aligners.py:
  * Comparison of Multiple Programs for Accurate Splice Site discovery (COMPASS) core program
  * Adjustment of ambiguous junctions to most likely splice sites based on species-specific splice signals

### c) analyze_exonic_and_intronic_sequence.py:
  * Extracts sequences surrounding the detected splice sites to examine motifs.
  * Selects the most likely branchpoint based on location and similarity to the consensus branchpoint sequence.

### d) create_splice_site_bed.py
  * Writes a bed file from the obtained splice sites for extracting total read depth from the unspliced bam files (samtools depth).

### e) add_unspliced_read_counts_to_junctions.py
  * Modifies the junction table created in (b) with the unspliced read counts for each junction from (d) for calculating splicing efficiency (SE).

## 3) COMPASS_combine_junction_tables_from_multiple_samples.py
  * Performs alternative splice junction calling and junction quality filtering.
