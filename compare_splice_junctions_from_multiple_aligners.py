# -*- coding: utf-8 -*-
'''
@author: kevinroy

This script selects the optimal alignment based on a score calculated from user-defined penalties for each CIGAR operation. 
-Ties are broken in various scenarios: 
    -ungapped alignment favored over gapped alignment (to be conservative)
    -annotated junctions favored over unannotated junctions (to be conservative)
    -smallest intron without junction-proximal mismatches favored (for parsimony)

output files:

-optimal COMPASS alignment as a single BAM file, along with separate bam files for ungapped (UGA), 
annotated spliced junctions (ASJ), and unannotated spliced junctions (USJ) files.

-alignment info at the level of each read pair * aligner, unannotated_gapped_alignments_outfile
    columns:
    'read_pair_ID', 'aligner_has_min_score', 'min_score_has_annotated_intron', 'min_score_has_no_intron', \
    'aligner', 'score',  'num_introns_found', 'perfect_gapped_alignment', 'chrom', 'amb_start', \
    'amb_stop', 'five_SS', 'three_SS', 'RNA_strand', 'annotated_junction', 'ann_5SS', 'ann_3SS', \
    'canonical_5SS', 'canonical_3SS', 'intron_size', 'intron_coords_adjusted', 'mismatch_near_junction', \
    'US_perfect_matches', 'DS_perfect_matches', 'read1_flag', 'read1_chrom', 'read1_coord', \
    'read1_cigar', 'read1_NH', 'read2_flag', 'read2_chrom', 'read2_coord', 'read2_cigar', 'read2_NH'

-alignment info at the level of each junction identified by COMPASS, an 'integrated' junction file with \
    summary statistics on agreement between aligners for each junction, including the read counts for each aligner agree on the selected junction

    columns:
    sample_info_header = ['sample_name', 'total_reads']

    junction_header = 'chrom, start, stop, five_SS, three_SS, RNA_strand, num_amb_junctions, annotated_junction, \
        ann_5SS, ann_3SS, canonical_5SS, canonical_3SS, intron_size'.split(', ')
    counts_header = ['COMPASS_counts'] + [e + '_counts' for e in aligners]

    stats_header_part_1 = 'min_score_unannotated_intron_junctions_disagree, intron_coords_adjusted_counts, \
        mismatch_near_junction_counts, perfect_gapped_alignment_counts'.split(', ')

    stats_header_part_2 = 'US_perfect_matches, DS_perfect_matches, five_SS_Q_scores, three_SS_Q_scores, \
        alignment_scores, num_introns_found, read_1_coord_cigars, read_2_coord_cigars'.split(', ')

    stats_header_part_3 = 'num_unique_read_pair_cigars, max_US_perfect_matches, max_DS_perfect_matches, \
        unique_US_perfect_matches, unique_DS_perfect_matches, best_alignment_score, mean_alignment_score, median_alignment_score'.split(', ')

    stats_header_part_4 = 'mean_five_SS_Q_score, mean_three_SS_Q_score'.split(', ')

Notes:
-BAM for all aligners must be name-sorted, where read names are chronologically numbered prior to mapping

-CIGAR format must be N for intron gap, M for match, S for soft-clip, I for insertion, D for deletion, 
    and X for mismatch (version 1.4 and above, see: https://samtools.github.io/hts-specs/SAMv1.pdf)
    For some aligners, this requires fixing the CIGAR strings for each aligner to be consistent.
    COMPASS uses samfixcigar, which requires that alignments are first coord sorted. 
    All of this is handled in the script: COMPASS_process_reads_and_align.sh

-The UGA bam file is needed for tallying the unspliced reads at each COMPASS junction. The ASJ and USJ are
useful for examining the alignments of specific junctions. It may be desirable to isolate all the reads mapping to specific
junctions. This can be done with the script: extract_alignments_for_specific_COMPASS_junctions_from_multiple_aligners.py

-The read pair * aligner ASJ vs. USJ files are useful for exploring the technical \
    performance of each aligner for different types of junctions. The annotated file is typically very large.

-The COMPASS junction file is useful for exploring the biology of the splice junctions.
'''

import os
import pandas as pd
import pysam
import statistics
import timeit
import random 
import COMPASS_functions_human
import importlib
importlib.reload(COMPASS_functions)
from COMPASS_functions import *

start_time = timeit.default_timer()
random.seed(1234)
# When there are ties between aligners, one alignment must be randomly chosen. 
# Each aligner still gets credit for having the best score.

from sys import argv
COMPASS_DIR = argv[1]
GENOME_REF_DIR = argv[2]
sample_name = argv[3]
FASTA = argv[4]
INTRONS_FILE = argv[5]
NUM_THREADS = int(argv[6])
MIN_INTRON_LENGTH = int(argv[7])
MAX_INTRON_LENGTH = int(argv[8])
ALIGNERS_FILE = argv[9]
reads_to_process = int(argv[10])

if reads_to_process < 0:
    reads_to_process = 10**10

READ_NUM_PROGRESS = 10**3

LOG_DIR = COMPASS_DIR + 'log/'
SEPARATE_ALIGNMENTS_DIR = COMPASS_DIR + 'separate_alignments/'
COMPASS_ALIGNMENTS_DIR = COMPASS_DIR + 'COMPASS_alignments/'
JUNCTION_ALIGNMENTS_DIR = COMPASS_DIR + 'junction_read_alignments/'
COMPASS_JUNCTIONS_DIR = COMPASS_DIR + 'COMPASS_junctions/'

sample_suffix = '_name_sorted.bam'

WRITE_ALL_COMPASS_BAM = True
WRITE_UNANNOTATED_JUNCTIONS_BAM = True

MISMATCH_DIST_FROM_JUNCTION_DISALLOWED = 10 
# Mismatches near the junction have a chance of impacting alignment accuracy.
# A single mismatch within 10 bp of a junction raises a flag. These are reported and can be filtered later.

# Splice site penalties are only used to resolve ambiguous junction alignments. 
# These do not impact alignment scores and are not used for selecting the optimal aligner for a given read.

CONSENSUS_5SS = ['G', 'T', 'AG', 'TAC', 'G', 'T'] # the 6th position in human 5'SS does not seem to matter as much
PENALTIES_5SS = [6,   4,   1,   1,    2,  0.5]

CONSENSUS_3SS = ['CT', 'A', 'G']
PENALTIES_3SS = [1,     3,   4]

NUM_ITEMS_TO_REPORT = 10

genome_fasta = pysam.FastaFile(FASTA)

aligner_to_current_read_num = {}
aligners_to_bam = {}

aligner_df = pd.read_csv(ALIGNERS_FILE, sep = '\t')
aligner_df.columns
aligners = []
for index, row in aligner_df.iterrows():
    aligner = row['aligner']
    aligners.append(aligner)
    bam_filename = SEPARATE_ALIGNMENTS_DIR + aligner + '/' + sample_name + sample_suffix
    aligners_to_bam[aligner] = pysam.AlignmentFile(bam_filename, 'rb')
    aligner_to_current_read_num[aligner] = 0

# This file looks like this:
# NC_000001.11	12226	12612	+
# NC_000001.11	12720	13220	+

annotated_intron_df, ambiguous_junction_to_annotated_junction, ambiguous_annotated_junctions, annotated_introns, junction_to_intron_type = get_ambiguous_junctions_in_annotated_introns(INTRONS_FILE, genome_fasta)
# annotated_intron_df columns: chrom	start	end	strand	type	Name	intron_type, adjusted_start, adjusted_stop
# annotated_intron = (chrom, adjusted_start, adjusted_stop)
# ambiguous_junction_to_annotated_junction[ambiguous_junction] = annotated_intron
# ambiguous_annotated_junctions.add(annotated_intron)
# junction_to_intron_type[ambiguous_junction]  = row['intron_type'] e.g. spliceosomal_intron, tRNA_intron, or mitochondrial_mRNA_intron

###################### OPEN OUTPUT FILES FOR WRITING ##########################

# open files: COMPASS_bam_outfile, COMPASS_UGA_plus_strand_bam_outfile, COMPASS_UGA_minus_strand_bam_outfile, COMPASS_ASJ_bam_outfile, COMPASS_USJ_bam_outfile, annotated_gapped_alignments_outfile, unannotated_gapped_alignments_outfile, splice_junction_counts_outfile

# all reads written to COMPASS bam file
# write UGA (ungapped alignments), ASJ (annotated splice junction alignments) and USJ (unannotated splice junction alignments) separately
# UGA is written to separate strands for extracting unspliced coverage at intron junctions in a later step
# pysam pileup lacks the ability to distinguish between strands, so it is easier to write UGA to separate strands
# USJ enables viewing alignments for unannotated introns in IGV

# TODO: make these files optional to write, to enable saving space
# WRITE_ALL_COMPASS_BAM = True
# WRITE_UNANNOTATED_JUNCTIONS = True
     
log_outfilename = LOG_DIR + sample_name + '_compare_splice_junctions_from_multiple_aligners_python_output.log'
log_outfile = open(log_outfilename, 'w')

log_message = ('Starting COMPASS script compare_splice_junctions_from_multiple_aligners...\n')
print(log_message)
log_outfile.write(log_message)

COMPASS_bamfilename = COMPASS_ALIGNMENTS_DIR + sample_name + '_COMPASS.bam'
COMPASS_bam_outfile = pysam.AlignmentFile(COMPASS_bamfilename, 'wb', template = aligners_to_bam[aligners[0]])

COMPASS_UGA_plus_strand_bam_filename = COMPASS_ALIGNMENTS_DIR + sample_name + '_COMPASS_UGA_plus_strand.bam'
COMPASS_UGA_plus_strand_bam_outfile = pysam.AlignmentFile(COMPASS_UGA_plus_strand_bam_filename, 'wb', template = aligners_to_bam[aligners[0]])
COMPASS_UGA_minus_strand_bam_filename = COMPASS_ALIGNMENTS_DIR + sample_name + '_COMPASS_UGA_minus_strand.bam'
COMPASS_UGA_minus_strand_bam_outfile = pysam.AlignmentFile(COMPASS_UGA_minus_strand_bam_filename, 'wb', template = aligners_to_bam[aligners[0]])

COMPASS_ASJ_bam_filename = COMPASS_ALIGNMENTS_DIR + sample_name + '_COMPASS_ASJ.bam'
COMPASS_ASJ_bam_outfile = pysam.AlignmentFile(COMPASS_ASJ_bam_filename, 'wb', template = aligners_to_bam[aligners[0]])

COMPASS_USJ_bam_filename = COMPASS_ALIGNMENTS_DIR + sample_name + '_COMPASS_USJ.bam'
COMPASS_USJ_bam_outfile = pysam.AlignmentFile(COMPASS_USJ_bam_filename, 'wb', template = aligners_to_bam[aligners[0]])

COMPASS_sorted_bam_filename = COMPASS_ALIGNMENTS_DIR + sample_name + '_COMPASS_sorted.bam'
COMPASS_UGA_plus_strand_sorted_bam_filename = COMPASS_ALIGNMENTS_DIR + sample_name + '_COMPASS_UGA_plus_strand_sorted.bam'
COMPASS_UGA_minus_strand_sorted_bam_filename = COMPASS_ALIGNMENTS_DIR + sample_name + '_COMPASS_UGA_minus_strand_sorted.bam'
COMPASS_USJ_sorted_bam_filename = COMPASS_ALIGNMENTS_DIR + sample_name + '_COMPASS_USJ_sorted.bam'
COMPASS_ASJ_sorted_bam_filename = COMPASS_ALIGNMENTS_DIR + sample_name + '_COMPASS_ASJ_sorted.bam'

# read pair X aligner outfiles
annotated_gapped_alignments_outfilename = JUNCTION_ALIGNMENTS_DIR + sample_name + '_annotated_gapped.tsv'
annotated_gapped_alignments_outfile = open(annotated_gapped_alignments_outfilename, 'w')
# unnanoted junctions are written separately from annotated junctions, 
# as the annotated junction files can be too large to be loaded into memory in a data table
unannotated_gapped_alignments_outfilename = JUNCTION_ALIGNMENTS_DIR + sample_name + '_unannotated_gapped.tsv'
unannotated_gapped_alignments_outfile = open(unannotated_gapped_alignments_outfilename, 'w')   
# if an intron is found, info will be written to either annotated or unannotated outfiles, 
# depending on whether the alignment with the best score maps to an annotated or unannotated intron, respectively
intron_found_to_alignments_outfiles = {True: annotated_gapped_alignments_outfile, False: unannotated_gapped_alignments_outfile}

# This file is written for troubleshooting purposes,
# in order to see why aligners disagree when each has the best score.
junction_min_score_disagreement_outfilename = JUNCTION_ALIGNMENTS_DIR + sample_name + '_min_score_alignments_disagree_on_junction.tsv'
junction_min_score_disagreement_outfile = open(junction_min_score_disagreement_outfilename, 'w')   

# COMPASS junction outfile
splice_junction_counts_outfilename = COMPASS_JUNCTIONS_DIR + sample_name + '_COMPASS_splice_junctions.tsv'
splice_junction_counts_outfile = open(splice_junction_counts_outfilename, 'w')

read_num = 0
aligner_to_current_aligned_segment = {}
# Initially, when iterating through each alignment, 
# it is not known whether it is read1 or read2,
# so this is stored as a 'segment'.
aligner_to_current_read1_alignment = {}
aligner_to_current_read2_alignment = {}

# keep a dictionary for adjustments to ambiguous junctions based on splice site motifs to avoid reanalyzing the same junction
junctions_to_ambiguous_junctions = {}
# keep a dictionary for junctions to aligners to counts
junction_counts = {}
# keep a dictionary for aggregating info at the level of each junction across all supporting reads, including info on:
# junctions to counts, intron_coords_adjusted, count mismatch_near_junction, list US_perfect_matches,  list DS_perfect_matches, etc.
COMPASS_junction_statistics = {}

concordant_ungapped_counts = 0 # concordant means complete agreement across all aligners on the CIGAR
concordant_gapped_counts = 0
discordant_ungapped_counts = 0 # no aligner has an intron, and cigars don't agree
discordant_gapped_counts = 0 # at least one aligner has an intron, but cigars don't agree

min_score_has_no_intron_cigars_equal = {True: 0, False: 0}
min_score_has_annotated_intron_cigars_equal = {True: 0, False: 0}
min_score_has_annotated_intron_junctions_equal = {True: 0, False: 0}
min_score_has_unannotated_intron_cigars_equal = {True: 0, False: 0}
min_score_has_unannotated_intron_junctions_equal = {True: 0, False: 0}

read_pair_by_aligner_columns = [
    'read_pair_ID', 'aligner_has_min_score', 'min_score_has_annotated_intron',
    'min_score_has_no_intron', 'aligner', 'score',  'num_introns_found', 'perfect_gapped_alignment',
    'chrom', 'start', 'stop', 'five_SS', 'three_SS', 'RNA_strand', 'annotated_junction', 'ann_5SS',
    'ann_3SS', 'canonical_5SS', 'canonical_3SS', 'intron_size', 'intron_coords_adjusted',
    'mismatch_near_junction', 'US_perfect_matches', 'DS_perfect_matches', 'read1_flag',
                    'read1_chrom', 'read1_coord', 'read1_cigar', 'read1_NH', 'read2_flag', 'read2_chrom',
                        'read2_coord', 'read2_cigar', 'read2_NH']
header = '\t'.join(read_pair_by_aligner_columns) + '\n'
[e.write(header) for e in intron_found_to_alignments_outfiles.values()]
eof = False

alignment_keys = 'flag, chrom, coord, cigar, NH, adjusted_introns, perfect_gapped_alignment, num_introns_found, alignment_score, RNA_strand'.split(', ')          
# NH refers to number of hits, i.e. number of distinct alignments for a given read.
# NH will be 1 for a uniquely mapped read.
# for reads mapping to multiple locations, NH will be greater than 1 
# and the aligner will designate one of the reads as the 'primary' alignment.
# The primary alignment is the only one which is analyzed by COMPASS.

COMPASS_common_intron_info = 'chrom, adj_start, adj_stop, five_SS, three_SS, RNA_strand, num_amb_junctions, annotated_junction, ann_5SS, ann_3SS, canonical_5SS, canonical_3SS, intron_size'.split(', ')
COMPASS_sample_specific_intron_info = 'intron_coords_adjusted, mismatch_near_junction, US_perfect_matches, DS_perfect_matches, five_SS_Q_score, three_SS_Q_score'.split(', ')
COMPASS_info = COMPASS_common_intron_info + COMPASS_sample_specific_intron_info

adjusted_intron_keys = 'chrom, amb_start, amb_stop, five_SS, three_SS, RNA_strand, num_amb_junctions, annotated_junction, ann_5SS, ann_3SS, canonical_5SS, canonical_3SS, intron_size, intron_coords_adjusted, mismatch_near_junction , US_perfect_matches, DS_perfect_matches, five_SS_Q_score, three_SS_Q_score'.split(', ')
adjusted_intron_keys_length = len(adjusted_intron_keys)

def equivalent_alignments(read1_alignments, read2_alignments, idx_lst):
    return len(set([read1_alignments[i] for i in idx_lst])) == len(set([read2_alignments[i] for i in idx_lst])) == 1


elapsed_time = timeit.default_timer() - start_time
log_message = ('Genome fasta and intron annotations processed in ' + str(elapsed_time/60) + ' minutes\nStarting COMPASS alignment comparisons...\n')
print(log_message)
log_outfile.write(log_message)

start_time = timeit.default_timer()

while read_num <= reads_to_process and not eof:
    # If alignment file reaches the end during iteration, all reads have been processed.
    if read_num > 0 and read_num % READ_NUM_PROGRESS == 0:
        elapsed_time = timeit.default_timer() - start_time
        reads_per_minute = read_num / elapsed_time * 60
        minutes_for_100_million_reads = 10**8 / reads_per_minute
        log_message = (str(read_num) + ' read pairs processed by COMPASS in ' + str(elapsed_time/60) + ' minutes.\n' + \
        'This corresponds to ' + str(reads_per_minute) + ' reads per minute.\n' + \
        'Estimated time for 100 million reads is ' + str(minutes_for_100_million_reads) + ' minutes.\n')
        print(log_message)
        log_outfile.write(log_message)
        
        log_message = ('The COMPASS alignment comparison step is currently at read number: ' + str(read_num) + '\n' + \
        'concordant_ungapped_counts: ' + str(concordant_ungapped_counts) + '\n' + \
        'discordant_ungapped_counts: ' + str(discordant_ungapped_counts) + '\n' + \
        'concordant_gapped_counts: ' + str(concordant_gapped_counts) + '\n' + \
        'discordant_gapped_counts: ' + str(discordant_gapped_counts) + '\n' + \
        'equivalent alignments for all min score aligners when a min score has ungapped alignment: ' + str(min_score_has_no_intron_cigars_equal) + '\n' + \
        'equivalent alignments for all min score aligners when a min score has annotated junction: ' + str(min_score_has_annotated_intron_cigars_equal) + '\n' + \
        'equivalent junctions for all min score aligners when a min score has annotated junction: ' + str(min_score_has_annotated_intron_junctions_equal) + '\n' + \
        # Aligners can disagree on the details of the complete alignment \
        # but have equivalent junctions, for example in the representation \ 
        # of an indel or MNP near the end of a read, where the splice site is identical.
        'equivalent alignments for all min score aligners when a min score has unannotated junction: ' + str(min_score_has_unannotated_intron_cigars_equal) + '\n' + \
        'equivalent junctions found for all min score aligners when a min score has unannotated junction: ' + str(min_score_has_unannotated_intron_junctions_equal) + '\n')
        print(log_message)
        log_outfile.write(log_message)
    read_num += 1
    gapped_alignment_present = False
    quality_scores_of_current_read = {'R1':[], 'R2':[]}
    # for each aligner, gather the alignments for the current read number, if present, otherwise None will be used
    for aligner in aligners: 
        secondary_alignment = False
        aligner_to_current_read1_alignment[aligner] = None
        aligner_to_current_read2_alignment[aligner] = None
        while aligner_to_current_read_num[aligner] < read_num or secondary_alignment: 
            try:
                aligned_segment = next(aligners_to_bam[aligner])
                secondary_alignment = aligned_segment.is_secondary
                # only consider primary alignments from each aligner
            except StopIteration:
                eof = True
                break
            aligner_bam_read_num = int(aligned_segment.query_name.split('_')[0])
            if aligner_bam_read_num < aligner_to_current_read_num[aligner] :
                log_message = ('lower read number ' + str(aligner_bam_read_num) + ' found after higher read num ' + str(aligner_to_current_read_num[aligner]) + ' for aligner: ' + aligner)
                print(log_message)
                log_outfile.write(log_message)
                raise ValueError
            if not secondary_alignment:
                aligner_to_current_read_num[aligner] = aligner_bam_read_num
                aligner_to_current_aligned_segment[aligner] = aligned_segment
        if aligner_to_current_read_num[aligner] == read_num: 
        # Process only the current read number, 
        # as there is a chance during the last while loop that an aligner didn't report a mapping for a read.
        # It would therefore jump to the next read, which would need to be kept in memory for the next round of COMPASS.
            aligned_segment = aligner_to_current_aligned_segment[aligner]
            read1_first = aligned_segment.is_read1
            if read1_first:
                aligner_to_current_read1_alignment[aligner] = aligned_segment
                if aligned_segment.query_qualities != None:
                    quality_scores_of_current_read['R1'] = aligned_segment.query_qualities
            else:
                aligner_to_current_read2_alignment[aligner] = aligned_segment
                if aligned_segment.query_qualities != None:
                    quality_scores_of_current_read['R2'] = aligned_segment.query_qualities
            # After getting the first read in a pair for a given read num, it is necessary to get the other read from the pair.
            # possibly requiring going through several other secondary mappings before getting to the primary mapping of the other read.
            while aligner_to_current_read_num[aligner] == read_num and (aligner_to_current_read1_alignment[aligner] == None or aligner_to_current_read2_alignment[aligner] == None):    
                try:
                    aligned_segment = next(aligners_to_bam[aligner])
                except StopIteration:
                    eof = True
                    break
                aligner_to_current_read_num[aligner] = int(aligned_segment.query_name.split('_')[0])
                aligner_to_current_aligned_segment[aligner] = aligned_segment
            # This code works under the assumption for a given read pair, 
            # there are no more than two primary alignments, one for each read in the pair.  
            # This code allows for the possibility that one read of a pair has an unreported alignment,
            # in which case that read will be assigned a None value.
                if aligner_to_current_read_num[aligner] == read_num and not aligned_segment.is_secondary:
                    if read1_first and aligned_segment.is_read1:
                        # this should never print, except in cases where a read is split onto two lines 
                        # due to the first part of the read aligning 'on the wrong side'
                        log_message = (str(aligner) + ' has two adjacent read 1s for read pair number: ' + str(read_num) + '.')
                        print(log_message)
                        log_outfile.write(log_message + '\n')
                    if aligned_segment.is_read1:
                        aligner_to_current_read1_alignment[aligner] = aligned_segment
                        if aligned_segment.query_qualities != None:
                            quality_scores_of_current_read['R1'] = aligned_segment.query_qualities
                    else:
                        aligner_to_current_read2_alignment[aligner] = aligned_segment
                        if aligned_segment.query_qualities != None:
                            quality_scores_of_current_read['R2'] = aligned_segment.query_qualities

    # first process alignment scores to see if best score meets threshold for further processing
    R1_lst = []
    R2_lst = []
    alignment_scores = []
    num_intron_found_lst = []
    for aligner in aligners:
        R1, junctions_to_ambiguous_junctions = process_alignment(aligner_to_current_read1_alignment[aligner], \
            quality_scores_of_current_read['R1'], junctions_to_ambiguous_junctions, \
            genome_fasta, annotated_intron_df, annotated_introns, \
            CONSENSUS_5SS, PENALTIES_5SS, CONSENSUS_3SS, PENALTIES_3SS, \
            MISMATCH_DIST_FROM_JUNCTION_DISALLOWED, MIN_INTRON_LENGTH)
        R2, junctions_to_ambiguous_junctions = process_alignment(aligner_to_current_read2_alignment[aligner], \
            quality_scores_of_current_read['R2'], junctions_to_ambiguous_junctions, \
            genome_fasta, annotated_intron_df, annotated_introns, \
            CONSENSUS_5SS, PENALTIES_5SS, CONSENSUS_3SS, PENALTIES_3SS, \
            MISMATCH_DIST_FROM_JUNCTION_DISALLOWED, MIN_INTRON_LENGTH)
        R1_lst.append(R1)
        R2_lst.append(R2)
        alignment_scores.append(R1['alignment_score'] + R2['alignment_score'])
        num_intron_found_lst.append(len(set(R1['splice_sites'] + R2['splice_sites'])))
    best_score = min(alignment_scores)
    min_score_has_annotated_intron = False
    min_score_has_no_intron = False

    # initialize lists for info on each alignment, access later by idx using this construct: for idx in range(len(aligners))
    perfect_gapped_alignment_lst = []
    read1_alignments = []
    read2_alignments = []
    splice_junction_lst = []
    intron_lengths_lst = [] # for each read, take the maximum intron length
    comment_tag_lst = []
    alignment_type_lst = []
    num_annotated_introns_lst = []
    idx_without_mismatch_near_junction = []
    idx_with_best_score = []
    idx_best_score_with_no_intron = []
    idx_best_score_with_annotated_intron = []
    # need to loop through the aligners first to pre-compute:  
    # min_score_has_annotated_intron, min_score_has_no_intron, perfect_gapped_alignment_lst
    for idx in range(len(aligners)):
        aligner = aligners[idx]
        R1 = R1_lst[idx]
        R2 = R2_lst[idx]
        splice_junctions_in_read_pair = []
        intron_lengths = [0]
        annotated_introns_found = 0
        for e in R1, R2:
            if e['adjusted_introns'] != []:
            # values in adjusted_introns list: chrom, adj_start, adj_stop, five_SS, three_SS, \
            # RNA_strand, num_amb_junctions, annotated_junction, ann_5SS, ann_3SS, canonical_5SS, \
            # canonical_3SS, intron_size, intron_coords_adjusted, mismatch_near_junction, \
            # US_perfect_matches, DS_perfect_matches, five_SS_Q_score, three_SS_Q_score

            # chrom, adj_start, adj_stop, five_SS, three_SS, RNA_strand, num_amb_junctions, \
            # annotated_junction, ann_5SS, ann_3SS, canonical_5SS, canonical_3SS = adjusted_intron  

            # read 1 and read 2 info entered on separate lines if both contain a gapped alignment, 
            # later on these can be further processed and integrated by groupby operations

            # if a read contains more than one intron (quite rare), it is entered on multiple lines
                for intron in e['adjusted_introns']:
                    splice_junctions_in_read_pair.append(tuple([intron[e] for e in  'chrom, adj_start, adj_stop, five_SS, three_SS, RNA_strand'.split(', ')] ))
                    if intron['annotated_junction']:
                        annotated_introns_found += 1
                    if not intron['mismatch_near_junction']:
                        idx_without_mismatch_near_junction.append(idx)  
                    intron_lengths.append(intron['intron_size'])
        intron_lengths_lst.append(max(intron_lengths))  
        splice_junction_lst.append(tuple(splice_junctions_in_read_pair))
        num_annotated_introns_lst.append(annotated_introns_found)
        if alignment_scores[idx] == best_score:
            idx_with_best_score.append(idx)
            if annotated_introns_found > 0:
                min_score_has_annotated_intron = True
                idx_best_score_with_annotated_intron.append(idx)
            if splice_junctions_in_read_pair == []:
                min_score_has_no_intron = True
                idx_best_score_with_no_intron.append(idx)
        # if introns are present, is there a perfect gapped alignment present in the read pair?
        perfect_gapped_alignment_lst.append(R1['perfect_gapped_alignment'] or R2['perfect_gapped_alignment'])  

        if num_intron_found_lst[idx] == 0:
            alignment_type = 'UGA'  
        else:
            if num_annotated_introns_lst[idx] > 0:
                alignment_type = 'ASJ'
            else:
                alignment_type = 'USJ'
        alignment_type_lst.append(alignment_type)

        comment_tag = [aligner, R1['chrom'], R1['coord'], R1['cigar'], R2['chrom'], R2['coord'], R2['cigar'], alignment_scores[idx], alignment_type_lst[idx]]
        comment_tag_lst.append(','.join([str(e) for e in comment_tag]) )

        # write the complete comparison table for all read alignments (big table = num reads X num aligners) 
        for idx in range(len(aligners)):            
            read_pair_info = [read_num, alignment_scores[idx] == best_score, \
                min_score_has_annotated_intron, min_score_has_no_intron] + [f[idx] \
                    for f in (aligners, alignment_scores, num_intron_found_lst, perfect_gapped_alignment_lst)]     
            R1 = R1_lst[idx]
            R2 = R2_lst[idx]  
            # alignment_keys = 'flag, chrom, coord, cigar, NH, adjusted_introns, perfect_gapped_alignment, num_introns_found, alignment_score'.split(', ')          
            read_info_keys = 'flag, chrom, coord, cigar, NH'.split(', ')
            read1_info = [R1[e] for e in read_info_keys]
            read2_info = [R2[e] for e in read_info_keys]
            # remove flag and NH for the alignment comparison
            read1_alignments.append(tuple(read1_info[1:4]))
            read2_alignments.append(tuple(read2_info[1:4]))
            if max(num_intron_found_lst) > 0: 
                for e in R1, R2:
                    if e['adjusted_introns'] != []:
                    # chrom, adj_start, adj_stop, five_SS, three_SS, RNA_strand, num_amb_junctions, annotated_junction, ann_5SS, ann_3SS, canonical_5SS, canonical_3SS = adjusted_intron  
                    # read 1 and read 2 info entered on separate lines if both contain a gapped alignment, 
                    # later on these can be further processed and integrated by groupby operations
                    # if a read contains more than intron (quite rare), it is entered on multiple lines
                        for intron in e['adjusted_introns']:
                            # # # adjusted_introns --> chrom, amb_start, amb_stop, five_SS, three_SS, RNA_strand, num_amb_junctions, annotated_junction, ann_5SS, ann_3SS, canonical_5SS, canonical_3SS, intron_size, intron_coords_adjusted, mismatch_near_junction , US_perfect_matches, DS_perfect_matches, five_SS_Q_score, three_SS_Q_score
                            output = read_pair_info + list(intron.values()) + read1_info + read2_info
                            output = [str(e) for e in output]
                            # when a read pair had its best scoring alignment contain an annotated intron, 
                            # all alignments are written to the annotated intron comparison table
                            intron_found_to_alignments_outfiles[min_score_has_annotated_intron].write('\t'.join(output) + '\n')          
                if num_intron_found_lst[idx] == 0: # this line will only execute if the preceding for loop does nothing
                    # if no gapped alignments found, enter NA values for the intron info
                    output = read_pair_info + ['NA']*adjusted_intron_keys_length + read1_info + read2_info
                    output = [str(e) for e in output]
                    intron_found_to_alignments_outfiles[min_score_has_annotated_intron].write('\t'.join(output) + '\n')
                    # ungapped alignments are only written to the gapped outfiles when the best scoring aligner had a gapped alignment,
                    # and other aligners have ungapped alignmnets with inferior scores
        
        # write the optimal alignment to the COMPASS bam file   
        all_aligners_agree_on_cigar = False
        if len(set(read1_alignments)) == len(set(read2_alignments)) == 1:  
            # check if all alignments are identical based on: flag, chrom, coord, cigar
            all_aligners_agree_on_cigar = True
            if min(num_intron_found_lst) > 0:
                concordant_gapped_counts += 1
            else:
                concordant_ungapped_counts += 1
        elif max(num_intron_found_lst) > 0:
            discordant_gapped_counts += 1  # at least one aligner has an intron, but cigars don't agree
        else:
            discordant_ungapped_counts += 1 # no aligner has an intron, and cigars don't agree
            
        num_aligners_with_best_score = len(idx_with_best_score)
        min_score_unannotated_intron_junctions_agree = True
        if min_score_has_no_intron:
            # need to consider what happens when an intron greater than max allowed length is present   
            aligner_idx_selected = random.choice(idx_best_score_with_no_intron)
            min_score_has_no_intron_cigars_equal[equivalent_alignments(read1_alignments, read2_alignments, idx_best_score_with_no_intron)]  += 1
        elif min_score_has_annotated_intron:
            aligner_idx_selected = random.choice(idx_best_score_with_annotated_intron)
            min_score_has_annotated_intron_cigars_equal[equivalent_alignments(read1_alignments, read2_alignments, idx_best_score_with_annotated_intron)]  += 1
            min_score_has_annotated_intron_junctions_equal[len(set([splice_junction_lst[i] for i in idx_best_score_with_annotated_intron])) == 1] += 1  
        else:
            min_score_has_unannotated_intron_cigars_equal[equivalent_alignments(read1_alignments, read2_alignments, idx_with_best_score)] += 1
            min_score_unannotated_intron_junctions_agree = (len(set([splice_junction_lst[i] for i in idx_with_best_score])) == 1)
            min_score_has_unannotated_intron_junctions_equal[min_score_unannotated_intron_junctions_agree] += 1
            if not min_score_unannotated_intron_junctions_agree:
                junction_min_score_disagreement_outfile.write(str(aligners) + '\n' + str(read1_alignments) + '\n' + str(read2_alignments) + '\n' + str(splice_junction_lst) + '\n' + str(alignment_scores) + '\n'*2)
            min_intron_length_with_best_score = min([intron_lengths_lst[i] for i in idx_with_best_score])
            idx_with_best_score_and_min_intron_size = [idx for idx in range(len(intron_lengths_lst)) if intron_lengths_lst[idx] == min_intron_length_with_best_score and idx in idx_with_best_score]
            tie_breaker_lst = [i for i in idx_without_mismatch_near_junction if i in idx_with_best_score_and_min_intron_size]
            if tie_breaker_lst != []:
                aligner_idx_selected = random.choice(tie_breaker_lst)
            else:
                aligner_idx_selected = random.choice(idx_with_best_score_and_min_intron_size)
        aligner_selected = aligners[aligner_idx_selected]    
        alignment_type = alignment_type_lst[aligner_idx_selected]
        comment_tag_out = ';'.join(comment_tag_lst) + '/' + str(num_aligners_with_best_score) + ',' + alignment_type
        aligners_with_best_score = [aligners[idx] for idx in idx_with_best_score]
        for read_alignment in aligner_to_current_read1_alignment, aligner_to_current_read2_alignment:
            # if read_alignment != None:
            # None type error here???
            read_alignment[aligner_selected].set_tag('PG', ','.join(aligners_with_best_score), 'Z')
            read_alignment[aligner_selected].set_tag('CO', comment_tag_out, 'Z')
        RNA_strand = R1_lst[aligner_idx_selected]['RNA_strand']
        [COMPASS_bam_outfile.write(e[aligner_selected]) for e in (aligner_to_current_read1_alignment, aligner_to_current_read2_alignment) if e != None]
        if alignment_type == 'UGA':
            if RNA_strand == '+':
                [COMPASS_UGA_plus_strand_bam_outfile.write(e[aligner_selected]) for e in (aligner_to_current_read1_alignment, aligner_to_current_read2_alignment) if e != None]
            else:
                [COMPASS_UGA_minus_strand_bam_outfile.write(e[aligner_selected]) for e in (aligner_to_current_read1_alignment, aligner_to_current_read2_alignment) if e != None]
        if alignment_type == 'USJ':
            [COMPASS_USJ_bam_outfile.write(e[aligner_selected]) for e in (aligner_to_current_read1_alignment, aligner_to_current_read2_alignment) if e != None]
        if alignment_type == 'ASJ':
            [COMPASS_ASJ_bam_outfile.write(e[aligner_selected]) for e in (aligner_to_current_read1_alignment, aligner_to_current_read2_alignment) if e != None]

        # write the optimal alignment to the COMPASS splice junction counts file  
        if alignment_type == 'USJ' or alignment_type == 'ASJ':
            R1_selected = R1_lst[aligner_idx_selected]
            R2_selected = R2_lst[aligner_idx_selected]
            perfect_gapped_alignment = perfect_gapped_alignment_lst[aligner_idx_selected]
            num_introns_found = num_intron_found_lst[aligner_idx_selected]
            
            introns_to_attributes = {}
            for e in R1_selected, R2_selected:
                # # # adjusted_introns --> chrom, amb_start, amb_stop, five_SS, three_SS, RNA_strand, num_amb_junctions, annotated_junction, ann_5SS, ann_3SS, canonical_5SS, canonical_3SS, intron_size, intron_coords_adjusted, mismatch_near_junction, US_perfect_matches, DS_perfect_matches, five_SS_Q_score, three_SS_Q_score
                # chrom, amb_start, amb_stop, five_SS, three_SS, RNA_strand, num_amb_junctions, annotated_junction, ann_5SS, ann_3SS, canonical_5SS, canonical_3SS, intron_size,
                # first 10 are invariate attributes of the coords, last 4 could be different (last 2 almost certainly different)
                # intron_coords_adjusted unlikely to be different for read pairs overlappign the same junction
                # intron_coords_adjusted, mismatch_near_junction, US_perfect_matches, DS_perfect_matches, five_SS_Q_score, three_SS_Q_score
                for intron in e['adjusted_introns']:
                    chrom_coords = tuple(intron[e] for e in COMPASS_common_intron_info)
                    if chrom_coords not in introns_to_attributes:
                        introns_to_attributes[chrom_coords] = {k:v for (k,v) in intron.items() if k in COMPASS_sample_specific_intron_info}
                        # intron_coords_adjusted, mismatch_near_junction, US_perfect_matches, DS_perfect_matches, five_SS_Q_score, three_SS_Q_score
                    # if read pairs overlap on same junction, integrate the data
                    else:
                        # If both read 1 and 2 have the same junction, they should agree on whether the intron coords were adjusted.
                        # If both read 1 and 2 have the same junction, and only one shows mismatches near the junction, assume those are sequencing errors.
                        introns_to_attributes[chrom_coords]['intron_coords_adjusted'] = intron['intron_coords_adjusted'] and introns_to_attributes[chrom_coords]['intron_coords_adjusted']
                        introns_to_attributes[chrom_coords]['mismatch_near_junction'] = intron['mismatch_near_junction'] or introns_to_attributes[chrom_coords]['mismatch_near_junction']
                        introns_to_attributes[chrom_coords]['US_perfect_matches'] = max(intron['US_perfect_matches'],introns_to_attributes[chrom_coords]['US_perfect_matches'])
                        introns_to_attributes[chrom_coords]['DS_perfect_matches'] = max(intron['DS_perfect_matches'], introns_to_attributes[chrom_coords]['DS_perfect_matches'])
                        introns_to_attributes[chrom_coords]['five_SS_Q_score'] = max(intron['five_SS_Q_score'],introns_to_attributes[chrom_coords]['five_SS_Q_score'])
                        introns_to_attributes[chrom_coords]['three_SS_Q_score'] = max(intron['three_SS_Q_score'], introns_to_attributes[chrom_coords]['three_SS_Q_score'])
            for chrom_coords in introns_to_attributes:    
                # update COMPASS statistics dict
                #  values = read.flag, read.reference_name, read.reference_start, read.cigarstring.replace('D', 'N').replace('M', 'X'), read.get_tag('NH'), adjusted_introns, perfect_gapped_alignment, num_introns_found, alignment_score
                if chrom_coords not in COMPASS_junction_statistics:
                    COMPASS_junction_statistics[chrom_coords] = {}
                    COMPASS_junction_statistics[chrom_coords]['US_perfect_matches'] = []
                    COMPASS_junction_statistics[chrom_coords]['DS_perfect_matches'] = []
                    COMPASS_junction_statistics[chrom_coords]['five_SS_Q_score'] = []
                    COMPASS_junction_statistics[chrom_coords]['three_SS_Q_score'] = []
                    COMPASS_junction_statistics[chrom_coords]['alignment_scores'] = []
                    COMPASS_junction_statistics[chrom_coords]['num_introns_found'] = []                   
                    COMPASS_junction_statistics[chrom_coords]['read_1_coord_cigars'] = []
                    COMPASS_junction_statistics[chrom_coords]['read_2_coord_cigars'] = []
                    COMPASS_junction_statistics[chrom_coords]['read_pair_cigar_combos'] = set([])
                    COMPASS_junction_statistics[chrom_coords]['intron_coords_adjusted_counts'] = 0
                    COMPASS_junction_statistics[chrom_coords]['mismatch_near_junction_counts'] = 0
                    COMPASS_junction_statistics[chrom_coords]['perfect_gapped_alignment_counts'] = 0
                    COMPASS_junction_statistics[chrom_coords]['min_score_unannotated_intron_junctions_disagree'] = 0
                COMPASS_junction_statistics[chrom_coords]['US_perfect_matches'].append(introns_to_attributes[chrom_coords]['US_perfect_matches'])
                COMPASS_junction_statistics[chrom_coords]['DS_perfect_matches'].append(introns_to_attributes[chrom_coords]['DS_perfect_matches'])
                COMPASS_junction_statistics[chrom_coords]['five_SS_Q_score'].append(introns_to_attributes[chrom_coords]['five_SS_Q_score'])
                COMPASS_junction_statistics[chrom_coords]['three_SS_Q_score'].append(introns_to_attributes[chrom_coords]['three_SS_Q_score'])
                COMPASS_junction_statistics[chrom_coords]['alignment_scores'].append(best_score)
                COMPASS_junction_statistics[chrom_coords]['num_introns_found'].append(num_introns_found)
                COMPASS_junction_statistics[chrom_coords]['read_1_coord_cigars'].append((R1_selected['coord'], R1_selected['cigar']))
                COMPASS_junction_statistics[chrom_coords]['read_2_coord_cigars'].append((R2_selected['coord'], R2_selected['cigar']))
                COMPASS_junction_statistics[chrom_coords]['read_pair_cigar_combos'].add((R1_selected['coord'], R1_selected['cigar'], R2_selected['coord'], R2_selected['cigar']))
                if introns_to_attributes[chrom_coords]['intron_coords_adjusted']:
                    COMPASS_junction_statistics[chrom_coords]['intron_coords_adjusted_counts'] += 1
                if introns_to_attributes[chrom_coords]['mismatch_near_junction']:
                    COMPASS_junction_statistics[chrom_coords]['mismatch_near_junction_counts'] += 1
                if perfect_gapped_alignment:
                    COMPASS_junction_statistics[chrom_coords]['perfect_gapped_alignment_counts'] += 1
                if not min_score_unannotated_intron_junctions_agree:
                    COMPASS_junction_statistics[chrom_coords]['min_score_unannotated_intron_junctions_disagree'] += 1
                
                # update counts for each aligner on each junction
                # even if an aligner is not selected by COMPASS, give credit for finding the COMPASS-selected intron
                if chrom_coords not in junction_counts:
                    junction_counts[chrom_coords] = {}
                    for aligner in aligners:
                        junction_counts[chrom_coords][aligner] = 0
                junction_counts[chrom_coords]['COMPASS'] = junction_counts[chrom_coords].get('COMPASS', 0) + 1
                chrom, start, stop = chrom_coords[:3]
                for idx in range(len(aligners)):
                    aligner = aligners[idx]
                    R1 = R1_lst[idx]
                    R2 = R2_lst[idx]
                    intron_found_by_this_aligner = False
                    for e in R1, R2:
                        for intron in e['adjusted_introns']:
                            chrom_coords_to_check = tuple([intron[e] for e in COMPASS_common_intron_info])
                            if chrom_coords_to_check == chrom_coords:
                                intron_found_by_this_aligner = True
                            # if chrom == 'chrX' and start > 73800 and stop < 74300:
                            #     print(aligner, intron)
                    if intron_found_by_this_aligner:
                        junction_counts[chrom_coords][aligner] += 1

[aligners_to_bam[aligner].close() for aligner in aligners]

[e.close() for e in (COMPASS_bam_outfile, COMPASS_UGA_plus_strand_bam_outfile, COMPASS_UGA_minus_strand_bam_outfile, COMPASS_ASJ_bam_outfile, COMPASS_USJ_bam_outfile, annotated_gapped_alignments_outfile, unannotated_gapped_alignments_outfile, junction_min_score_disagreement_outfile)] 

elapsed_time = timeit.default_timer() - start_time
reads_per_minute = read_num / elapsed_time * 60
minutes_for_100_million_reads = 10*8 / reads_per_minute
log_message = (str(read_num) + ' read pairs processed by COMPASS in ' + str(elapsed_time/60) + ' minutes.\n\
    This corresponds to ' + str(reads_per_minute) + ' read pairs per minute.\n')
print(log_message)
log_outfile.write(log_message)

log_message = ('COMPASS finished comparing alignments for ' + str(read_num) + ' reads.\n' + \
'concordant_ungapped_counts: ' + str(concordant_ungapped_counts) + '\n' + \
'discordant_ungapped_counts: ' + str(discordant_ungapped_counts) + '\n' + \
'concordant_gapped_counts: ' + str(concordant_gapped_counts) + '\n' + \
'discordant_gapped_counts: ' + str(discordant_gapped_counts) + '\n' + \
'equivalent alignments for all min score aligners when a min score has ungapped alignment: ' + str(min_score_has_no_intron_cigars_equal) + '\n' + \
'equivalent alignments for all min score aligners when a min score has annotated junction: ' + str(min_score_has_annotated_intron_cigars_equal) + '\n' + \
'equivalent junctions for all min score aligners when a min score has annotated junction: ' + str(min_score_has_annotated_intron_junctions_equal) + '\n' + \
'equivalent alignments for all min score aligners when a min score has unannotated junction: ' + str(min_score_has_unannotated_intron_cigars_equal) + '\n' + \
'equivalent junctions found for all min score aligners when a min score has unannotated junction: ' + str(min_score_has_unannotated_intron_junctions_equal) + '\nStarting the final junction analysis step...\n')
print(log_message)
log_outfile.write(log_message)

sample_info_header = ['sample_name', 'total_reads']
junction_header = 'chrom, adj_start, adj_stop, five_SS, three_SS, RNA_strand, num_amb_junctions, annotated_junction, ann_5SS, ann_3SS, canonical_5SS, canonical_3SS, intron_size'.split(', ')
counts_header = ['COMPASS_counts'] + [e + '_counts' for e in aligners]

stats_header_part_1 = 'min_score_unannotated_intron_junctions_disagree, intron_coords_adjusted_counts, mismatch_near_junction_counts, perfect_gapped_alignment_counts'.split(', ')
stats_header_part_2 = 'US_perfect_matches, DS_perfect_matches, five_SS_Q_score, three_SS_Q_score, alignment_scores, num_introns_found, read_1_coord_cigars, read_2_coord_cigars'.split(', ')
stats_header_part_3 = 'num_unique_read_pair_cigars, max_US_perfect_matches, max_DS_perfect_matches, unique_US_perfect_matches, unique_DS_perfect_matches, best_alignment_score, mean_alignment_score, median_alignment_score'.split(', ')
stats_header_part_4 = 'mean_five_SS_Q_score, mean_three_SS_Q_score'.split(', ')
COMPASS_header = '\t'.join(sample_info_header + junction_header + counts_header + stats_header_part_1 + stats_header_part_2 + stats_header_part_3 + stats_header_part_4) + '\n'
splice_junction_counts_outfile.write(COMPASS_header)

start_time = timeit.default_timer()
junctions_processed = 0
sorted_chrom_coords = sorted(COMPASS_junction_statistics.keys())
total_junctions_to_process = len(sorted_chrom_coords)
for chrom_coords in sorted_chrom_coords:
    junctions_processed += 1
    if junctions_processed % 100 == 0:
        elapsed_time = timeit.default_timer() - start_time
        junctions_per_minute = junctions_processed / elapsed_time * 60
        remaining_junctions_to_process = total_junctions_to_process - junctions_processed
        estimated_minutes_remaining = remaining_junctions_to_process / junctions_per_minute
        log_message = (str(junctions_processed) + ' junctions processed by COMPASS in ' + str(elapsed_time/60) + ' minutes.\n' + \
        'This corresponds to ' + str(junctions_per_minute) + ' junctions per minute.\n' +  \
        'Finishing the remaining set of ' + str(remaining_junctions_to_process) + \
        ' junctions should take ' + str(estimated_minutes_remaining) + ' minutes.\n')
        print(log_message)
        log_outfile.write(log_message)
    chrom, adj_start, adj_stop, five_SS, three_SS, RNA_strand = chrom_coords[:6]

    sample_info_output = [sample_name, str(read_num)]
    junction_output = [str(e) for e in chrom_coords]
    counts_output = [str(junction_counts[chrom_coords]['COMPASS'])] + [str(junction_counts[chrom_coords][e]) for e in aligners]

    stats_output_part_1 = [str(COMPASS_junction_statistics[chrom_coords][e]) for e in stats_header_part_1]
    stats_output_part_2 = [str(count_frequency(COMPASS_junction_statistics[chrom_coords][e], NUM_ITEMS_TO_REPORT)) for e in stats_header_part_2]    
    
    US_perfect_matches_set = set(COMPASS_junction_statistics[chrom_coords]['US_perfect_matches'])
    DS_perfect_matches_set = set(COMPASS_junction_statistics[chrom_coords]['DS_perfect_matches'])

    stats_output_part_3_raw = [len(COMPASS_junction_statistics[chrom_coords]['read_pair_cigar_combos']), \
        max(US_perfect_matches_set), max(DS_perfect_matches_set), len(US_perfect_matches_set), \
            len(DS_perfect_matches_set), min(COMPASS_junction_statistics[chrom_coords]['alignment_scores']), \
                statistics.mean(COMPASS_junction_statistics[chrom_coords]['alignment_scores']), \
                    statistics.median(COMPASS_junction_statistics[chrom_coords]['alignment_scores']) ]
    stats_output_part_3 = [str(e) for e in stats_output_part_3_raw]
    
    stats_output_part_4_raw = [statistics.mean(COMPASS_junction_statistics[chrom_coords]['five_SS_Q_score']), \
        statistics.mean(COMPASS_junction_statistics[chrom_coords]['three_SS_Q_score']) ]
    stats_output_part_4 = [str(e) for e in stats_output_part_4_raw]
    
    COMPASS_output = '\t'.join(sample_info_output + junction_output + counts_output + stats_output_part_1 + stats_output_part_2 + stats_output_part_3 + stats_output_part_4)  + '\n'
    splice_junction_counts_outfile.write(COMPASS_output)

splice_junction_counts_outfile.close()

elapsed_time = timeit.default_timer() - start_time
junctions_per_minute = junctions_processed / elapsed_time * 60
remaining_junctions_to_process = total_junctions_to_process - junctions_processed
estimated_minutes_remaining = remaining_junctions_to_process / junctions_per_minute
log_message = (str(junctions_processed) + ' junctions processed by COMPASS in ' + str(elapsed_time/60) + ' minutes.\n' + \
'This corresponds to ' + str(junctions_per_minute) + ' junctions per minute.\n' +  \
'COMPASS junction level analysis is complete.\nStarting the BAM sorting and indexing step...\n')
print(log_message)
log_outfile.write(log_message)

start_time = timeit.default_timer()

pysam.sort('-o', COMPASS_UGA_plus_strand_sorted_bam_filename, '-@', str(NUM_THREADS), COMPASS_UGA_plus_strand_bam_filename)
pysam.sort('-o', COMPASS_UGA_minus_strand_sorted_bam_filename, '-@', str(NUM_THREADS), COMPASS_UGA_minus_strand_bam_filename)
pysam.sort('-o', COMPASS_USJ_sorted_bam_filename, '-@', str(NUM_THREADS), COMPASS_USJ_bam_filename)
pysam.sort('-o', COMPASS_ASJ_sorted_bam_filename, '-@', str(NUM_THREADS), COMPASS_ASJ_bam_filename)
pysam.sort('-o', COMPASS_sorted_bam_filename, '-@', str(NUM_THREADS), COMPASS_bamfilename)

for sorted_bam in (COMPASS_UGA_plus_strand_sorted_bam_filename, COMPASS_UGA_minus_strand_sorted_bam_filename, COMPASS_USJ_sorted_bam_filename, COMPASS_ASJ_sorted_bam_filename, COMPASS_sorted_bam_filename):
    pysam.index(sorted_bam)

elapsed_time = timeit.default_timer() - start_time
log_message = ('BAM sorted and indexed in ' + str(elapsed_time/60) + ' minutes\n')
print(log_message)
log_outfile.write(log_message)

log_message = 'The COMPASS program compare_splice_junctions_from_multiple_aligners.py finished successfully.\n'
print(log_message)
log_outfile.write(log_message)
log_outfile.close()
