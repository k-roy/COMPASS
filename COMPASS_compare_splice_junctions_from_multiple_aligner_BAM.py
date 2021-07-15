# -*- coding: utf-8 -*-
"""
updated on Friday, July 2, 2021
@author: kevinroy

command line usage: python compare_splice_junctions_from_multiple_aligner_BAM.py GENOME_FASTA GFF_ANNOTATIONS ALIGNMENTS_ALIGNERS JUNCTION_STATISTICS_OUTFILENAME INTEGRATED_ALIGNMENT_SUMMARY_OUTFILENAME INTEGRATED_ALIGNMENT_SAM_OUTFILENAME  FILTERED_SAM_DIR SAMPLE_PREFIX

takes arbitrary number of BAM files containing alignments for a single sample, and outputs the optimal alignments as a single BAM file

arguments:
-GENOME_FASTA
-GFF_ANNOTATIONS
-OUT_DIR
-ALIGNMENT_DIR
-comma separated aligner names

SRR5582776_subsampled   bbmap
SRR5582776_subsampled   STAR_default
SRR5582776_subsampled   STAR_noncanonical
SRR5582776_subsampled   HISAT2_default
SRR5582776_subsampled   HISAT2_noncanonical

output files:
COMPASS table, with separate rows for each read pair * aligner combination, and columns describing info on each read pair
separate bam files for ungapped, annotated, and unannotated splicing junctions

requirements:
-BAM must be name sorted where read names are chronologically numbered prior to mapping

-CIGAR format must be N for intron gap, M for match, S for soft-clip, 
I for insertion, D for deletion, and X for mismatch

what this script does:
-selects the optimal alignment based on user-defined criteria for penalties 
for each CIGAR operation, with the option of breaking ties in two scenarios: 
    -reads map without gapped alignment favored over gapped alignment
    -reads mapped to annotated junctions favored over unannotated junctions
-indicates which aligners agree on the selected junction for each read
-writes an "integrated" junction file with summary statistics on agreement between aligners for each junction

cd /mnt/mindrinos/kevinroy/projects/COMPASS/
echo '' > COMPASS_compare_splice_junctions_from_multiple_aligner_BAM.py
nano COMPASS_compare_splice_junctions_from_multiple_aligner_BAM.py

"""

import pysam
import pandas as pd
from scipy import stats
import statistics
import os
import random
import operator
random.seed(1234)

# from sys import argv
# ALIGNMENTS_DIR = argv[1]
# OUT_DIR = argv[2]
# FASTA = argv[3]
# INTRONS_FILE = argv[4]
# ALIGNERS_FILE = argv[5]
# sample_name = argv[6]
# sample_suffix = argv[7]
# reads_to_process = int(argv[8])
# regions_to_mask = int(argv[9]) ## block rRNA region from processing as the computation appears to stall there often

COMPASS_DIR = '/u/project/guillom/kevinh97/COMPASS/' # '/mnt/mindrinos/kevinroy/projects/COMPASS/' #
GENOME_DIR = COMPASS_DIR + 'S288C_reference_genome_R64-2-1_20150113/'
PROCESSED_DIR = COMPASS_DIR + 'processed_data/'
ALIGNMENTS_DIR = PROCESSED_DIR + 'alignments/'
OUT_DIR = ALIGNMENTS_DIR + 'COMPASS_integration_full_output/'# 'COMPASS_integration_from_USJ/'
mkdir $OUT_DIR
# os.mkdir(OUT_DIR)
FASTA ='/u/project/guillom/kevinh97/COMPASS/S288C_reference_genome_R64-2-1_20150113/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa'
INTRONS_FILE = '/u/project/guillom/kevinh97/COMPASS/S288C_reference_genome_R64-2-1_20150113/saccharomyces_cerevisiae_R64-2-1_20150113_introns_no_chr.tsv'
ALIGNERS_FILE = '/u/project/guillom/kevinh97/sample_aligner_info.tsv'
sample_name = 'RRP6_1SC_subsampled' # 'SRR5582776'
sample_suffix = '_COMPASS_name_sorted.bam' # '_COMPASS_USJ.bam'
reads_to_process = 10**7

MINIMUM_INTRON_LENGTH = 20
MAXIMUM_INTRON_LENGTH = 2000
MAX_EDIT_DIST_TO_CONSIDER = 8 ## this is the max sum of bp involved in mismatches and indels in both reads (not including introns) to consider
MISMATCH_DIST_FROM_JUNCTION_DISALLOWED = 10

CONSENSUS_5SS = ['G', 'T', 'A', 'TAC', 'G', 'TAC']
PENALTIES_5SS = [ 4,   3,   1,   1,    2,   1]

CONSENSUS_3SS = ['CT', 'A', 'G']
PENALTIES_3SS = [1,     3,   3]

READ_NUM_PROGRESS = 10**5
NUM_ITEMS_TO_REPORT = 10

try:
    R64
except:    
    R64 =  pysam.FastaFile(FASTA) ## load_genome(GENOME_FASTA)

aligner_to_current_read_num = {}
aligners_to_bam = {}

aligner_df = pd.read_csv(ALIGNERS_FILE, sep = '\t') 
aligner_df.columns
aligners = []
for index, row in aligner_df.iterrows():
    aligner = row['aligner']
    aligners.append(aligner)
    bam_filename = ALIGNMENTS_DIR + aligner + '/' + sample_name + sample_suffix
    aligners_to_bam[aligner] = pysam.AlignmentFile(bam_filename, "rb")
    aligner_to_current_read_num[aligner] = 0

annotated_intron_df = pd.read_csv(INTRONS_FILE, sep = '\t')
## adjust 1-based introns to 0-based python coords
annotated_intron_df['adjusted_start'] = annotated_intron_df['start'] - 1
annotated_intron_df['adjusted_end'] = annotated_intron_df['end'] - 1
annotated_junctions = set([])
ambiguous_annotated_junctions = set([])
junctions_to_ambiguous_junctions = {}
annotated_junc_to_type = {}
for index, row in annotated_intron_df.iterrows():
    chrom = row['chrom']
    start = row['adjusted_start']
    end = row['adjusted_end']
    intron = (chrom, start, end)
    annotated_junctions.add( intron )
    annotated_junc_to_type[intron] = row['intron_type']
    idx = 1
    ambiguous = False
  #  while R64[chrom][start-idx] == R64[chrom][stop-idx]:
    while R64.fetch(chrom, start-idx, start-idx+1) == R64.fetch(chrom, end-idx+1, end-idx+2):
        ambiguous = True
        ambiguous_intron = (chrom, start-idx, end-idx)
        annotated_junctions.add( ambiguous_intron )
        annotated_junc_to_type[ ambiguous_intron]  = row['intron_type']
        idx += 1
    idx = 0
   # while R64[chrom][start-idx] == R64[chrom][stop-idx]:
    while R64.fetch(chrom, start-idx, start-idx+1) == R64.fetch(chrom, end-idx+1, end-idx+2):
        ambiguous = True
        ambiguous_intron = (chrom, start-idx, end-idx)
        annotated_junctions.add( ambiguous_intron )
        annotated_junc_to_type[ ambiguous_intron] = row['intron_type']
        idx -= 1    
    if ambiguous:
        ambiguous_annotated_junctions.add(intron)
print('number of annotated introns:', len(annotated_intron_df.index))
print('number of annotated introns with identical nt upstream or downstream of junctions:', len(ambiguous_annotated_junctions))
print('total number of intron coords mapping to annotated junctions:',len(annotated_junctions) )

keys = 'flag, chrom, coord, cigar, NH, adjusted_introns, perfect_gapped_alignment, num_introns_found, alignment_score'.split(', ')          

def parse_cigar(cigar):
    '''
    input: CIGAR string from SAM alignment
    output: list of 2 item tuples: mapped segment length, CIGAR operation
    '''
    total_cigar_operation_bp = {'N':0, 'X':0, 'I':0, 'D':0, 'S':0, '=':0, 'M':0 }
    segment_length = ''
    cigartuples = []
    for char in cigar:
        if char in '0123456789':
            segment_length += char
        else:
            operation = char
            mapped_segment_length = int(segment_length)
            ## convert D to N for BBMap, as it does not use N for introns
            if mapped_segment_length >= MINIMUM_INTRON_LENGTH and operation == 'D':
                operation = 'N'
            elif mapped_segment_length < MINIMUM_INTRON_LENGTH and operation == 'N':
                operation = 'D'
            segment_length = ''
            cigartuples.append( (mapped_segment_length, operation) )
            total_cigar_operation_bp[operation] += mapped_segment_length
    return cigartuples, total_cigar_operation_bp

def hamming(a, b):
    return len([i for i in filter(lambda x: x[0] != x[1], zip(a, b))])
    
def five_SS_score(query):
    return sum(PENALTIES_5SS[i] for i in range(len(query)) if query[i] not in CONSENSUS_5SS[i])

def three_SS_score(query):
    return sum(PENALTIES_3SS[i] for i in range(len(query)) if query[i] not in CONSENSUS_3SS[i])

def rc(seq):
    return seq.translate( str.maketrans("ACTG", "TGAC") )[::-1]

def count_frequency(lst, NUM_ITEMS_TO_REPORT):
    counts = {}
    for i in lst: counts[i] = counts.get(i, 0) + 1
    sort_counts= sorted(counts.items(), key=operator.itemgetter(1))
    return sort_counts[:NUM_ITEMS_TO_REPORT]

def adjust_ambiguous_junctions(chrom, start, stop, RNA_strand):
    '''
    takes intron coordinates
    checks for ambiguous junction based on nt upstream and downstream of exon-intron and intron-exon junctions
    if ambiguous, adjusts based on closest match to 5SS consensus with 3SS adding smaller scoring component
    ## max mismatch for 5SS = 10, max mismatch for 3SS = 6
    returns adjusted coordinates and other info on the intron: 
    five_SS, three_SS, RNA_strand, num_amb_junctions, annotated_junction, ann_5SS, ann_3SS, canonical_5SS, canonical_3SS, intron_size, intron_coords_adjusted
    '''
    junctions = [ (chrom, start, stop) ]
    idx = 1
    # while start-idx > 0 and R64[chrom][start-idx] == R64[chrom][stop-idx+1]:
    ## don't consider junctions with 10 identical nt upstream or downstream
    while idx < 10 and start-idx > 0 and R64.fetch(chrom, start-idx, start-idx+1) == R64.fetch(chrom, stop-idx+1, stop-idx+2):
        #print(chrom, start, stop, idx, R64.fetch(chrom, start-idx, start-idx+1) , R64.fetch(chrom, stop-idx+1, stop-idx+2))
        ambiguous_intron = (chrom, start-idx, stop-idx)
        junctions.append(ambiguous_intron)
        idx += 1
    idx = 0
    # while stop-idx+1 < len(R64[chrom]) and R64[chrom][start-idx] == R64[chrom][stop-idx+1]:
    while idx > -10 and stop-idx+1 < R64.get_reference_length(chrom) and R64.fetch(chrom, start-idx, start-idx+1) == R64.fetch(chrom, stop-idx+1, stop-idx+2):
        #print(chrom, start, stop, idx, R64.fetch(chrom, start-idx, start-idx+1) , R64.fetch(chrom, stop-idx+1, stop-idx+2))
        ambiguous_intron = (chrom, start-idx+1, stop-idx+1)
        junctions.append(ambiguous_intron)
        idx -= 1
    potential_5SS = []
    potential_3SS = []
    SS_scores = []
    junctions = sorted(junctions) ## this is critical to ensure that different ambiguous junction inputs yield identical adjusted junctions
    ## in case there are two potential amb junctions with identical scores
    for junction in junctions:
        chrom, amb_start, amb_stop = junction
        if RNA_strand == '+':
            fiveSS = R64.fetch(chrom, amb_start, amb_start+6) # R64[chrom][amb_start:amb_start+6 ]
            threeSS = R64.fetch(chrom, amb_stop-2, amb_stop+1)
            potential_5SS.append( fiveSS )
            potential_3SS.append( threeSS ) # R64[chrom][amb_stop-2:amb_stop+1 ] )
            SS_scores.append( five_SS_score(fiveSS) + three_SS_score(threeSS) )
        if RNA_strand == '-':
            fiveSS = rc( R64.fetch(chrom, amb_stop-5, amb_stop+1) )# R64[chrom][amb_stop-5:amb_stop+1] )
            threeSS = rc( R64.fetch(chrom, amb_start, amb_start+3) )
            potential_5SS.append( fiveSS )
            potential_3SS.append( threeSS )#  R64[chrom][amb_start:amb_start+3 ] ) )
            SS_scores.append( five_SS_score(fiveSS) + three_SS_score(threeSS) )
    
    most_likely_intron = junctions[ SS_scores.index( min(SS_scores) ) ] 
    chrom, amb_start, amb_stop  = most_likely_intron
    intron_coords_adjusted = (amb_start == start)
    if RNA_strand == '+':
        ann_5SS = ((annotated_intron_df['strand'] == RNA_strand) & (annotated_intron_df['start'] == amb_start)).any()
        ann_3SS = ((annotated_intron_df['strand'] == RNA_strand) & (annotated_intron_df['end'] == amb_start)).any()
    else:
        ann_3SS = ((annotated_intron_df['strand'] == RNA_strand) & (annotated_intron_df['start'] == amb_start)).any()
        ann_5SS = ((annotated_intron_df['strand'] == RNA_strand) & (annotated_intron_df['end'] == amb_start)).any()
    five_SS = potential_5SS[ SS_scores.index( min(SS_scores) ) ] 
    three_SS = potential_3SS[ SS_scores.index( min(SS_scores) ) ] 
    canonical_5SS = (five_SS[:2] == 'GT')
    canonical_3SS = (three_SS[1:] == 'AG')
    intron_size = abs(stop - start)
    annotated_junction = (most_likely_intron in annotated_junctions)
    num_amb_junctions = len(potential_5SS)
    ## required later in this script: chrom, amb_start, amb_stop, annotated_junction, 
    ## the rest of the returned values in the list could be generated later in the pipeline, but it is convenient to already process this here 
    return [chrom, amb_start, amb_stop, five_SS, three_SS, RNA_strand, num_amb_junctions, annotated_junction, ann_5SS, ann_3SS, canonical_5SS, canonical_3SS, intron_size, intron_coords_adjusted]

common_intron_info = 'chrom, amb_start, amb_stop, five_SS, three_SS, RNA_strand, num_amb_junctions, annotated_junction, ann_5SS, ann_3SS, canonical_5SS, canonical_3SS, intron_size'.split(', ')
sample_specific_intron_info = 'intron_coords_adjusted, mismatch_near_junction, left_perfect_matches, right_perfect_matches, quality_score_at_left_junction, quality_score_at_right_junction'.split(', ')
adjusted_intron_keys = common_intron_info + sample_specific_intron_info

def process_alignment(read):
    '''
    input: a pysam AlignedSegment object
    '''
    keys = 'flag, chrom, coord, cigar, NH, adjusted_introns, perfect_gapped_alignment, num_introns_found, alignment_score'.split(', ')          
    
    if type(read) != pysam.libcalignedsegment.AlignedSegment or read.cigartuples == None:
        ## values = read.flag, read.reference_name, read.reference_start, read.cigarstring.replace('N', 'D'), read.get_tag('NH'), adjusted_introns, perfect_gapped_alignment, num_introns_found, alignment_score
        values = [None] * 5 + [ [], False, 0, 1000]  ## alignment score is set to 1000 for unmapped reads
        return dict(zip(keys, values))
    
    quality_scores = read.query_qualities
    cigartuples, total_cigar_operation_bp = parse_cigar(read.cigarstring)  
    coordinate = read.reference_start
    chrom = read.reference_name
    RNA_strand = '+' if ( (read.is_read1 and read.is_reverse) or (read.is_read2 and not read.is_reverse) ) else '-'
    ## parsing the cigar has three main functions: 
    ## (1) assess whether alignment has a gap within the allowable range for an intron
    ## (2) calculate the putative intron boundaries (i.e. 5' and 3' splice sites)
    ## (3) gather metrics on the alignment, including the total mismatches, 
    ## lengths of left and right mapped segments around the gap,
    ## and the number of bp matching perfectly on the left and right of the gap
    mapped_segment_length = 0
    total_gapped_bp = 0
    alignment_score = 0
    splice_sites = []
    quality_scores_at_splice_sites = []
    perfect_matches_flanking_splice_sites  = []
    mapped_segment_lengths = []
    perfect_gapped_alignment = True  
    mismatch_near_junction = []
    num_introns_found = 0 
    for cigar_idx in range(len(cigartuples)):  ## cigartuples.append( (mapped_segment_length, operation) )
        bp, operation = cigartuples[cigar_idx]    
        if operation in ('=', 'X', 'D', 'M'):
            mapped_segment_length += bp
            coordinate += bp
        elif operation == 'N':
            total_gapped_bp += bp
            mapped_segment_lengths.append(mapped_segment_length)
            mapped_segment_length = 0
            ## the cigar operation flanking N should always be =, but need to check this just in case
            if cigartuples[cigar_idx-1][1] != '=': 
                left_perfect_matches = 0
            else:
                left_perfect_matches = cigartuples[cigar_idx-1][0] 
            if cigartuples[cigar_idx+1][1] != '=' :
                right_perfect_matches = 0
            else:
                right_perfect_matches = cigartuples[cigar_idx+1][0]            
            if bp >= MINIMUM_INTRON_LENGTH:
                num_introns_found += 1
                splice_sites.append( (coordinate, coordinate + bp - 1) ) # a single read may have more than one splice site
                quality_score_at_left_junction = quality_scores[mapped_segment_length - 1]
                quality_score_at_right_junction = quality_scores[mapped_segment_length]
                quality_scores_at_splice_sites.append( (quality_score_at_left_junction, quality_score_at_right_junction) )
                coordinate = coordinate + bp
                perfect_matches_flanking_splice_sites.append( (left_perfect_matches, right_perfect_matches) )
                mismatch_near_junction.append( (cigar_idx >= 2 and left_perfect_matches < MISMATCH_DIST_FROM_JUNCTION_DISALLOWED) or (cigar_idx+2 < len(cigartuples) and right_perfect_matches < MISMATCH_DIST_FROM_JUNCTION_DISALLOWED) )
            else:
                perfect_gapped_alignment = False
                alignment_score += bp
        if operation in ('X', 'I', 'S', 'D', 'M'):
            perfect_gapped_alignment = False
            alignment_score += bp 
    mapped_segment_lengths.append(mapped_segment_length)
    
    adjusted_introns = []
    for idx in range( num_introns_found ) :
        start, end = splice_sites[idx]
        left_perfect_matches, right_perfect_matches = perfect_matches_flanking_splice_sites[idx]
        quality_score_at_left_junction, quality_score_at_right_junction = quality_scores_at_splice_sites[idx]
        intron = chrom, start, end
       #  if end - start < MAXIMUM_INTRON_LENGTH:
        if intron not in junctions_to_ambiguous_junctions:
            junctions_to_ambiguous_junctions[intron] = adjust_ambiguous_junctions(chrom, start, end, RNA_strand)
            ## adjust_ambiguous_junctions returns: chrom, amb_start, amb_stop, five_SS, three_SS, RNA_strand, num_amb_junctions, annotated_junction, ann_5SS, ann_3SS, canonical_5SS, canonical_3SS, intron_size, intron_coords_adjusted, 
            ## additional info added: mismatch_near_junction, left_perfect_matches, right_perfect_matches
            ## values in adjusted_introns list: 
            ## chrom, amb_start, amb_stop, five_SS, three_SS, RNA_strand, num_amb_junctions, annotated_junction, ann_5SS, ann_3SS, canonical_5SS, canonical_3SS, intron_size, intron_coords_adjusted, mismatch_near_junction, left_perfect_matches, right_perfect_matches
        adjusted_intron_values = junctions_to_ambiguous_junctions[intron] + [mismatch_near_junction[idx], left_perfect_matches, right_perfect_matches, quality_score_at_left_junction, quality_score_at_right_junction ]
        adjusted_introns.append( dict(zip(adjusted_intron_keys, adjusted_intron_values))  )
    values = read.flag, read.reference_name, read.reference_start, read.cigarstring.replace('D', 'N').replace('M', 'X'), read.get_tag('NH'), adjusted_introns, perfect_gapped_alignment, num_introns_found, alignment_score
    return dict(zip(keys, values))  

## open files: COMPASS_bam_outfile, COMPASS_UGA_bam_outfile, COMPASS_ASJ_bam_outfile, COMPASS_USJ_bam_outfile, annotated_gapped_alignments_outfile, unannotated_gapped_alignments_outfile, splice_junction_counts_outfile
## all reads written to COMPASS bam file
COMPASS_bamfilename = OUT_DIR + sample_name + "_COMPASS.bam"
COMPASS_bam_outfile = pysam.AlignmentFile(COMPASS_bamfilename, "wb", template = aligners_to_bam[aligners[0]] )
## write UGA, ASJ and USJ separately
## these are written to enable easily extracting coverage at intron junctions for calculating splicing efficiency      
COMPASS_UGA_bam_filename = OUT_DIR + sample_name + "_COMPASS_UGA.bam"
COMPASS_UGA_bam_outfile = pysam.AlignmentFile(COMPASS_UGA_bam_filename, "wb", template = aligners_to_bam[aligners[0]] )
COMPASS_ASJ_bam_filename = OUT_DIR + sample_name + "_COMPASS_ASJ.bam"
COMPASS_ASJ_bam_outfile = pysam.AlignmentFile(COMPASS_ASJ_bam_filename, "wb", template = aligners_to_bam[aligners[0]] )
COMPASS_USJ_bam_filename = OUT_DIR + sample_name + "_COMPASS_USJ.bam"
COMPASS_USJ_bam_outfile = pysam.AlignmentFile(COMPASS_USJ_bam_filename, "wb", template = aligners_to_bam[aligners[0]] )

COMPASS_bam_sorted_filename = OUT_DIR + sample_name + "_COMPASS_sorted.bam"
COMPASS_UGA_bam_sorted_filename = OUT_DIR + sample_name + "_COMPASS_UGA_sorted.bam"
COMPASS_USJ_bam_sorted_filename = OUT_DIR + sample_name + "_COMPASS_USJ_sorted.bam"
COMPASS_ASJ_bam_sorted_filename = OUT_DIR + sample_name + "_COMPASS_ASJ_sorted.bam"

annotated_gapped_alignments_outfilename = OUT_DIR + sample_name + "_annotated_gapped.tsv"
annotated_gapped_alignments_outfile = open(annotated_gapped_alignments_outfilename, 'w')
## unnanoted junctions are written separately from annotated junctions, 
## as the annotated junction files can be too large to be loaded into memory in a data table
unannotated_gapped_alignments_outfilename = OUT_DIR + sample_name + "_unannotated_gapped.tsv"
unannotated_gapped_alignments_outfile = open(unannotated_gapped_alignments_outfilename, 'w')   
## if an intron is found, info will be written to either annotated or unannotated outfiles, 
## depending on whether the alignment with the best score maps to an annotated or unannotated intron, respectively

JUNCTION_MIN_SCORE_DISAGREEMENT_OUTFILENAME = OUT_DIR + sample_name + "_min_score_alignments_disagree_on_junction.tsv"
JUNCTION_MIN_SCORE_DISAGREEMENT_OUTFILE = open(JUNCTION_MIN_SCORE_DISAGREEMENT_OUTFILENAME, 'w')   

intron_found_to_alignments_outfiles = {True: annotated_gapped_alignments_outfile, False: unannotated_gapped_alignments_outfile}

splice_junction_counts_outfilename = OUT_DIR + sample_name + "_COMPASS_splice_junctions.tsv"
splice_junction_counts_outfile = open(splice_junction_counts_outfilename, 'w')

read_num = 0
aligner_to_current_aligned_segment = {}
aligner_to_current_read1_alignment = {}
aligner_to_current_read2_alignment = {}

## keep a dictionary for junctions to aligners to counts
junction_counts = {}
## keep a dictionary for junctions to count intron_coords_adjusted, count mismatch_near_junction, list left_perfect_matches, and list right_perfect_matches
COMPASS_junction_statistics = {}

concordant_ungapped_counts = 0
concordant_gapped_counts = 0
discordant_ungapped_counts = 0
discordant_gapped_counts = 0
min_score_has_no_intron_cigars_equal = {True: 0, False: 0}
min_score_has_annotated_intron_cigars_equal = {True: 0, False: 0}
min_score_has_annotated_intron_junctions_equal = {True: 0, False: 0}
min_score_has_unannotated_intron_cigars_equal = {True: 0, False: 0}
min_score_has_unannotated_intron_junctions_equal = {True: 0, False: 0}

columns = ['read_pair_ID', 'aligner_has_min_score', 'min_score_has_annotated_intron', 'min_score_has_no_intron', 'aligner', 'score',  'num_introns_found', 'perfect_gapped_alignment', 'chrom', 'amb_start', 'amb_stop', 'five_SS', 'three_SS', 'RNA_strand', 'annotated_junction', 'ann_5SS', 'ann_3SS', 'canonical_5SS', 'canonical_3SS', 'intron_size', 'intron_coords_adjusted', 'mismatch_near_junction', 'left_perfect_matches', 'right_perfect_matches', 'read1_flag', 'read1_chrom', 'read1_coord', 'read1_cigar', 'read1_NH', 'read2_flag', 'read2_chrom', 'read2_coord', 'read2_cigar', 'read2_NH' ]
header = '\t'.join(columns) + '\n'
[e.write(header) for e in intron_found_to_alignments_outfiles.values()]
eof = False

while read_num <= reads_to_process and not eof:
    if read_num % READ_NUM_PROGRESS == 0:
        print('COMPASS currently at read number:', read_num)
        print('concordant_ungapped_counts:', concordant_ungapped_counts)
        print('discordant_ungapped_counts:', discordant_ungapped_counts)
        print('concordant_gapped_counts:', concordant_gapped_counts)
        print('discordant_gapped_counts:', discordant_gapped_counts)
        print('equivalent alignments for all min score aligners when a min score has ungapped alignment, ', min_score_has_no_intron_cigars_equal)
        print('equivalent alignments for all min score aligners when a min score has annotated junction, ', min_score_has_annotated_intron_cigars_equal)
        print('equivalent junctions for all min score aligners when a min score has annotated junction, ', min_score_has_annotated_intron_junctions_equal)
        print('equivalent alignments for all min score aligners when a min score has unannotated junction, ', min_score_has_unannotated_intron_cigars_equal)
        print('equivalent junctions found for all min score aligners when a min score has unannotated junction, ', min_score_has_unannotated_intron_junctions_equal, '\n')
    read_num += 1
    gapped_alignment_present = False

    ## for each aligner, gather the alignments for the current read number, if present, otherwise None will be used
    for aligner in aligners: 
        secondary_alignment = False
        aligner_to_current_read1_alignment[aligner] = None
        aligner_to_current_read2_alignment[aligner] = None
     ## if alignment file reaches the end during iteration, that indicates all reads have been processed
        while aligner_to_current_read_num[aligner] < read_num or secondary_alignment: 
            try:
                aligned_segment = next( aligners_to_bam[aligner] )
                secondary_alignment = aligned_segment.is_secondary
                ## only consider primary alignments from each aligner
            except StopIteration:
                eof = True
                break
                break
                break
     
            aligner_bam_read_num = int(aligned_segment.query_name.split('_')[0])
            if aligner_bam_read_num < aligner_to_current_read_num[aligner] :
                print('lower read number', aligner_bam_read_num, 'found after higher read num', aligner_to_current_read_num[aligner] , 'for aligner:', aligner)
                raise ValueError
            if not secondary_alignment:
                aligner_to_current_read_num[aligner] = aligner_bam_read_num
                aligner_to_current_aligned_segment[aligner] = aligned_segment

        if aligner_to_current_read_num[aligner] == read_num: ## process only the current read number, 
        ## as there is a chance during the last while loop that an aligner didn't report a mapping for a read
        ## and it would therefore jump to the next read, which would need to be kept in memory for the next round of COMPASS
        ## after getting the first read in a pair for a given read num, need to get the other read from the pair,
        ## and possibly need to go through several other secondary mappings before getting to the primary mapping of the other read
            aligned_segment = aligner_to_current_aligned_segment[aligner]
            read1_first = aligned_segment.is_read1
            if read1_first:
                aligner_to_current_read1_alignment[aligner] = aligned_segment 
            else:
                aligner_to_current_read2_alignment[aligner] = aligned_segment
                
            while aligner_to_current_read_num[aligner] == read_num and (aligner_to_current_read1_alignment[aligner] == None or aligner_to_current_read2_alignment[aligner] == None):    
                try:
                    aligned_segment = next( aligners_to_bam[aligner] )
                except StopIteration:
                    eof = True
                    break
                    break
                    break
                aligner_to_current_read_num[aligner] = int(aligned_segment.query_name.split('_')[0])
                aligner_to_current_aligned_segment[aligner] = aligned_segment
        ## this code works under the assumption for a given read pair, there are no more than two primary alignments,
        ## one for each read in the pair.  This code allows for the possibility that one read of a pair has an unreported alignment
        ## in which case that read will be assigned a None value

                if aligner_to_current_read_num[aligner] == read_num and not aligned_segment.is_secondary:
                    if read1_first and aligned_segment.is_read1:
                        print('two adjacent read 1s') ## this should never print
                    if aligned_segment.is_read1:
                        aligner_to_current_read1_alignment[aligner] = aligned_segment 
                    else:
                        aligner_to_current_read2_alignment[aligner] = aligned_segment

    ## initialize lists for info on each alignment, access later by idx using this construct: for idx in range(len(aligners))
    ## first process alignment scores and confirm if best score meets threshold for further processing
    alignment_scores = []
    num_intron_found_lst = []
    R1_lst = []
    R2_lst = []
    for aligner in aligners:
        R1 = process_alignment( aligner_to_current_read1_alignment[aligner] )
        R2 = process_alignment( aligner_to_current_read2_alignment[aligner] )
        R1_lst.append(R1)
        R2_lst.append(R2)
        alignment_scores.append( R1['alignment_score'] + R2['alignment_score'] )
        num_intron_found_lst.append( max( R1['num_introns_found'], R2['num_introns_found'] ))
    best_score = min(alignment_scores)

    if best_score < MAX_EDIT_DIST_TO_CONSIDER:
        min_score_has_annotated_intron = False  ## does any best alignment have annotated intron
        min_score_has_no_intron = False  ## does any best alignment have no intron
        perfect_lst = []
        read1_alignments = []
        read2_alignments = []
        splice_junction_lst = []
        idx_best_score_with_annotated_intron = []
        idx_best_score_with_no_intron = []
        idx_with_best_score = []
        aligners_with_best_score = []
        intron_lengths_lst = []
        idx_without_mismatch_near_junction = []
        comment_tag_lst = []
        alignment_type_lst = []
        
## need to loop through the aligners first to pre-compute:  
## min_score_has_annotated_intron, min_score_has_no_intron, num_intron_found_lst, perfect_lst
        for idx in range(len(aligners)):
            aligner = aligners[idx]
            
            R1 = R1_lst[idx]
            R2 = R2_lst[idx]
            
            annotated_junction = False
            mismatch_near_junction = False
            splice_junctions_in_read_pair = []
            
            intron_lengths = [0]
            for e in R1, R2:
                if e['adjusted_introns'] != []:
                ## values in adjusted_introns list: chrom, amb_start, amb_stop, five_SS, three_SS, RNA_strand, num_amb_junctions, annotated_junction, ann_5SS, ann_3SS, canonical_5SS, canonical_3SS, intron_size, intron_coords_adjusted, mismatch_near_junction, left_perfect_matches, right_perfect_matches, quality_score_at_left_junction, quality_score_at_right_junction
                # chrom, amb_start, amb_stop, five_SS, three_SS, RNA_strand, num_amb_junctions, annotated_junction, ann_5SS, ann_3SS, canonical_5SS, canonical_3SS = adjusted_intron  
                ## read 1 and read 2 info entered on separate lines if both contain a gapped alignment, 
                ## later on these can be further processed and integrated by groupby operations
                ## if a read contains more than intron (quite rare), it is entered on multiple lines
                    for intron in e['adjusted_introns']:
                        splice_junctions_in_read_pair.append( tuple( [intron[e] for e in  'chrom, amb_start, amb_stop, five_SS, three_SS, RNA_strand'.split(', ')]  ) )
                        if intron['annotated_junction']:
                            annotated_junction = True
                        if intron['mismatch_near_junction']:
                            mismatch_near_junction = True
                        if alignment_scores[idx] == best_score and annotated_junction:
                            min_score_has_annotated_intron = True ## does any best alignment have annotated intron
                            idx_best_score_with_annotated_intron.append(idx)
                        if not mismatch_near_junction:
                            idx_without_mismatch_near_junction.append(idx)  
                        intron_lengths.append( intron['intron_size'] )
            intron_lengths_lst.append( max(intron_lengths) )
            splice_junction_lst.append( tuple(splice_junctions_in_read_pair) )
            if alignment_scores[idx] == best_score:
                idx_with_best_score.append(idx)
                aligners_with_best_score.append(aligners[idx])
                if num_intron_found_lst[idx] == 0:
                    min_score_has_no_intron = True ## is any best alignment ungapped
                    idx_best_score_with_no_intron.append(idx)
            ## need to consider what happens when an intron greater than max allowed length is present    
            ## if introns are present, is there a perfect gapped alignment present in the read pair?
            
            if R1['num_introns_found'] > 0 and R2['num_introns_found'] > 0:
                perfect_gapped_alignment = R1['perfect_gapped_alignment'] or R2['perfect_gapped_alignment']
            elif R1['num_introns_found'] > 0:
                perfect_gapped_alignment = R1['perfect_gapped_alignment']
            elif R2['num_introns_found'] > 0:
                perfect_gapped_alignment = R2['perfect_gapped_alignment']
            else:
                perfect_gapped_alignment = False
            perfect_lst.append( perfect_gapped_alignment )  
                
            if (R1['num_introns_found'] + R2['num_introns_found']) > 0:
                if annotated_junction:
                    alignment_type = 'ASJ'
                else:
                    alignment_type = 'USJ'
            else:
                alignment_type = 'UGA'
            alignment_type_lst.append(alignment_type)
            
            comment_tag = [aligner, R1['chrom'], R1['coord'], R1['cigar'], R2['chrom'], R2['coord'], R2['cigar'], alignment_scores[idx], alignment_type_lst[idx] ]
            comment_tag_lst.append( ','.join([str(e) for e in comment_tag])  )
              
        ## write the complete COMPASS table for all read alignments (big table = num reads X num aligners)  
        for idx in range(len(aligners)):            
            read_pair_info = [read_num, alignment_scores[idx] == best_score, min_score_has_annotated_intron, min_score_has_no_intron ] + [ f[idx] for f in (aligners, alignment_scores, num_intron_found_lst, perfect_lst) ]     
            R1 = R1_lst[idx]
            R2 = R2_lst[idx]  
            read1_info = [R1[e] for e in keys[:5]] 
            read2_info = [R2[e] for e in keys[:5]] 
            read1_alignments.append(tuple(read1_info))
            read2_alignments.append(tuple(read2_info))
            for e in R1, R2:
                if e['adjusted_introns'] != []:
                # chrom, amb_start, amb_stop, five_SS, three_SS, RNA_strand, num_amb_junctions, annotated_junction, ann_5SS, ann_3SS, canonical_5SS, canonical_3SS = adjusted_intron  
                ## read 1 and read 2 info entered on separate lines if both contain a gapped alignment, 
                ## later on these can be further processed and integrated by groupby operations
                ## if a read contains more than intron (quite rare), it is entered on multiple lines
                    for intron in e['adjusted_introns']:
                        # # # adjusted_introns --> chrom, amb_start, amb_stop, five_SS, three_SS, RNA_strand, num_amb_junctions, annotated_junction, ann_5SS, ann_3SS, canonical_5SS, canonical_3SS, intron_size, intron_coords_adjusted, mismatch_near_junction , left_perfect_matches, right_perfect_matches, quality_score_at_left_junction, quality_score_at_right_junction
                        output = read_pair_info + list(intron.values()) + read1_info + read2_info
                        output = [str(e) for e in output]
                        if max(num_intron_found_lst) > 0:
                            intron_found_to_alignments_outfiles[ min_score_has_annotated_intron ].write('\t'.join(output) + '\n')           
            if R1['adjusted_introns'] == R2['adjusted_introns'] == []: ## this line will only execute if the preceding for loop does nothing
                ## if no gapped alignments found, enter 13 NA values
                output = read_pair_info + ['NA']*15 + read1_info + read2_info
                output = [str(e) for e in output]
                if max(num_intron_found_lst) > 0:
                    intron_found_to_alignments_outfiles[ min_score_has_annotated_intron ].write('\t'.join(output) + '\n')
                    ## ungapped alignments are only written to the gapped outfiles if one aligner had a gapped alignment
                    ## need to consider what happens when an intron greater than max allowed length is present   
                    ## it will be written to this outfile as NA's
        
        ## write the optimal alignments to bam files   
        all_aligners_agree_on_cigar = False
        if len(set(read1_alignments)) == len(set(read2_alignments)) == 1:  ## check if all alignments are identical based on: flag, chrom, coord, cigar, NH
            all_aligners_agree_on_cigar = True
            if min(num_intron_found_lst) > 0:
                concordant_gapped_counts += 1
            else:
                concordant_ungapped_counts += 1
        elif max(num_intron_found_lst) > 0:
            discordant_gapped_counts += 1  ## at least one aligner has an intron, but cigars don't agree
        else:
            discordant_ungapped_counts += 1 ## no aligner has an intron, and cigars don't agree
            
        idx_with_best_score = [index for index, value in enumerate(alignment_scores) if value == min(alignment_scores)]
        num_aligners_with_best_score = len(idx_with_best_score)
        min_score_unannotated_intron_junctions_disagree = False
        if all_aligners_agree_on_cigar:
            aligner_idx_used = 0
        elif min_score_has_no_intron:                     ## need to consider what happens when an intron greater than max allowed length is present   
            aligner_idx_used = random.choice(idx_best_score_with_no_intron)
            min_score_has_no_intron_cigars_equal[ len(set([read1_alignments[i] for i in idx_best_score_with_no_intron])) == len(set([read2_alignments[i] for i in idx_best_score_with_no_intron])) == 1] += 1
        elif min_score_has_annotated_intron:
            aligner_idx_used = random.choice(idx_best_score_with_annotated_intron)
            min_score_has_annotated_intron_cigars_equal[ len(set([read1_alignments[i] for i in idx_best_score_with_annotated_intron])) == len(set([read2_alignments[i] for i in idx_best_score_with_annotated_intron])) == 1] += 1
            min_score_has_annotated_intron_junctions_equal[ len(set([splice_junction_lst[i] for i in idx_best_score_with_annotated_intron])) == 1] += 1  
            # if not (len(set([read1_alignments[i] for i in idx_best_score_with_annotated_intron])) == len(set([read2_alignments[i] for i in idx_best_score_with_annotated_intron])) == 1):
            #     print(alignment_type, alignment_scores)
            #     print(idx_with_best_score)
            #     print(idx_best_score_with_annotated_intron)
            #     print(read1_alignments, read2_alignments)
            #     for read_1_alignment_test in read1_alignments:
            #         print(read_1_alignment_test)
            #     for read_2_alignment_test in read2_alignments:
            #         print(read_2_alignment_test)
            #     for SJ in splice_junction_lst:
            #         print(SJ)
            #     print('\n')
            # if (99, 'chrVII', 364736, '99=1X8=1X7=1X32=', 1) not in read1_alignments and (83, 'chrXIII', 226330, '150=', 1) not in read1_alignments and (99, 'chrIX', 98721, '150=', 1) not in read1_alignments:
            #     raise ValueError
        else:
            min_score_has_unannotated_intron_cigars_equal[ len(set([read1_alignments[i] for i in idx_with_best_score])) == len(set([read2_alignments[i] for i in idx_with_best_score])) == 1] += 1
            # if len(set([read1_alignments[i] for i in idx_with_best_score])) == len(set([read2_alignments[i] for i in idx_with_best_score])) != 1:
            #     print('read1_alignments', read1_alignments, idx_with_best_score, '\n' )
            min_score_has_unannotated_intron_junctions_equal[ len(set([splice_junction_lst[i] for i in idx_with_best_score])) == 1] += 1
            if len(set([splice_junction_lst[i] for i in idx_with_best_score])) != 1:
                min_score_unannotated_intron_junctions_disagree = True
                JUNCTION_MIN_SCORE_DISAGREEMENT_OUTFILE.write(str(aligners) + '\n' + str(read1_alignments) + '\n' + str(read2_alignments) + '\n' + str(splice_junction_lst) + '\n' + str(alignment_scores) + '\n'*2)
            min_intron_length_with_best_score = min([intron_lengths_lst[i] for i in idx_with_best_score ])
            idx_with_best_score_and_min_intron_size = [idx for idx in range(len(intron_lengths_lst)) if intron_lengths_lst[idx] == min_intron_length_with_best_score and idx in idx_with_best_score]
            tie_breaker_lst = [i for i in idx_without_mismatch_near_junction if i in idx_with_best_score_and_min_intron_size]
            if tie_breaker_lst != []:
                aligner_idx_used = random.choice(tie_breaker_lst)
            else:
                aligner_idx_used = random.choice(idx_with_best_score_and_min_intron_size)
        aligner_used = aligners[aligner_idx_used]    
        alignment_type = alignment_type_lst[aligner_idx_used]
        comment_tag_out = ';'.join(comment_tag_lst) + '/' + str(num_aligners_with_best_score) + ',' + alignment_type
        for read_alignment in aligner_to_current_read1_alignment, aligner_to_current_read2_alignment:
            # if read_alignment != None:
            ## None type error here???
            read_alignment[aligner_used].set_tag('PG', ','.join(aligners_with_best_score), 'Z')
            read_alignment[aligner_used].set_tag('CO', comment_tag_out, 'Z')
            
        [ COMPASS_bam_outfile.write(e[aligner_used]) for e in (aligner_to_current_read1_alignment, aligner_to_current_read2_alignment) if e != None ]
        if alignment_type == 'UGA':
            [COMPASS_UGA_bam_outfile.write(e[aligner_used]) for e in (aligner_to_current_read1_alignment, aligner_to_current_read2_alignment) if e != None]
        if alignment_type == 'USJ':
            [ COMPASS_USJ_bam_outfile.write(e[aligner_used]) for e in (aligner_to_current_read1_alignment, aligner_to_current_read2_alignment) if e != None]
        if alignment_type == 'ASJ':
            [ COMPASS_ASJ_bam_outfile.write(e[aligner_used]) for e in (aligner_to_current_read1_alignment, aligner_to_current_read2_alignment) if e != None]

        ## write the optimal alignment to the COMPASS splice junction counts file  
        if alignment_type == 'USJ' or alignment_type == 'ASJ':
            R1_selected = R1_lst[aligner_idx_used]
            R2_selected = R2_lst[aligner_idx_used]
            perfect_gapped_alignment = perfect_lst[aligner_idx_used]
            num_introns_found = num_intron_found_lst[aligner_idx_used]
            
            introns_to_attributes = {}
            for e in R1_selected, R2_selected:
                # # # adjusted_introns --> chrom, amb_start, amb_stop, five_SS, three_SS, RNA_strand, num_amb_junctions, annotated_junction, ann_5SS, ann_3SS, canonical_5SS, canonical_3SS, intron_size, intron_coords_adjusted, mismatch_near_junction, left_perfect_matches, right_perfect_matches, quality_score_at_left_junction, quality_score_at_right_junction
                # chrom, amb_start, amb_stop, five_SS, three_SS, RNA_strand, num_amb_junctions, annotated_junction, ann_5SS, ann_3SS, canonical_5SS, canonical_3SS, intron_size,
                # first 10 are invariate attributes of the coords, last 4 could be different (last 2 almost certainly different)
                # intron_coords_adjusted unlikely to be different for read pairs overlappign the same junction
                # intron_coords_adjusted, mismatch_near_junction, left_perfect_matches, right_perfect_matches, quality_score_at_left_junction, quality_score_at_right_junction
                for intron in e['adjusted_introns']:
                    chrom_coords = tuple(intron[e] for e in common_intron_info)
                    if chrom_coords not in introns_to_attributes:
                        introns_to_attributes[chrom_coords] = {k:v for (k,v) in intron.items() if k in sample_specific_intron_info}
                        # intron_coords_adjusted, mismatch_near_junction, left_perfect_matches, right_perfect_matches, quality_score_at_left_junction, quality_score_at_right_junction
                    ## if read pairs overlap on same junction, integrate the data
                    else:
                        introns_to_attributes[chrom_coords]['intron_coords_adjusted'] = intron['intron_coords_adjusted'] and introns_to_attributes[chrom_coords]['intron_coords_adjusted']
                        introns_to_attributes[chrom_coords]['mismatch_near_junction'] = intron['mismatch_near_junction'] and introns_to_attributes[chrom_coords]['mismatch_near_junction']
                        introns_to_attributes[chrom_coords]['left_perfect_matches'] = max(intron['left_perfect_matches'],introns_to_attributes[chrom_coords]['left_perfect_matches'])
                        introns_to_attributes[chrom_coords]['right_perfect_matches'] = max(intron['right_perfect_matches'], introns_to_attributes[chrom_coords]['right_perfect_matches'])
                        introns_to_attributes[chrom_coords]['quality_score_at_left_junction'] = max(intron['quality_score_at_left_junction'],introns_to_attributes[chrom_coords]['quality_score_at_left_junction'])
                        introns_to_attributes[chrom_coords]['quality_score_at_right_junction'] = max(intron['quality_score_at_right_junction'], introns_to_attributes[chrom_coords]['quality_score_at_right_junction'])
            for chrom_coords in introns_to_attributes:    
                ## update COMPASS statistics dict
                #  values = read.flag, read.reference_name, read.reference_start, read.cigarstring.replace('D', 'N').replace('M', 'X'), read.get_tag('NH'), adjusted_introns, perfect_gapped_alignment, num_introns_found, alignment_score
                if chrom_coords not in COMPASS_junction_statistics:
                    COMPASS_junction_statistics[chrom_coords] = {}
                    COMPASS_junction_statistics[chrom_coords]['intron_coords_adjusted_counts'] = 0
                    COMPASS_junction_statistics[chrom_coords]['mismatch_near_junction_counts'] = 0
                    COMPASS_junction_statistics[chrom_coords]['perfect_gapped_alignment_counts'] = 0
                    COMPASS_junction_statistics[chrom_coords]['left_perfect_matches'] = []
                    COMPASS_junction_statistics[chrom_coords]['right_perfect_matches'] = []
                    COMPASS_junction_statistics[chrom_coords]['quality_score_at_left_junction'] = []
                    COMPASS_junction_statistics[chrom_coords]['quality_score_at_right_junction'] = []
                    COMPASS_junction_statistics[chrom_coords]['alignment_scores'] = []
                    COMPASS_junction_statistics[chrom_coords]['num_introns_found'] = []                   
                    COMPASS_junction_statistics[chrom_coords]['read_1_coord_cigars'] = []
                    COMPASS_junction_statistics[chrom_coords]['read_2_coord_cigars'] = []
                    COMPASS_junction_statistics[chrom_coords]['read_pair_cigar_combos'] = set([])
                    COMPASS_junction_statistics[chrom_coords]['min_score_unannotated_intron_junctions_disagree'] = 0
                if introns_to_attributes[chrom_coords]['intron_coords_adjusted']:
                    COMPASS_junction_statistics[chrom_coords]['intron_coords_adjusted_counts'] += 1
                if introns_to_attributes[chrom_coords]['mismatch_near_junction']:
                    COMPASS_junction_statistics[chrom_coords]['mismatch_near_junction_counts'] += 1
                if perfect_gapped_alignment:
                    COMPASS_junction_statistics[chrom_coords]['perfect_gapped_alignment_counts'] += 1
                COMPASS_junction_statistics[chrom_coords]['left_perfect_matches'].append(introns_to_attributes[chrom_coords]['left_perfect_matches'])
                COMPASS_junction_statistics[chrom_coords]['right_perfect_matches'].append(introns_to_attributes[chrom_coords]['right_perfect_matches'])
                COMPASS_junction_statistics[chrom_coords]['quality_score_at_left_junction'].append(introns_to_attributes[chrom_coords]['quality_score_at_left_junction'])
                COMPASS_junction_statistics[chrom_coords]['quality_score_at_right_junction'].append(introns_to_attributes[chrom_coords]['quality_score_at_right_junction'])
                COMPASS_junction_statistics[chrom_coords]['alignment_scores'].append(best_score)
                COMPASS_junction_statistics[chrom_coords]['num_introns_found'].append(num_introns_found)
                COMPASS_junction_statistics[chrom_coords]['read_1_coord_cigars'].append( ( R1_selected['coord'], R1_selected['cigar'] ) )
                COMPASS_junction_statistics[chrom_coords]['read_2_coord_cigars'].append( ( R2_selected['coord'], R2_selected['cigar'] ) )
                COMPASS_junction_statistics[chrom_coords]['read_pair_cigar_combos'].add( (R1_selected['coord'], R1_selected['cigar'], R2_selected['coord'], R2_selected['cigar']) )
                if min_score_unannotated_intron_junctions_disagree:
                    COMPASS_junction_statistics[chrom_coords]['min_score_unannotated_intron_junctions_disagree'] += 1
                ## update junction counts
                if chrom_coords not in junction_counts:
                    junction_counts[chrom_coords] = {}
                    for aligner in aligners:
                        junction_counts[chrom_coords][aligner] = 0
                junction_counts[chrom_coords]['COMPASS'] = junction_counts[chrom_coords].get('COMPASS', 0) + 1
                chrom, start, end = chrom_coords[:3]
                for idx in range(len(aligners)):
                    aligner = aligners[idx]
                    R1 = R1_lst[idx]
                    R2 = R2_lst[idx]
                    intron_found_by_this_aligner = False
                    for e in R1, R2:
                        for intron in e['adjusted_introns']:
                            chrom_coords_to_check = tuple([intron[e] for e in common_intron_info])
                            if chrom_coords_to_check == chrom_coords:
                                intron_found_by_this_aligner = True
                            if chrom == 'chrX' and start > 73800 and end < 74300:
                                print(aligner, intron)
                    if intron_found_by_this_aligner:
                        junction_counts[chrom_coords][aligner] += 1
                        
[aligners_to_bam[aligner].close() for aligner in aligners]

[e.close() for e in (COMPASS_bam_outfile, COMPASS_UGA_bam_outfile, COMPASS_ASJ_bam_outfile, COMPASS_USJ_bam_outfile, annotated_gapped_alignments_outfile, unannotated_gapped_alignments_outfile, JUNCTION_MIN_SCORE_DISAGREEMENT_OUTFILE) ] 

## pysam sort bam

pysam.sort("-o", COMPASS_UGA_bam_sorted_filename, COMPASS_UGA_bam_filename)
pysam.sort("-o", COMPASS_USJ_bam_sorted_filename, COMPASS_USJ_bam_filename)
pysam.sort("-o", COMPASS_ASJ_bam_sorted_filename, COMPASS_ASJ_bam_filename)
pysam.sort("-o", COMPASS_bam_sorted_filename, COMPASS_bamfilename)

for sorted_bam in (COMPASS_UGA_bam_sorted_filename, COMPASS_USJ_bam_sorted_filename, COMPASS_ASJ_bam_sorted_filename):
    pysam.index(sorted_bam)

COMPASS_UGA_sorted_bam = pysam.AlignmentFile(COMPASS_UGA_bam_sorted_filename, "rb")
       
##  use pybedtools to get strand specific coverage at 5 and 3 SS       

# chrom, amb_start, amb_stop, five_SS, three_SS, RNA_strand, num_amb_junctions, annotated_junction, ann_5SS, ann_3SS, canonical_5SS, canonical_3SS, intron_size, intron_coords_adjusted, mismatch_near_junction, left_perfect_matches, right_perfect_matches, quality_score_at_left_junction, quality_score_at_right_junction                        
## columns are composed of 3 components: sample invariate info on the junction, aggregate COMPASS statistics across all junctions, and junction counts

sample_info_header = ['sample_name', 'total_reads', 'fiveSS_unspliced_reads', 'threeSS_unspliced_reads']
junction_header = 'chrom, amb_start, amb_stop, five_SS, three_SS, RNA_strand, num_amb_junctions, annotated_junction, ann_5SS, ann_3SS, canonical_5SS, canonical_3SS, intron_size'.split(', ')
counts_header = ['COMPASS_counts'] + [e + '_counts' for e in aligners]

stats_header_part_1 = 'min_score_unannotated_intron_junctions_disagree, intron_coords_adjusted_counts, mismatch_near_junction_counts, perfect_gapped_alignment_counts'.split(', ')
stats_header_part_2 = 'US_perfect_matches, DS_perfect_matches, five_SS_Q_scores, three_SS_Q_scores, alignment_scores, num_introns_found, read_1_coord_cigars, read_2_coord_cigars'.split(', ')
stats_header_part_3 = 'num_unique_read_pair_cigars, max_US_perfect_matches, max_DS_perfect_matches, unique_US_perfect_matches, unique_DS_perfect_matches, best_alignment_score, mean_alignment_score, median_alignment_score'.split(', ')
stats_header_part_4 = 'mean_five_SS_Q_score, mean_three_SS_Q_score'.split(', ')
COMPASS_header = '\t'.join( sample_info_header + junction_header + counts_header + stats_header_part_1 + stats_header_part_2 + stats_header_part_3 + stats_header_part_4) + '\n'
splice_junction_counts_outfile.write(COMPASS_header)

for chrom_coords in sorted(COMPASS_junction_statistics.keys()):
 
    chrom, amb_start, amb_stop, five_SS, three_SS, RNA_strand = chrom_coords[:6]
    fiveSS_unspliced_reads, threeSS_unspliced_reads = 0, 0
    if RNA_strand == '+':
        US_perfect_matches = COMPASS_junction_statistics[chrom_coords]['left_perfect_matches']
        DS_perfect_matches = COMPASS_junction_statistics[chrom_coords]['right_perfect_matches']        
        five_SS_Q_scores = COMPASS_junction_statistics[chrom_coords]['quality_score_at_left_junction']
        three_SS_Q_scores = COMPASS_junction_statistics[chrom_coords]['quality_score_at_right_junction']
        for pileupcolumn in COMPASS_UGA_sorted_bam.pileup(chrom, amb_start, amb_start + 1, truncate = True):
            fiveSS_unspliced_reads = pileupcolumn.n  ## this should be modified to obtain strand specific unspliced reads
        for pileupcolumn in COMPASS_UGA_sorted_bam.pileup(chrom, amb_stop - 1, amb_stop, truncate = True):
            threeSS_unspliced_reads = pileupcolumn.n        
    else:
        DS_perfect_matches = COMPASS_junction_statistics[chrom_coords]['left_perfect_matches']
        US_perfect_matches = COMPASS_junction_statistics[chrom_coords]['right_perfect_matches']
        five_SS_Q_scores = COMPASS_junction_statistics[chrom_coords]['quality_score_at_right_junction']
        three_SS_Q_scores = COMPASS_junction_statistics[chrom_coords]['quality_score_at_left_junction']
        for pileupcolumn in COMPASS_UGA_sorted_bam.pileup(chrom, amb_start, amb_start + 1, truncate = True):
             threeSS_unspliced_reads = pileupcolumn.n
        for pileupcolumn in COMPASS_UGA_sorted_bam.pileup(chrom, amb_stop - 1, amb_stop, truncate = True):
             fiveSS_unspliced_reads = pileupcolumn.n
    
    sample_info_output = [ sample_name, str(read_num), str(fiveSS_unspliced_reads), str(threeSS_unspliced_reads)]
    junction_output = [str(e) for e in chrom_coords]
    counts_output = [ str(junction_counts[chrom_coords]['COMPASS']) ] + [ str(junction_counts[chrom_coords][e]) for e in aligners ]

    stats_output_part_1 = [ str(COMPASS_junction_statistics[chrom_coords][e]) for e in stats_header_part_1]
    
    COMPASS_junction_statistics[chrom_coords]['US_perfect_matches'] = US_perfect_matches
    COMPASS_junction_statistics[chrom_coords]['DS_perfect_matches'] = DS_perfect_matches
    COMPASS_junction_statistics[chrom_coords]['five_SS_Q_scores'] = five_SS_Q_scores
    COMPASS_junction_statistics[chrom_coords]['three_SS_Q_scores'] = three_SS_Q_scores

    stats_output_part_2 = [ str(count_frequency(COMPASS_junction_statistics[chrom_coords][e], NUM_ITEMS_TO_REPORT)) for e in stats_header_part_2]    
    
    US_perfect_matches_set = set( US_perfect_matches )
    DS_perfect_matches_set = set( DS_perfect_matches )
    stats_output_part_3_raw = [len(COMPASS_junction_statistics[chrom_coords]['read_pair_cigar_combos']), max(US_perfect_matches_set), max(DS_perfect_matches_set), len(US_perfect_matches_set), len(DS_perfect_matches_set), min(COMPASS_junction_statistics[chrom_coords]['alignment_scores']), statistics.mean(COMPASS_junction_statistics[chrom_coords]['alignment_scores']), statistics.median(COMPASS_junction_statistics[chrom_coords]['alignment_scores'])  ]
    stats_output_part_3 = [str(e) for e in stats_output_part_3_raw]
    
    stats_output_part_4_raw = [ statistics.mean(five_SS_Q_scores), statistics.mean(three_SS_Q_scores)  ]
    stats_output_part_4 = [str(e) for e in stats_output_part_4_raw]
    
    COMPASS_output = '\t'.join(sample_info_output + junction_output + counts_output + stats_output_part_1 + stats_output_part_2 + stats_output_part_3 + stats_output_part_4)  + '\n'
    splice_junction_counts_outfile.write(COMPASS_output)

splice_junction_counts_outfile.close()

print('COMPASS finished at read number:', read_num)
print('concordant_ungapped_counts:', concordant_ungapped_counts)
print('discordant_ungapped_counts:', discordant_ungapped_counts)
print('concordant_gapped_counts:', concordant_gapped_counts)
print('discordant_gapped_counts:', discordant_gapped_counts)
print('equivalent alignments for all min score aligners when a min score has ungapped alignment, ', min_score_has_no_intron_cigars_equal)
print('equivalent alignments for all min score aligners when a min score has annotated junction, ', min_score_has_annotated_intron_cigars_equal)
print('equivalent junctions for all min score aligners when a min score has annotated junction, ', min_score_has_annotated_intron_junctions_equal)
print('equivalent alignments for all min score aligners when a min score has unannotated junction, ', min_score_has_unannotated_intron_cigars_equal)
print('equivalent junctions found for all min score aligners when a min score has unannotated junction, ', min_score_has_unannotated_intron_junctions_equal, '\n')

print('COMPASS program finished successfully')

