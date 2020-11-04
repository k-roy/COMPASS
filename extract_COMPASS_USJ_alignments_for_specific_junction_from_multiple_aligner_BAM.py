#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 13 20:24:55 2020

@author: kevinroy

goes back to each individual aligners name_sorted bam file 
grabs reads which map to unannotated junctions as specific unannotated junctions
writes this subset to bam files 

cd /mnt/mindrinos/kevinroy/projects/COMPASS/
echo '' > extract_COMPASS_USJ_alignments_for_specific_junction_from_multiple_aligner_BAM.py
nano extract_COMPASS_USJ_alignments_for_specific_junction_from_multiple_aligner_BAM.py

"""

import pandas as pd
import pysam
import operator

from sys import argv
script, ALIGNMENTS_DIR, OUT_DIR, FASTA, INTRONS_FILE, ALIGNERS_FILE, SAMPLE_NAMES, SAMPLE_SUFFIX, COMBINED_SAMPLES_NAME, QUERY_CHROM, QUERY_START, QUERY_END, QUERY_STRAND = argv
for arg in argv:
    print(arg)

# COMPASS_DIR = '/Volumes/SPxDrive/COMPASS/' # '/mnt/mindrinos/kevinroy/projects/COMPASS/' # 
# GENOME_DIR = COMPASS_DIR + 'S288C_reference_genome_R64-2-1_20150113/'
# OUT_DIR = COMPASS_DIR + 'processed_data/'
# ALIGNMENTS_DIR = OUT_DIR + 'alignments/'
# COMPASS_JUNCTION_DIR = ALIGNMENTS_DIR + 'COMPASS_integration_USJ_only/'
# FASTA = GENOME_DIR + 'S288C_reference_sequence_R64-2-1_20150113_reformatted_chromosome_names.fasta'
# INTRONS_FILE = GENOME_DIR + 'S_cerevisiae_all_introns.tsv'
# ALIGNERS_FILE = COMPASS_DIR + 'sample_aligner_info.tsv'
# SAMPLE_NAMES = 'SRR5582776' #SRR5582777 SRR5582778 SRR5582779 SRR5582780 SRR5582781'
# COMBINED_SAMPLES_NAME = 'Aslanzadeh_combined'
# SAMPLE_SUFFIX = '_COMPASS_USJ_coord_sorted.bam'
# QUERY_CHROM = 'chrXVI'
# QUERY_START = 218647 - 1# 218647 
# QUERY_END = 218724 -1  # 218724
# QUERY_STRAND = '+'

query_start = int(QUERY_START)
query_end = int(QUERY_END)
sample_names = SAMPLE_NAMES.split(' ')

READ_NUM_PROGRESS = 1000

try:
    R64
except:    
    R64 =  pysam.FastaFile(FASTA) ## load_genome(GENOME_FASTA)

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

def get_ambiguous_junctions(chrom, start, stop, RNA_strand):
    '''
    takes intron coordinates
    checks for ambiguous junction based on nt upstream and downstream of exon-intron and intron-exon junctions
    returns list of ambiguous junctions
    '''
    junctions = [ (chrom, start, stop) ]
    idx = 1
    while start-idx > 0 and R64.fetch(chrom, start-idx, start-idx+1) == R64.fetch(chrom, stop-idx+1, stop-idx+2):
        ambiguous_intron = (chrom, start-idx, stop-idx)
        junctions.append(ambiguous_intron)
        idx += 1
    idx = 0
    # while stop-idx+1 < len(R64[chrom]) and R64[chrom][start-idx] == R64[chrom][stop-idx+1]:
    while stop-idx+1 < R64.get_reference_length(chrom) and R64.fetch(chrom, start-idx, start-idx+1) == R64.fetch(chrom, stop-idx+1, stop-idx+2):
        ambiguous_intron = (chrom, start-idx+1, stop-idx+1)
        junctions.append(ambiguous_intron)
        idx -= 1
    return junctions

CONSENSUS_5SS = ['G', 'T', 'A', 'TAC', 'G', 'TAC']
PENALTIES_5SS = [ 4,   3,   1,   1,    2,   1]

CONSENSUS_3SS = ['CT', 'A', 'G']
PENALTIES_3SS = [1,     3,   3]

    
def five_SS_score(query):
    return sum(PENALTIES_5SS[i] for i in range(len(query)) if query[i] not in CONSENSUS_5SS[i])

def three_SS_score(query):
    return sum(PENALTIES_3SS[i] for i in range(len(query)) if query[i] not in CONSENSUS_3SS[i])

def rc(seq):
    return seq.translate( str.maketrans("ACTG", "TGAC") )[::-1]

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


MISMATCH_DIST_FROM_JUNCTION_DISALLOWED = 10
MINIMUM_INTRON_LENGTH = 20
keys = 'flag, chrom, coord, cigar, NH, adjusted_introns, perfect_gapped_alignment, num_introns_found, alignment_score'.split(', ')          
common_intron_info = 'chrom, amb_start, amb_stop, five_SS, three_SS, RNA_strand, num_amb_junctions, annotated_junction, ann_5SS, ann_3SS, canonical_5SS, canonical_3SS, intron_size'.split(', ')
sample_specific_intron_info = 'intron_coords_adjusted, mismatch_near_junction, left_perfect_matches, right_perfect_matches, quality_score_at_left_junction, quality_score_at_right_junction'.split(', ')
adjusted_intron_keys = common_intron_info + sample_specific_intron_info

junctions_to_ambiguous_junctions = {}

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

ambiguous_junctions = get_ambiguous_junctions(QUERY_CHROM, query_start, query_end, QUERY_STRAND)

sample_name_to_USJ_read_pairs = {}
sample_name_to_USJ_start_pos = {}
USJ_reads_processed = 0
eof = False

## first get all read nums mapping to the junction

## need to fetch over the pile up region, check if the junction matches the exact junction,
## or any ambiguous junctions thereof, and then add the read nums to a set

aligner_df = pd.read_csv(ALIGNERS_FILE, sep = '\t') 
for index, row in aligner_df.iterrows():
    USJ_reads_processed = 0
    total_reads = 0
    aligner = row['aligner']
    print(aligner)
    USJ_start_pos_added = set([])
    for sample_name in sample_names:
        if sample_name not in sample_name_to_USJ_read_pairs:
            sample_name_to_USJ_read_pairs[sample_name] = set([])
        bam_infilename = ALIGNMENTS_DIR + aligner + '/' + sample_name + SAMPLE_SUFFIX
        bam_in = pysam.AlignmentFile(bam_infilename, "rb")
        print(bam_in)
        for read in bam_in.fetch(QUERY_CHROM, query_start - 500, query_end + 500):
            total_reads += 1
            read_num = int(read.query_name.split('_')[0])
            processed_read = process_alignment( read )
            USJ_start_pos = processed_read['coord']
            
            adjusted_introns = processed_read['adjusted_introns']
            for adjusted_intron in adjusted_introns:
                chrom, amb_start, amb_stop, RNA_strand = [adjusted_intron[e] for e in  'chrom, amb_start, amb_stop, RNA_strand'.split(', ')]
                # print (chrom, amb_start, amb_stop)
                if (chrom, amb_start, amb_stop) in ambiguous_junctions:
                    if USJ_start_pos not in USJ_start_pos_added:
                        sample_name_to_USJ_read_pairs[sample_name].add(read_num)
                        USJ_start_pos_added.add(USJ_start_pos)
                        USJ_reads_processed += 1
                   
                
        print('COMPASS has identified', USJ_reads_processed, 'read pairs for', ambiguous_junctions)
        print('from a total of', total_reads, 'from aligner:', aligner, 'and sample:', sample_name)
        bam_in.close()

for index, row in aligner_df.iterrows():
    USJ_reads_processed = 0
    
    total_reads = 0
    
    aligner = row['aligner']
    print(aligner)
    bam_infilename = ALIGNMENTS_DIR + aligner + '/' + sample_names[0] + SAMPLE_SUFFIX
    bam_in = pysam.AlignmentFile(bam_infilename, "rb")
        
    bam_prefix = ALIGNMENTS_DIR + aligner + '/' + aligner + '_' + COMBINED_SAMPLES_NAME + '_' +  QUERY_CHROM + '_' + str(query_start) + '_' + str(query_end)
    bam_outfilename = bam_prefix + '.bam'
    bam_out = pysam.AlignmentFile(bam_outfilename, "wb", template = bam_in)
    bam_sorted_outfilename = bam_prefix + '_sorted.bam'
    
    bam_in.close()
    
    for sample_name in sample_names:
        bam_infilename = ALIGNMENTS_DIR + aligner + '/' + sample_name + SAMPLE_SUFFIX
        bam_in = pysam.AlignmentFile(bam_infilename, "rb")
    
   
        for read in bam_in.fetch(QUERY_CHROM, query_start - 500, query_end + 500):
            total_reads += 1
            read_num = int(read.query_name.split('_')[0])
            read.query_name = read.query_name + '_' + sample_name
            if read_num in sample_name_to_USJ_read_pairs[sample_name]:
                print(aligner, read_num, 'read_num in USJ_read_pairs')
                bam_out.write(read)
                USJ_reads_processed += 1
        
        print('For the read recovery step:')
        print('COMPASS has retrieved', USJ_reads_processed, 'read pairs for', ambiguous_junctions)
        print('from a total of', total_reads, 'from aligner:', aligner)
        bam_in.close()
    bam_out.close()
        
    pysam.sort("-o", bam_sorted_outfilename, bam_outfilename)
    pysam.index(bam_sorted_outfilename)
