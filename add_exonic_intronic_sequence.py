#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  6 16:04:55 2020

@author: kevinroy

For info on 5SS and 3SS splice signals in humans:

A comprehensive survey of non-canonical splice sites in the human transcriptome 
Guillermo E. Parada, Roberto Munita, Cledi A. Cerda, Katia Gysling Author Notes
Nucleic Acids Research, Volume 42, Issue 16, 15 September 2014, Pages 10564–10578

Ranking noncanonical 5′ splice site usage by genome-wide RNA-seq analysis and splicing reporter assays
PMID: 30355602 PMCID: PMC6280755 DOI: 10.1101/gr.235861.118
noncanonical 5'ss usage ranking:
GC > TT > AT > GA > GG > CT

"""

import pybedtools
import pandas as pd
from itertools import product
import pysam

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

LOG_DIR = COMPASS_DIR + 'log/'
SEPARATE_ALIGNMENTS_DIR = COMPASS_DIR + 'separate_alignments/'
COMPASS_ALIGNMENTS_DIR = COMPASS_DIR + 'COMPASS_alignments/'
JUNCTION_ALIGNMENTS_DIR = COMPASS_DIR + 'junction_read_alignments/'
COMPASS_JUNCTIONS_DIR = COMPASS_DIR + 'COMPASS_junctions/'

BP_ADENOSINE_MISMATCH_PENALTY = 4
MIN_DISTANCE_BETWEEN_3SS_AND_BP = 1
MAX_DISTANCE_BETWEEN_3SS_AND_BP = 200 
# this is important to prevent excessive branchpoint searching on long introns
CONSENSUS_BP = 'NNCTNAN' # UACUAAC

genome_fasta = pysam.FastaFile(FASTA)

def list_ambiguous_junctions(chrom, start, stop, RNA_strand, genome_fasta):
    '''
    takes intron coordinates
    checks for ambiguous junction based on nt upstream and downstream of exon-intron and intron-exon junctions
    if ambiguous, adjusts based on closest match to 5SS consensus
    returns adjusted coordinates and other info on the intron: 
    five_SS, three_SS, RNA_strand, annotated_junction, ann_5SS, ann_3SS, canonical_5SS, canonical_3SS, intron_size, intron_coords_adjusted
    '''
    junctions = [(chrom, start, stop)]
    idx = 1
    num_left_amb_nt = 0
    while start-idx > 0 and genome_fasta.fetch(chrom, start - idx, start - idx + 1).upper() == genome_fasta.fetch(chrom, stop - idx + 1, stop - idx + 2).upper():
        ambiguous_intron = (chrom, start - idx, stop - idx)
        junctions.append(ambiguous_intron)
        idx += 1
        num_left_amb_nt += 1
    left_amb_nt = genome_fasta.fetch(chrom, start - num_left_amb_nt, start).upper()
    num_right_amb_nt = 0
    idx = 0
    while stop-idx+1 < genome_fasta.get_reference_length(chrom) and genome_fasta.fetch(chrom, start - idx, start - idx + 1).upper() == genome_fasta.fetch(chrom, stop - idx + 1, stop - idx + 2).upper():
        ambiguous_intron = (chrom, start - idx + 1, stop - idx + 1)
        junctions.append(ambiguous_intron)
        idx -= 1
        num_right_amb_nt += 1
    right_amb_nt = genome_fasta.fetch(chrom, start, start + num_right_amb_nt).upper()
    potential_5SS = []
    potential_3SS = []
    for junction in junctions:
        chrom, amb_start, amb_stop = junction
        if RNA_strand == '+':
            fiveSS_2nt = genome_fasta.fetch(chrom, amb_start, amb_start + 2).upper()
            threeSS_2nt = genome_fasta.fetch(chrom, amb_stop - 1, amb_stop + 1).upper()
        if RNA_strand == '-':
            fiveSS_2nt = rc(genome_fasta.fetch(chrom, amb_stop - 1, amb_stop + 1)).upper()
            threeSS_2nt = rc(genome_fasta.fetch(chrom, amb_start, amb_start + 2)).upper()
        potential_5SS.append(fiveSS_2nt)
        potential_3SS.append(threeSS_2nt)
    # The 5SS and 3SS dinucleotides can first be processed as pairs
    # favored in this order: 
    # ('GT','AG'), ('GC','AG'), ('AT','AC'), ('GA','AG'), ('GT','TG'), ('GG', 'AG'), ('AT', 'AG')
    # 
    # and then
    # be processed independently according to the order below.
    splice_site_pair_priority = ('GT','AG'), ('GC','AG'), ('AT','AC'), ('GA','AG'), ('GT','TG'), ('GG', 'AG'), ('AT', 'AG')
    splice_site_pairs = []
    prioritized_splice_site_pair_found = False
    for idx in range(len(potential_5SS)):
        splice_site_pairs.append( (potential_5SS[idx], potential_3SS[idx]) )
    for splice_site_pair in splice_site_pairs:
        if splice_site_pair in splice_site_pair_priority:
            prioritized_splice_site_pair_found = True
            adj_5SS, adj_3SS = splice_site_pair
            break
    if not prioritized_splice_site_pair_found:
        priority_for_5SS = 'GT', 'GC', 'AT', 'TT', 'GA', 'GG', 'CT'
        for fiveSS_priority in priority_for_5SS:
            best_5SS_indices = [i for i, fiveSS in enumerate(potential_5SS) if fiveSS == fiveSS_priority]
            if best_5SS_indices != []:
                break
        if len(best_5SS_indices) > 1:
            priority_for_3SS = 'AG', 'AC', 'TG', 'CG', 'GG', 'AT', 'AA'
            for threeSS_priority in priority_for_3SS:
                if threeSS_priority in [potential_3SS[i] for i in best_5SS_indices]:
                    break
            adj_5SS = fiveSS_priority
            adj_3SS = threeSS_priority
        elif len(best_5SS_indices) == 1:
            adj_5SS = potential_5SS[best_5SS_indices[0]]
            adj_3SS = potential_3SS[best_5SS_indices[0]]
        else:
            adj_5SS = potential_5SS[0]
            adj_3SS = potential_3SS[0]
    output = [potential_5SS, potential_3SS, adj_5SS, adj_3SS, left_amb_nt, right_amb_nt]
    return output

def rc(seq):
    return seq.translate(str.maketrans('ACGTacgtRYMKrymkVBHDvbhdNn ', 'TGCAtgcaYRKMyrkmBVDHbvdhNn '))[::-1]

def mismatches_from_branchpoint_consensus(consensus_branchpoint, query_string, BP_ADENOSINE_MISMATCH_PENALTY):
    '''
    input: two strings of equal length
    output: number of mismatches
    '''
    if len(consensus_branchpoint) != len(query_string):
        print(consensus_branchpoint, query_string)
        print(len(consensus_branchpoint), len(query_string))
        raise ValueError
    mismatches = 0
    if query_string[-2] != 'A':
        mismatches += BP_ADENOSINE_MISMATCH_PENALTY  ## branchpoint adenosine gets a significant penalty
    for idx in range(len(consensus_branchpoint)):
        if consensus_branchpoint[idx] != query_string[idx]:
            mismatches += 1
    return mismatches

def get_sequence_flanking_SS_and_likely_BP(chrom, start, end, strand, \
    max_US_perfect_matches, max_DS_perfect_matches, genome_fasta, \
        seqs_to_branchpoint_scores, MAX_DISTANCE_BETWEEN_3SS_AND_BP = MAX_DISTANCE_BETWEEN_3SS_AND_BP, \
            MIN_DISTANCE_BETWEEN_3SS_AND_BP = MIN_DISTANCE_BETWEEN_3SS_AND_BP):
    '''
    input: intron coordinates and the reference sequence as a dictionary of chroms
    output: a list of nucleotide content metrics surround the splice sites
    '''
    if strand == '+':
        # the max perfect matches for all junction-supporting reads,
        # rather than a fixed number of bp,
        # is extracted to enable calculating the edit distance.
        upstream_5SS = genome_fasta.fetch(chrom, start - max_US_perfect_matches, start)
        fiveSS_seq = genome_fasta.fetch(chrom, start, start + 6)
        downstream_5SS = genome_fasta.fetch(chrom, start + 6,start + 6 + max_DS_perfect_matches)
        
        upstream_3SS = genome_fasta.fetch(chrom, end - max_US_perfect_matches - 2, end - 2)
        threeSS_seq = genome_fasta.fetch(chrom, end - 2, end + 1)
        downstream_3SS = genome_fasta.fetch(chrom, end + 1, end + max_DS_perfect_matches + 1)
        seqs = [upstream_5SS, fiveSS_seq,  downstream_5SS, upstream_3SS, threeSS_seq, downstream_3SS ]

        US_3SS_motif_10nt = genome_fasta.fetch(chrom, end - 12, end - 2).upper()
        
        US_5SS_10nt = genome_fasta.fetch(chrom, max(start - 10, 1), start).upper()
        DS_5SS_10nt = genome_fasta.fetch(chrom, start, start + 10).upper()
        
        US_3SS_10nt = genome_fasta.fetch(chrom, end - 9, end + 1).upper()
        DS_3SS_10nt = genome_fasta.fetch(chrom, end + 1, min(end + 11, genome_fasta.get_reference_length(chrom)) ).upper()
        
    if strand == '-':
        downstream_3SS = genome_fasta.fetch(chrom, start - max_DS_perfect_matches, start)
        threeSS_seq = genome_fasta.fetch(chrom, start, start + 3)
        upstream_3SS = genome_fasta.fetch(chrom, start + 3, start + 3 + max_US_perfect_matches)
        
        downstream_5SS = genome_fasta.fetch(chrom, end - max_DS_perfect_matches - 5, end - 5)
        fiveSS_seq = genome_fasta.fetch(chrom, end - 5, end + 1)
        upstream_5SS = genome_fasta.fetch(chrom, end + 1, end + 1 + max_US_perfect_matches)
        seqs = [ rc(e) for e in (upstream_5SS, fiveSS_seq,  downstream_5SS, upstream_3SS, threeSS_seq, downstream_3SS  ) ]
        
        US_3SS_motif_10nt = rc( genome_fasta.fetch(chrom, start + 3, start + 13) ).upper()
        
        US_5SS_10nt = rc(genome_fasta.fetch(chrom, end + 1, min(end + 11, genome_fasta.get_reference_length(chrom))) ).upper()
        DS_5SS_10nt = rc(genome_fasta.fetch(chrom, end - 9, end + 1) ).upper()
        
        US_3SS_10nt = rc(genome_fasta.fetch(chrom, start, start + 10) ).upper()
        DS_3SS_10nt = rc(genome_fasta.fetch(chrom, max(start - 10, 1), start) ).upper()
    seqs = [e.upper() for e in seqs]            
    upstream_5SS, fiveSS_seq, downstream_5SS, upstream_3SS, threeSS_seq, downstream_3SS  =  seqs
    poly_U_count = US_3SS_motif_10nt.count('T')
    poly_Y_count = poly_U_count + US_3SS_motif_10nt.count('C')
    US_5SS_1 = US_5SS_10nt[-1]
    US_5SS_2 = US_5SS_10nt[-2]
    US_5SS_3 = US_5SS_10nt[-3]
    DS_3SS_1 = DS_3SS_10nt[0]
    DS_3SS_2 = DS_3SS_10nt[1]
    DS_3SS_3 = DS_3SS_10nt[2]
    if fiveSS_seq == 'GTATGT':
        fiveSS_type = 'GUAUGU'
    elif fiveSS_seq[:2] == 'GT':
        fiveSS_type = 'GU'
    elif fiveSS_seq[:2] == 'GC':
        fiveSS_type = 'GC'
    elif fiveSS_seq[:2] == 'AT':
        fiveSS_type = 'AU'
    else:
        fiveSS_type = 'non-GU'
    if threeSS_seq[-2:] == 'AG':
        threeSS_type = threeSS_seq
    elif threeSS_seq[-1] == 'G':
        threeSS_type = 'BG'
    elif threeSS_seq in ('AAT', 'CAT','TAT'):
        threeSS_type = 'HAU'
    elif threeSS_seq[-2:] == 'AC':
        threeSS_type = 'AC'
    else:
        threeSS_type = 'other'
    all_flanking_seq_info = [poly_U_count, poly_Y_count, US_5SS_1, US_5SS_2, \
        US_5SS_3, DS_3SS_1, DS_3SS_2, DS_3SS_3, US_5SS_10nt, DS_5SS_10nt, \
            US_3SS_10nt, DS_3SS_10nt] + seqs + [fiveSS_type, threeSS_type]
    mismatches_for_best_match_to_BP_consensus = 20
    best_BP_sequence = None
    if strand == '+':
        start_BP_search_coord = max(start + 6, end - MAX_DISTANCE_BETWEEN_3SS_AND_BP)
        end_BP_search_coord = end - 10 - MIN_DISTANCE_BETWEEN_3SS_AND_BP
        for coord in range(start_BP_search_coord, end_BP_search_coord):
            possible_BP_sequence = genome_fasta.fetch(chrom, coord, coord + 7).upper()
            mismatches = seqs_to_branchpoint_scores.get(possible_BP_sequence, 20)
            if mismatches < mismatches_for_best_match_to_BP_consensus:
                mismatches_for_best_match_to_BP_consensus = mismatches
                coord_for_best_match_to_BP_consensus = coord + 6
                BP_3SS_dist = end - coord_for_best_match_to_BP_consensus
                best_BP_sequence = possible_BP_sequence
                BP_5SS_dist = coord_for_best_match_to_BP_consensus - start
    else:
        start_BP_search_coord = min(end - 13, start + MAX_DISTANCE_BETWEEN_3SS_AND_BP)
        end_BP_search_coord = start + 3 + MIN_DISTANCE_BETWEEN_3SS_AND_BP
        for coord in range(start_BP_search_coord, end_BP_search_coord, -1):
            possible_BP_sequence = rc(genome_fasta.fetch(chrom, coord, coord + 7)).upper()
            mismatches = seqs_to_branchpoint_scores.get(possible_BP_sequence, 20)
            if mismatches < mismatches_for_best_match_to_BP_consensus:
                mismatches_for_best_match_to_BP_consensus = mismatches
                coord_for_best_match_to_BP_consensus = coord + 2
                BP_3SS_dist = coord_for_best_match_to_BP_consensus - start 
                best_BP_sequence = possible_BP_sequence
                BP_5SS_dist = end - coord_for_best_match_to_BP_consensus   
    try: 
        BP_3SS_dist
    except:
        print(chrom, start, end, strand, genome_fasta.get_reference_length(chrom))
    
    # canonical yeast (S. cerevisiae) motif (CUAAC)
    # 'CUCAC', 'CUGAC', 'CUAAC', 'CUAAU', 'UUAAC'
    # 102 branchpoint sites that conform to the branchpoint sequence UCCUURAY 
    # and splice sites of RNU12-dependent introns are spliced by the minor spliceosome
    if best_BP_sequence == None:
        CUNAN_branchpoint = False
        best_BP_sequence = 'NA'
        print('best_BP_sequence not found for:', chrom, start, end, strand)
        print('start_BP_search_coord:', start_BP_search_coord)
        print('end_BP_search_coord:', end_BP_search_coord)
    elif best_BP_sequence[2:7] in ['CUCAC', 'CUGAC', 'CUAAC', 'CUAAU', 'UUAAC']: 
        CUNAN_branchpoint = True 
    else: 
        CUNAN_branchpoint = False
    branchpoint_info =  [BP_3SS_dist, BP_5SS_dist, best_BP_sequence, coord_for_best_match_to_BP_consensus, CUNAN_branchpoint]
    return all_flanking_seq_info + branchpoint_info

all_7mers = [''.join(e) for e in list(product('ATCGN', repeat=7))]
seqs_to_branchpoint_scores = {}
for potential_branchpoint in all_7mers:
    seqs_to_branchpoint_scores[potential_branchpoint] = mismatches_from_branchpoint_consensus(CONSENSUS_BP, potential_branchpoint, BP_ADENOSINE_MISMATCH_PENALTY)

# for sample_name in 'Flag-PRPF18_2-9-D_1',: # SAMPLE="Flag-PRPF18_2-9-D_1"
for aligner in '', '_bbmap', '_HISAT2_default', '_HISAT2_noncanonical', '_STAR_default', '_STAR_noncanonical', '_MAGIC_BLAST', '_GSNAP':
    sample_aligner_name = sample_name + aligner
    print('starting ' + sample_aligner_name + '...')
    splice_junction_counts_infilename = COMPASS_JUNCTIONS_DIR + sample_aligner_name + '_COMPASS_splice_junctions.tsv'
    splice_junction_counts_outfilename = COMPASS_JUNCTIONS_DIR + sample_aligner_name + '_COMPASS_splice_junctions_with_sequence_info.tsv'

    junction_df = pd.read_csv(splice_junction_counts_infilename, sep = '\t')

    junction_df_filtered = junction_df.head(n=reads_to_process)
    # junction_df_filtered = junction_df.query('intron_size <= 2000 & intron_size >= 20') 
    # # ' & max_US_perfect_matches >= 10 & max_DS_perfect_matches >= 40 & COMPASS_counts >= 10')
    # junction_df_filtered = junction_df_filtered.reset_index()

    print('get_sequence_flanking_SS_and_likely_BP started...')
    new_cols = 'poly_U_count, poly_Y_count, US_5SS_1, US_5SS_2, US_5SS_3, DS_3SS_1, \
        DS_3SS_2, DS_3SS_3, US_5SS_10nt, DS_5SS_10nt, US_3SS_10nt, DS_3SS_10nt, upstream_5SS, \
            fiveSS_seq, downstream_5SS, upstream_3SS, threeSS_seq, downstream_3SS, fiveSS_type, \
                threeSS_type, BP_3SS_dist, BP_5SS_dist, best_BP_sequence, \
                    coord_for_best_match_to_BP_consensus, CUNAN_branchpoint'.split(', ')
    new_df = junction_df_filtered.apply(lambda x: get_sequence_flanking_SS_and_likely_BP(x.chrom, \
        x.adj_start, x.adj_stop, x.RNA_strand, x.max_US_perfect_matches, x.max_DS_perfect_matches, \
            genome_fasta, seqs_to_branchpoint_scores), axis=1)
    junction_df_filtered[new_cols] = pd.DataFrame(new_df.to_list(), index = junction_df_filtered.index)
    print('get_sequence_flanking_SS_and_likely_BP finished...')

    print('list_ambiguous_junctions started...')
    new_cols = 'potential_5SS, potential_3SS, adj_5SS, adj_3SS, left_amb_nt, right_amb_nt'.split(', ')
    new_df = junction_df_filtered.apply(lambda x: list_ambiguous_junctions(x.chrom, x.adj_start, x.adj_stop, x.RNA_strand, genome_fasta), axis=1)
    junction_df_filtered[new_cols] = pd.DataFrame(new_df.to_list(), index = junction_df_filtered.index)
    junction_df_filtered['sample_name'] = sample_name
    if aligner == '':
        junction_df_filtered['aligner'] = 'COMPASS'
    else:
        junction_df_filtered['aligner'] = aligner
    junction_df_filtered.to_csv(splice_junction_counts_outfilename, index=False, sep='\t')
    print('list_ambiguous_junctions finished...')
