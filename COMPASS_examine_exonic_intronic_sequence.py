#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  6 16:04:55 2020

@author: kevinroy
"""

import pandas as pd
import pysam
import pyranges as pr
import matplotlib as plt

COMPASS_DIR = '/u/project/guillom/kevinh97/COMPASS/' # '/mnt/mindrinos/kevinroy/projects/COMPASS/' #
GENOME_DIR = COMPASS_DIR + 'S288C_reference_genome_R64-2-1_20150113/'

OUT_DIR = COMPASS_DIR + 'processed_data/'
ALIGNMENTS_DIR = OUT_DIR + 'alignments/'
COMPASS_JUNCTION_DIR = ALIGNMENTS_DIR + 'COMPASS_integration_full_output/'
FASTA = GENOME_DIR + 'Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa'
INTRONS_FILE = '/u/project/guillom/kevinh97/COMPASS/S288C_reference_genome_R64-2-1_20150113/saccharomyces_cerevisiae_R64-2-1_20150113_introns_no_chr.tsv'
ALIGNERS_FILE = '/u/project/guillom/kevinh97/sample_aligner_info.tsv'
sample_suffix = '_COMPASS_splice_junctions.tsv'

try:
    R64
except:    
    R64 =  pysam.FastaFile(FASTA) ## load_genome(GENOME_FASTA)

def list_ambiguous_junctions(chrom, start, stop, RNA_strand):
    '''
    takes intron coordinates
    checks for ambiguous junction based on nt upstream and downstream of exon-intron and intron-exon junctions
    if ambiguous, adjusts based on closest match to 5SS consensus
    returns adjusted coordinates and other info on the intron: 
    five_SS, three_SS, RNA_strand, annotated_junction, ann_5SS, ann_3SS, canonical_5SS, canonical_3SS, intron_size, intron_coords_adjusted
    '''
    junctions = [ (chrom, start, stop) ]
    idx = 1
    # while start-idx > 0 and R64[chrom][start-idx] == R64[chrom][stop-idx+1]:
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
    potential_5SS = []
    potential_3SS = []
    for junction in junctions:
        chrom, amb_start, amb_stop = junction
        if RNA_strand == '+':
            fiveSS_2nt = R64.fetch(chrom, amb_start, amb_start+2) # R64[chrom][amb_start:amb_start+6 ]
            threeSS_2nt = R64.fetch(chrom, amb_stop-1, amb_stop+1)
        if RNA_strand == '-':
            fiveSS_2nt = rc( R64.fetch(chrom, amb_stop-1, amb_stop+1) )# R64[chrom][amb_stop-5:amb_stop+1] )
            threeSS_2nt = rc( R64.fetch(chrom, amb_start, amb_start+2) ) #  R64[chrom][amb_start:amb_start+3 ] ) )
        potential_5SS.append( fiveSS_2nt )
        potential_3SS.append( threeSS_2nt )
        
    priority_for_5SS = 'GT', 'GA', 'GC', 'GG', 'AT'
    for fiveSS_priority in priority_for_5SS:
        best_5SS_indices = [i for i, fiveSS in enumerate(potential_5SS) if fiveSS == fiveSS_priority]
        if best_5SS_indices != []:
            break
        
    if len(best_5SS_indices) > 1:
        priority_for_3SS = 'AG', 'CG', 'TG', 'GG', 'AT', 'AC'
        for threeSS_priority in priority_for_3SS:
            if threeSS_priority in [potential_3SS[i] for i in best_5SS_indices]:
                break
        adj_5SS = fiveSS_priority
        adj_3SS = threeSS_priority
    elif len(best_5SS_indices) == 1:
        adj_5SS = potential_5SS[ best_5SS_indices[0] ]
        adj_3SS = potential_3SS[ best_5SS_indices[0] ]
    else:
        adj_5SS = potential_5SS[0]
        adj_3SS = potential_3SS[0]
         
    return [potential_5SS, potential_3SS, adj_5SS, adj_3SS]

def poly_U_PWMS(upstream_3SS_seq):
    '''
    takes 10 bp sequences upstream of the 3'SS
    returns a score from the PWSM from Ma et al. doi:  10.1155/2011/212146
    '''
    U_score = 0
    U_score_by_position = [0.2, 0.7114, .8332, .9090, 1.2812, 1.0412, 0.6793, 0.4574, 0.7532, 0.0850]
    if len(U_score_by_position) != len(upstream_3SS_seq):
        raise ValueError
    for idx in range(len(upstream_3SS_seq)):
        nt = upstream_3SS_seq[idx]
        if nt == 'T':
            U_score += U_score_by_position[idx]
    #print(U_score, 'U_score', upstream_3SS_seq)
    return U_score

def rc(seq):
    return seq.translate( str.maketrans("ACTG", "TGAC") )[::-1]

def mismatches_from_branchpoint_consensus(consensus_branchpoint, query_string):
    '''
    input: two strings of equal length
    output: number of mismatches
    '''
    BP_ADENOSINE_MISMATCH_PENALTY = 4
    if len(consensus_branchpoint) != len(query_string):
        print(consensus_branchpoint, query_string)
        print(len(consensus_branchpoint), len(query_string) )
        raise ValueError
    mismatches = 0
    if query_string[-2] != 'A':
        mismatches += BP_ADENOSINE_MISMATCH_PENALTY  ## branchpoint adenosine gets a significant penalty
    for idx in range(len(consensus_branchpoint)):
        if consensus_branchpoint[idx] != query_string[idx]:
            mismatches += 1
    return mismatches

def get_sequence_flanking_SS(chromosome, start, end, strand, max_US_perfect_matches, max_DS_perfect_matches):
    '''
    input: intron coordinates and the reference sequence as a dictionary of chromosomes
    output: a list of nucleotide content metrics surround the splice sites
    '''
    if strand == '+':
        upstream_5SS = R64.fetch(chromosome, start - max_US_perfect_matches, start)
        fiveSS_seq = R64.fetch(chromosome, start, start + 6)
        downstream_5SS = R64.fetch(chromosome, start + 6,start + 6 + max_DS_perfect_matches)
        
        upstream_3SS = R64.fetch(chromosome, end - max_US_perfect_matches - 2, end - 2)
        threeSS_seq = R64.fetch(chromosome, end - 2, end + 1)
        downstream_3SS = R64.fetch(chromosome, end + 1, end + max_DS_perfect_matches + 1)
        seqs = [ upstream_5SS, fiveSS_seq,  downstream_5SS, upstream_3SS, threeSS_seq, downstream_3SS  ]
        
        US_3SS_motif_10nt = R64.fetch(chromosome, end - 12, end - 2)
        
        US_5SS_10nt = R64.fetch(chromosome, start - 10, start)
        DS_5SS_10nt = R64.fetch(chromosome, start, start + 10)
        
        US_3SS_10nt = R64.fetch(chromosome, end - 9, end + 1)
        DS_3SS_10nt = R64.fetch(chromosome, end + 1, end + 11)
        
    if strand == '-':
        downstream_3SS = R64.fetch(chromosome, start - max_DS_perfect_matches, start)
        threeSS_seq = R64.fetch(chromosome, start, start + 3)
        upstream_3SS = R64.fetch(chromosome, start + 3, start + 3 + max_US_perfect_matches)
        
        downstream_5SS = R64.fetch(chromosome, end - max_DS_perfect_matches - 5, end - 5)
        fiveSS_seq = R64.fetch(chromosome, end - 5, end + 1)
        upstream_5SS = R64.fetch(chromosome, end + 1, end + 1 + max_US_perfect_matches)
        seqs = [ rc(e) for e in (upstream_5SS, fiveSS_seq,  downstream_5SS, upstream_3SS, threeSS_seq, downstream_3SS  ) ]
        
        US_3SS_motif_10nt = rc( R64.fetch(chromosome, start + 3, start + 13) )
        
        US_5SS_10nt = rc(R64.fetch(chromosome, end + 1, end + 11) )
        DS_5SS_10nt = rc(R64.fetch(chromosome, end - 9, end + 1) )
        
        US_3SS_10nt = rc(R64.fetch(chromosome, start, start + 10) )
        DS_3SS_10nt = rc(R64.fetch(chromosome, start - 10, start) )
        
    upstream_5SS, fiveSS_seq,  downstream_5SS, upstream_3SS, threeSS_seq, downstream_3SS  =  seqs
    poly_U_count = US_3SS_motif_10nt.count('T')
    poly_Y_count = poly_U_count + US_3SS_motif_10nt.count('C')
    U_score = poly_U_PWMS(US_3SS_motif_10nt)
    US_5SS_1 = US_5SS_10nt[-1]
    US_5SS_2 = US_5SS_10nt[-2]
    US_5SS_3 = US_5SS_10nt[-3]
    DS_3SS_1 = DS_3SS_10nt[0]
    DS_3SS_2 = DS_3SS_10nt[1]
    DS_3SS_3 = DS_3SS_10nt[2]
    if fiveSS_seq == 'GTATGT':
        fiveSS_type = 'GTATGT'
    elif fiveSS_seq[:2] == 'GT':
        fiveSS_type = 'GT'
    elif fiveSS_seq[:2] == 'AT':
        fiveSS_type = 'AT'
    else:
        fiveSS_type = 'non-GT'
    if threeSS_seq[-2:] == 'AG':
        threeSS_type = threeSS_seq
    elif threeSS_seq[-1] == 'G':
        threeSS_type = 'BG'
    elif threeSS_seq in ('AAT', 'CAT','TAT'):
        threeSS_type = 'HAT'
    elif threeSS_seq[-2:] == 'AC':
        threeSS_type = 'AC'
    else:
        threeSS_type = 'other'
    all_flanking_seq_info = [poly_U_count, poly_Y_count, US_5SS_1, US_5SS_2, US_5SS_3, DS_3SS_1, DS_3SS_2, DS_3SS_3, US_5SS_10nt, DS_5SS_10nt, US_3SS_10nt, DS_3SS_10nt] + seqs + [fiveSS_type, threeSS_type]
    
    mismatches_for_best_match_to_BP_consensus = 20
    MIN_DISTANCE_BETWEEN_3SS_AND_BP = 1
    CONSENSUS_BP = 'TACTAAC'
    # 1.)ACTAA, 2.)PYTPAYP, 3.)YTPAY
    #  BP_priority = 'TACTAAC', 'ACTAAC', 'ACTAA', 'RYTRAYR', 'YTRAY'
    best_BP_sequence = None
    if strand == '+':
        for coord in range(start + 6, end - 10 - MIN_DISTANCE_BETWEEN_3SS_AND_BP):
            possible_BP_sequence = R64.fetch(chromosome, coord, coord + 7)
            mismatches = mismatches_from_branchpoint_consensus(CONSENSUS_BP, possible_BP_sequence)
            if mismatches < mismatches_for_best_match_to_BP_consensus:
                mismatches_for_best_match_to_BP_consensus = mismatches
                coord_for_best_match_to_BP_consensus = coord + 6
                BP_3SS_dist = end - coord_for_best_match_to_BP_consensus
                best_BP_sequence = possible_BP_sequence
                BP_5SS_dist = coord_for_best_match_to_BP_consensus - start
    else:
        for coord in range( end - 13, start + 3 + MIN_DISTANCE_BETWEEN_3SS_AND_BP, -1):
            possible_BP_sequence = rc(R64.fetch(chromosome, coord, coord + 7))
     
            mismatches = mismatches_from_branchpoint_consensus(CONSENSUS_BP, possible_BP_sequence)
            if mismatches < mismatches_for_best_match_to_BP_consensus:
                mismatches_for_best_match_to_BP_consensus = mismatches
                coord_for_best_match_to_BP_consensus = coord + 2
                BP_3SS_dist = coord_for_best_match_to_BP_consensus - start 
                best_BP_sequence = possible_BP_sequence
                BP_5SS_dist = end - coord_for_best_match_to_BP_consensus   
    intron_length = end - start
    try: 
        BP_3SS_dist
    except:
        print(chromosome, start, end, strand, R64.get_reference_length(chromosome) )
    
    if best_BP_sequence[2:7] in ['CTAAC', 'CTAAT', 'CTGAC', 'CTGAT', 'TTAAC', 'TTAAT', 'TTGAC', 'TTGAT']: 
        YTRAY_branchpoint = True 
    else: 
        YTRAY_branchpoint = False
    branchpoint_info =  [intron_length, BP_3SS_dist, BP_5SS_dist, best_BP_sequence, coord_for_best_match_to_BP_consensus, YTRAY_branchpoint]
    return all_flanking_seq_info + branchpoint_info

# U_score = poly_U_PWMS(upstream_3SS_seq):
# poly_U_count, poly_Y_count, US_5SS_1, US_5SS_2, US_5SS_3, DS_3SS_1, DS_3SS_2, DS_3SS_3, upstream_5SS, fiveSS_seq,  downstream_5SS, upstream_3SS, threeSS_seq, downstream_3SS, fiveSS_type, threeSS_type  =  
# mismatches = mismatches_from_branchpoint_consensus(consensus_branchpoint, query_string)
# intron_length, BP_3SS_dist, BP_5SS_dist, best_BP_sequence, coord_for_best_match_to_BP_consensus = get_closest_upstream_branchpoint(chromosome, start, end, strand)
# 'SRR5582776 SRR5582777 SRR5582778 SRR5582779 SRR5582780 SRR5582781'.split(' ')  # 
# sample_names = 'SRR9130287 SRR9130288 SRR9130289 SRR9130290 SRR9130291 SRR9130292'.split(' ') # Roy et al.
# sample_names = 'SRR5041706 SRR5041707 SRR5041708 SRR5041709'.split(' ')  # Talkish et al.
# sample_names = 'SRR5582776 SRR5582777 SRR5582778 SRR5582779 SRR5582780 SRR5582781'.split(' ')  # Aslanzadeh et al.
sample_names = 'RRP6_1SC__'.split(' ') # 'SRR5582778', 'SRR5582779', 'SRR5582780', 'SRR5582781'
for sample_name in sample_names:
    print(sample_name)
    junction_filename = COMPASS_JUNCTION_DIR + sample_name + sample_suffix
    junction_df = pd.read_csv(junction_filename, sep = '\t')
    junction_df.columns
    junction_df['intron_size'] = junction_df.amb_stop - junction_df.amb_start
    junction_df_filtered = junction_df.query('intron_size <= 2000 & intron_size >= 20') # ' & max_US_perfect_matches >= 10 & max_DS_perfect_matches >= 40 & COMPASS_counts >= 10')
    junction_df_filtered = junction_df_filtered.reset_index()
    new_cols = 'poly_U_count, poly_Y_count, US_5SS_1, US_5SS_2, US_5SS_3, DS_3SS_1, DS_3SS_2, DS_3SS_3, US_5SS_10nt, DS_5SS_10nt, US_3SS_10nt, DS_3SS_10nt, upstream_5SS, fiveSS_seq, downstream_5SS, upstream_3SS, threeSS_seq, downstream_3SS, fiveSS_type, threeSS_type, intron_length, BP_3SS_dist, BP_5SS_dist, best_BP_sequence, coord_for_best_match_to_BP_consensus, YTRAY_branchpoint'.split(', ')
    new_df = junction_df_filtered.apply( lambda x: get_sequence_flanking_SS(x.chrom, x.amb_start, x.amb_stop, x.RNA_strand, x.max_US_perfect_matches, x.max_DS_perfect_matches), axis=1)
    junction_df_filtered[new_cols] = pd.DataFrame(new_df.to_list(), index = junction_df_filtered.index)
    junction_df_filtered.columns
    
    new_cols = 'potential_5SS, potential_3SS, adj_5SS, adj_3SS'.split(', ')
    new_df = junction_df_filtered.apply( lambda x: list_ambiguous_junctions(x.chrom, x.amb_start, x.amb_stop, x.RNA_strand), axis=1)
    junction_df_filtered[new_cols] = pd.DataFrame(new_df.to_list(), index = junction_df_filtered.index)

    new_sample_suffix = '_COMPASS_splice_junctions_with_sequence_info.tsv'                                           
    junction_filename = COMPASS_JUNCTION_DIR + sample_name + new_sample_suffix
    junction_df_filtered.to_csv(junction_filename, sep = '\t' )
