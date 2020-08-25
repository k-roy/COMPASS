# -*- coding: utf-8 -*-
"""
Created on Thu Aug 29 07:55:58 2019

@author: kevinroy
"""

import pysam
import pandas as pd
#from sys import argv
GENOME_FASTA = '/Users/kevinroy/Google_Drive/yeast_genome_references/S288C_reference_sequence_R64-2-1_20150113_reformatted_chromosome_names.fasta'
#GENOME_FASTA = '/Volumes/SPxDrive/yeast_genome_references/S288C_reference_sequence_R64-2-1_20150113_reformatted_chromosome_names.fasta'
MINIMUM_INTRON_LENGTH = 20
MAXIMUM_INTRON_LENGTH = 2000
MAX_EDIT_DIST_TO_CONSIDER = 8 ## this is the max sum of bp involved in mismatches and indels (not including introns) to consider
MISMATCH_DIST_FROM_JUNCTION_DISALLOWED = 10
FIVE_SS_G1_MISMATCH_PENALTY = 2
FIVE_SS_T2_MISMATCH_PENALTY = 1
reads_to_process = 10**8
ALIGNMENTS_DIR = "/Volumes/TOSHIBA_EXT/alignments/"
OUT_DIR = "/Volumes/TOSHIBA_EXT/alignments/COMPASS_test2/"
#import os
#os.mkdir(OUT_DIR)
aligners = "bbmap", "STAR_default_annotated", "STAR_noncanonical_annotated", "HISAT2_noncanonical_annotated", "HISAT2_default_annotated", #  "HISAT2_noncanonical", "HISAT2_default"
suffix = "_name_sorted.bam"

try:
    R64
except:    
    R64 =  pysam.FastaFile(GENOME_FASTA) ## load_genome(GENOME_FASTA)

junctions_to_ambiguous_junctions = {}

annotated_intron_df = pd.read_csv('/Users/kevinroy/Google_Drive/yeast_genome_references/S_cerevisiae_all_introns.txt', sep = '\t') 
annotated_junctions = set([])
annotated_intron_df.columns
ambiguous_annotated_junctions = set([])
annotated_junc_to_type = {}
annotated_intron_df['start'] -= 1
annotated_intron_df['end'] -= 1
 
for index, row in annotated_intron_df.iterrows():
    chrom = row['chrom']
    if chrom == 'chrmt':
        chrom = 'chrMito'
    start = row['start']
    end = row['end']
    intron = (chrom, start, end)
    annotated_junctions.add( intron )
    annotated_junc_to_type[ intron] = row['intron_type']
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
print(len(ambiguous_annotated_junctions))
print(len(annotated_junctions) )

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
    
def five_SS_score(query, consensus):
    p1 = FIVE_SS_G1_MISMATCH_PENALTY if query[0] != 'G' else 0
    p2 = FIVE_SS_T2_MISMATCH_PENALTY if query[1] != 'T' else 0
    return hamming(query, consensus) + p1 + p2   

def rc(seq):
    return seq.translate( str.maketrans("ACTG", "TGAC") )[::-1]

def find_amb_junc(chrom, start, stop, RNA_strand):
    junctions = [ (chrom, start, stop) ]
    idx = 1
    # while start-idx > 0 and R64[chrom][start-idx] == R64[chrom][stop-idx+1]:
    while start-idx > 0 and R64.fetch(chrom, start-idx, start-idx+1) == R64.fetch(chrom, stop-idx+1, stop-idx+2):
        
        ambiguous_intron = (chrom, start-idx, stop-idx)
        junctions.append(ambiguous_intron)
        if start == 819326:
            print(ambiguous_intron)
        idx += 1
    idx = 0
    # while stop-idx+1 < len(R64[chrom]) and R64[chrom][start-idx] == R64[chrom][stop-idx+1]:
    while stop-idx+1 < R64.get_reference_length(chrom) and R64.fetch(chrom, start-idx, start-idx+1) == R64.fetch(chrom, stop-idx+1, stop-idx+2):
        ambiguous_intron = (chrom, start-idx+1, stop-idx+1)
        junctions.append(ambiguous_intron)
        if start == 819326:
            print(ambiguous_intron, R64.fetch(chrom, start-idx, start-idx+1), R64.fetch(chrom, stop-idx+1, stop-idx+2))
        idx -= 1
    if start == 819326:
        print(ambiguous_intron, R64.fetch(chrom, start-idx, start-idx+1), R64.fetch(chrom, stop-idx+1, stop-idx+2))
    potential_5SS = []
    potential_3SS = []
    hd = []
    for junction in junctions:
        chrom, amb_start, amb_stop = junction
        if RNA_strand == '+':
#            chrom = 'chrI' 
#            amb_start, amb_stop =  87387, 87499
            fiveSS = R64.fetch(chrom, amb_start, amb_start+6) # R64[chrom][amb_start:amb_start+6 ]
            potential_5SS.append( fiveSS )
            potential_3SS.append( R64.fetch(chrom, amb_stop-2, amb_stop+1) ) # R64[chrom][amb_stop-2:amb_stop+1 ] )
            hd.append(five_SS_score(fiveSS,  'GTATGT') )
        if RNA_strand == '-':
#            chrom = 'chrII'
#            amb_start, amb_stop =  47058,   47145
            fiveSS = rc( R64.fetch(chrom, amb_stop-5, amb_stop+1) )# R64[chrom][amb_stop-5:amb_stop+1] )
            potential_5SS.append( fiveSS )
            potential_3SS.append( rc( R64.fetch(chrom, amb_start, amb_start+3) ) )#  R64[chrom][amb_start:amb_start+3 ] ) )
            hd.append( five_SS_score(fiveSS, 'GTATGT') )
    
    if start == 819326:
        print(potential_5SS)
        print(potential_3SS)
    most_likely_intron = junctions[ hd.index( min(hd) ) ]  # return values.index(min(values))
    chrom, amb_start, amb_stop  = most_likely_intron
    if RNA_strand == '+':
        ann_5SS = ((annotated_intron_df['strand'] == RNA_strand) & (annotated_intron_df['start'] == amb_start)).any()
        ann_3SS = ((annotated_intron_df['strand'] == RNA_strand) & (annotated_intron_df['end'] == amb_start)).any()
    else:
        ann_3SS = ((annotated_intron_df['strand'] == RNA_strand) & (annotated_intron_df['start'] == amb_start)).any()
        ann_5SS = ((annotated_intron_df['strand'] == RNA_strand) & (annotated_intron_df['end'] == amb_start)).any()
    five_SS = potential_5SS[ hd.index( min(hd) ) ] 
    three_SS = potential_3SS[ hd.index( min(hd) ) ] 
    canonical_5SS = (five_SS[:2] == 'GT')
    canonical_3SS = (three_SS[1:] == 'AG')
    intron_size = abs(stop - start)
    annotated_junction = (most_likely_intron in annotated_junctions)
    return [chrom, amb_start, amb_stop, five_SS, three_SS, RNA_strand, annotated_junction, ann_5SS, ann_3SS, canonical_5SS, canonical_3SS, intron_size]

keys = 'flag, chrom, coord, cigar, NH, adjusted_introns, perfect_gapped_alignment, num_introns_found, alignment_score'.split(', ')          

def process_alignment(read):
    '''
    input: a list of tuples of (operation, length), where operation is one of “MIDNSHP=XB”
    '''
    if type(read) != pysam.libcalignedsegment.AlignedSegment or read.cigartuples == None:
        values = [None] * 5 + [ [], False, 0, 1000]
        return dict(zip(keys, values))
        
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
    perfect_matches_flanking_splice_sites  = []
    mapped_segment_lengths = []
    perfect_gapped_alignment = True  
    mismatch_near_junction = []
    num_introns_found = 0
    
    for cigar_idx in range(len(cigartuples)):
        bp, operation = cigartuples[cigar_idx]    
        if operation in ('=', 'X', 'D', 'M'):
            mapped_segment_length += bp
            coordinate += bp
        elif operation == 'N':
            total_gapped_bp += bp
            mapped_segment_lengths.append(mapped_segment_length)
            mapped_segment_length = 0
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
                splice_sites.append( (coordinate, coordinate + bp - 1) )
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
        intron = chrom, start, end
        if intron not in junctions_to_ambiguous_junctions:
            junctions_to_ambiguous_junctions[intron] = find_amb_junc(chrom, start, end, RNA_strand)
            ## chrom, amb_start, amb_stop, five_SS, three_SS, RNA_strand, annotated_junction, canonical_5SS, canonical_3SS, intron_size, mismatch_near_junction, left_perfect_matches, right_perfect_matches
        adjusted_introns.append( junctions_to_ambiguous_junctions[intron] + [mismatch_near_junction[idx], left_perfect_matches, right_perfect_matches ] )
    values = read.flag, read.reference_name, read.reference_start, read.cigarstring.replace('N', 'D'), read.get_tag('NH'), adjusted_introns, perfect_gapped_alignment, num_introns_found, alignment_score
    return dict(zip(keys, values))  

for sample in 'upf1', 'upf1_prp18':
    for replicate_num in '123':
        sample_name = sample + '_rep' + replicate_num
        print(sample_name)
        
        read_num = 0
        aligner_to_current_aligned_segment = {}
        aligner_to_current_read1_alignment = {}
        aligner_to_current_read2_alignment = {}
        aligner_to_current_read_num = {}
        aligners_to_bam = {}
        
        for aligner in aligners:
            bam_filename = ALIGNMENTS_DIR + aligner + '/' + sample_name + '_name_sorted.bam'
            aligners_to_bam[aligner] = pysam.AlignmentFile(bam_filename, "rb")
            aligner_to_current_read_num[aligner] = 0  
            
        complete_agreement_ungapped_bamfilename = OUT_DIR + sample_name + "_complete_agreement_ungapped.bam"
        complete_agreement_ungapped_bam_out = pysam.AlignmentFile(complete_agreement_ungapped_bamfilename, "wb", template = aligners_to_bam[aligners[0]] )
        
        complete_agreement_gapped_bamfilename = OUT_DIR + sample_name + "_complete_agreement_gapped.bam"
        complete_agreement_gapped_bam_out = pysam.AlignmentFile(complete_agreement_gapped_bamfilename, "wb", template = aligners_to_bam[aligners[0]] )
        
        annotated_gapped_alignments_outfilename = OUT_DIR + sample_name + "_annotated_gapped.tsv"
        annotated_gapped_alignments_outfile = open(annotated_gapped_alignments_outfilename, 'w')
        
        unannotated_gapped_alignments_outfilename = OUT_DIR + sample_name + "_unannotated_gapped.tsv"
        unannotated_gapped_alignments_outfile = open(unannotated_gapped_alignments_outfilename, 'w')
        
#        incongruent_ungapped_alignments_outfilename = OUT_DIR + sample_name + "_incongruent_ungapped.tsv"
#        incongruent_ungapped_alignments_outfile = open(incongruent_ungapped_alignments_outfilename, 'w')
#        
        intron_found_to_alignments_outfiles = {True: annotated_gapped_alignments_outfile, False: unannotated_gapped_alignments_outfile}
        
        #####

        complete_ungapped_agreement_counts = 0
        complete_gapped_agreement_counts = 0
        ungapped_agreement_counts = 0
        gapped_agreement_counts = 0
        
        columns = ['read_pair_ID', 'aligner_has_min_score', 'min_score_has_annotated_intron', 'min_score_has_no_intron', 'aligner', 'score',  'num_introns_found', 'perfect_gapped_alignment', 'chrom', 'amb_start', 'amb_stop', 'five_SS', 'three_SS', 'RNA_strand', 'annotated_junction', 'ann_5SS', 'ann_3SS', 'canonical_5SS', 'canonical_3SS', 'intron_size', 'mismatch_near_junction', 'left_perfect_matches', 'right_perfect_matches', 'read1_flag', 'read1_chrom', 'read1_coord', 'read1_cigar', 'read1_NH', 'read2_flag', 'read2_chrom', 'read2_coord', 'read2_cigar', 'read2_NH' ]
        header = '\t'.join(columns) + '\n'
        [e.write(header) for e in intron_found_to_alignments_outfiles.values()]
        eof = False
        
        while read_num <= reads_to_process and not eof:
            if read_num % 10000 == 0:
                print(read_num, 'reads processed by COMPASS')
            read_num += 1
            gapped_alignment_present = False
            for aligner in aligners:
                while aligner_to_current_read_num[aligner] < read_num:
                    try:
                        aligned_segment = next( aligners_to_bam[aligner] )
                    except StopIteration:
                        eof = True
                        break
                        break
                        break
                    if not aligned_segment.is_secondary:
                        aligner_to_current_read_num[aligner] = int(aligned_segment.query_name.split('_')[0])
                        aligner_to_current_aligned_segment[aligner] = aligned_segment
                        
                read1_first = True
                aligned_segment = aligner_to_current_aligned_segment[aligner]
                aligner_to_current_read1_alignment[aligner] = None
                aligner_to_current_read2_alignment[aligner] = None
                
                if aligner_to_current_read_num[aligner] == read_num:
                    if aligned_segment.is_read1:
                        aligner_to_current_read1_alignment[aligner] = aligned_segment 
                    else:
                        aligner_to_current_read2_alignment[aligner] = aligned_segment
                        read1_first = False
                    try:
                        aligned_segment = next( aligners_to_bam[aligner] )
                    except StopIteration:
                        eof = True
                        break
                        break
                        break
                    aligner_to_current_read_num[aligner] = int(aligned_segment.query_name.split('_')[0])
                
                    if aligner_to_current_read_num[aligner] == read_num and not aligned_segment.is_secondary:
                        if read1_first and aligned_segment.is_read1:
                            print('two adjacent read 1s')
                        if aligned_segment.is_read1:
                            aligner_to_current_read1_alignment[aligner] = aligned_segment 
                        else:
                            aligner_to_current_read2_alignment[aligner] = aligned_segment
        
                    
            read1_alignments = []
            read2_alignments = []
            perfect_lst = []
            num_intron_found_lst = []
            alignment_scores = []
            R1_lst = []
            R2_lst = []
            for aligner in aligners:
                R1 = process_alignment( aligner_to_current_read1_alignment[aligner] )
                R2 = process_alignment( aligner_to_current_read2_alignment[aligner] )
                R1_lst.append(R1)
                R2_lst.append(R2)
        #     values = read.flag, read.reference_name, read.reference_start, read.cigarstring.replace('N', 'D'), read.get_tag('NH'), adjusted_introns, perfect_gapped_alignment, num_introns_found, alignment_score
                # keys =  'flag, chrom, coord, cigar, NH, adjusted_introns, annotated_junction, mismatch_near_junction, perfect_gapped_alignment, intron_size, canonical_5SS, canonical_3SS, intron_found, score'.split(', ')
                read1_alignments.append( tuple([R1[e] for e in keys[:5]]) )
                read2_alignments.append( tuple([R2[e] for e in keys[:5]]) )
                # # adjusted_introns --> chrom, amb_start, amb_stop, five_SS, three_SS, RNA_strand, annotated_junction, canonical_5SS, canonical_3SS, intron_size, mismatch_near_junction , left_perfect_matches, right_perfect_matches
                alignment_scores.append( R1['alignment_score'] + R2['alignment_score'] )
                num_intron_found_lst.append( max( R1['num_introns_found'], R2['num_introns_found'] ))
                if R1['num_introns_found'] > 0 and R2['num_introns_found'] > 0:
                    perfect_gapped_alignment = R1['perfect_gapped_alignment'] or R2['perfect_gapped_alignment']
                elif R1['num_introns_found'] > 0:
                    perfect_gapped_alignment = R1['perfect_gapped_alignment']
                elif R2['num_introns_found'] > 0:
                    perfect_gapped_alignment = R2['perfect_gapped_alignment']
                else:
                    perfect_gapped_alignment = False
                perfect_lst.append( perfect_gapped_alignment )
            best_score = min(alignment_scores)
            min_score_has_annotated_intron = False  ## does any best alignment have annotated intron
            min_score_has_no_intron = False  ## does any best alignment have no intron
            if best_score < MAX_EDIT_DIST_TO_CONSIDER:
                if len(set(read1_alignments)) == len(set(read2_alignments)) == 1:
                    if min(num_intron_found_lst) > 0:
                        complete_gapped_agreement_counts += 1
                        [complete_agreement_gapped_bam_out.write(e[aligner]) for e in (aligner_to_current_read1_alignment, aligner_to_current_read2_alignment) ]
                    else:
                        complete_ungapped_agreement_counts += 1
                        [complete_agreement_ungapped_bam_out.write(e[aligner]) for e in (aligner_to_current_read1_alignment, aligner_to_current_read2_alignment) ]
                for idx in range(len(aligners)):            
                    R1 = R1_lst[idx]
                    R2 = R2_lst[idx]
                    for e in R1, R2:
                        if e['adjusted_introns'] != []:
                        # chrom, amb_start, amb_stop, five_SS, three_SS, RNA_strand, annotated_junction, canonical_5SS, canonical_3SS = adjusted_intron  
                        ## read 1 and read 2 info entered on separate lines if both contain a gapped alignment, 
                        ## later on these can be further processed and integrated by groupby operations
                        ## if a read contains more than intron (quite rare), it is entered on multiple lines
                            for intron in e['adjusted_introns']:
                                if alignment_scores[idx] == best_score and intron[6]:
                                    min_score_has_annotated_intron = True ## does any best alignment have annotated intron
                    if num_intron_found_lst[idx] == 0 and alignment_scores[idx] == best_score:
                        min_score_has_no_intron = True ## does any best alignment have no intron
                for idx in range(len(aligners)):            
                    read_pair_info = [read_num, alignment_scores[idx] == best_score, min_score_has_annotated_intron, min_score_has_no_intron ] + [ f[idx] for f in (aligners, alignment_scores, num_intron_found_lst, perfect_lst) ]
                    read1_info = list(read1_alignments[idx])
                    read2_info = list(read2_alignments[idx])
                    R1 = R1_lst[idx]
                    R2 = R2_lst[idx]
                    for e in R1, R2:
                        if e['adjusted_introns'] != []:
                        # chrom, amb_start, amb_stop, five_SS, three_SS, RNA_strand, annotated_junction, canonical_5SS, canonical_3SS = adjusted_intron  
                        ## read 1 and read 2 info entered on separate lines if both contain a gapped alignment, 
                        ## later on these can be further processed and integrated by groupby operations
                        ## if a read contains more than intron (quite rare), it is entered on multiple lines
                            for intron in e['adjusted_introns']:
                                # # # adjusted_introns --> chrom, amb_start, amb_stop, five_SS, three_SS, RNA_strand, annotated_junction, canonical_5SS, canonical_3SS, intron_size, mismatch_near_junction , left_perfect_matches, right_perfect_matches
        
                                output = read_pair_info + intron + read1_info + read2_info
                                output = [str(e) for e in output]
                                if max(num_intron_found_lst) > 0:
                                    intron_found_to_alignments_outfiles[ min_score_has_annotated_intron ].write('\t'.join(output) + '\n')
                                if max(num_intron_found_lst) > 0:
                                    gapped_agreement_counts += 1
                                else:
                                    ungapped_agreement_counts += 1
                    if R1['adjusted_introns'] == R2['adjusted_introns'] == []:
                        ## no gapped alignments found, 13 NA values
                        output = read_pair_info + ['NA']*15 + read1_info + read2_info
                        output = [str(e) for e in output]
                        if max(num_intron_found_lst) > 0:
                            intron_found_to_alignments_outfiles[ min_score_has_annotated_intron ].write('\t'.join(output) + '\n')
        
                    
        
        [e.close() for e in (complete_agreement_ungapped_bam_out, complete_agreement_gapped_bam_out, annotated_gapped_alignments_outfile, unannotated_gapped_alignments_outfile) ]
        
        [aligners_to_bam[aligner].close() for aligner in aligners]
                
        [print(e) for e in [complete_ungapped_agreement_counts, complete_gapped_agreement_counts, ungapped_agreement_counts, gapped_agreement_counts]]
