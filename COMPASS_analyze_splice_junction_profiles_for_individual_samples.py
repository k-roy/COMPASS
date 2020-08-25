# -*- coding: utf-8 -*-
"""
Created on Tue Dec 12 11:46:02 2017

@author: kevinroy

Integrates filtered splice junctions from multiple samples (processed with the script splice_junctions_from_BBMAP_and_STAR_SAM.py)
Analyzes splice junctions overlapping annotated introns, must pass filters for maximum 5'SS-alt5'SS distance and max intron size
Calculates splicing efficiency (SE) and fraction of annotated splicing (FAnS) for each junction
Checks for ambiguous junction assignments due to identical nucleotides upstream and/or downstream of mapped junction
Adjusts alignment to favor the 5´SS closest in hamming distance to the consensus GTATGT, with extra penalties for not matching the first G or T
Calculates the most likely branchpoint sequence for each 3´SS
Analyzes sequence flanking each splice site
Designates sample-specific vs. common attributes for each junction parameter
Writes a data table for each sample in the dataset with each row representing a different junction, and each column attributes of that junction
"""
import ast, operator

PARENT_DIR= '/Volumes/SPxDrive/' # '/Users/kevinroy/Google_Drive/' # 
DIR=PARENT_DIR + 'splice_junction_analysis/'
SAMPLES='BY_upf1del_prp18del_3 BY_upf1del_prp18del_2 upf1_8.14.15 upf1_8.17.15 WT prp18 upf1 upf1_prp18'.split(' ') 
genome_file = PARENT_DIR + 'yeast_genome_references/saccharomyces_cerevisiae_sequence.fasta'

#script_placeholder, GENOME_FASTA, GFF_ANNOTATIONS, ALIGNMENTS_ALIGNERS, JUNCTION_STATISTICS_OUTFILENAME, INTEGRATED_ALIGNMENT_OUTFILENAME, FILTERED_SAM_DIR, SAMPLE_PREFIX = argv
SAMPLE_NAME = 'SRR5582776'
GENOME_FASTA = '/Users/kevinroy/Google_Drive/yeast_genome_references/S288C_reference_sequence_R64-2-1_20150113_reformatted_chromosome_names.fasta'
GFF_ANNOTATIONS = '/Users/kevinroy/Google_Drive/yeast_genome_references/saccharomyces_cerevisiae_annotations.gff'
DIR = '/Volumes/MyPassportforMac_5TB/processed_sequencing_data/fast_slow_Pol_II_splicing/'
ALIGNMENTS_ALIGNERS = DIR + SAMPLE_NAME + '_unannotated_junctions.sam,STAR;' + DIR + SAMPLE_NAME + '_unannotated_junction_reads_bbmap.sam,BBMap;' + DIR + 'HISAT2_alignments/' + SAMPLE_NAME + '_output.sam,HISAT2'
JUNCTION_STATISTICS_OUTFILENAME = DIR + SAMPLE_NAME + '_junction_statistics.txt'
INTEGRATED_ALIGNMENT_OUTFILENAME = DIR + SAMPLE_NAME + '_integrated_alignments.txt'
FILTERED_SAM_DIR = DIR
SAMPLE_PREFIX = SAMPLE_NAME


MINIMUM_NUMBER_OF_UNIQUELY_MAPPED_READS = 1
MAX_BP_FOR_ALT_5SS_TO_ANNOTATED_5SS = 1000
MAXIMUM_INTRON_SIZE = 2000
MINIMUM_INTRON_SIZE = 20
MAX_MISMATCHES_FOR_FIVE_PRIME_SPLICE_SITE = 6 ## allow all mapped five prime splice sites
CONSENSUS_FIVE_PRIME_SPLICE_SITE = 'GTATGT'
FLANKING_SEQ_BP = 10
FIVE_SS_G1_MISMATCH_PENALTY = 4
FIVE_SS_T2_MISMATCH_PENALTY = 3
BP_ADENOSINE_MISMATCH_PENALTY = 4
CONSENSUS_BP = 'TACTAAC'

## load annotated splicing efficiency files
## these were generated from STAR mapping, followed by extracting only unspliced reads to obtain the unspliced intron reads
def load_annotated_splicing_efficiency(infilename):
    with open(infilename, 'r') as infile:
        header = infile.readline()
        annotated_introns = {}
        for line in infile:
            systematic_gene_name, common_gene_name, chromosome, start, end, strand, unique_junction_reads, intron_reads5, intron_reads3, SE5, SE3 = line.strip().split()
            start, end, unique_junction_reads, intron_reads5, intron_reads3 = [int(e) for e in [start, end, unique_junction_reads, intron_reads5, intron_reads3] ]
            intron_reads = (intron_reads5 + intron_reads3) / 2.0 ## intron reads come from aggregating unspliced reads from STAR mapping             
            if SE5 != 'NA':
                SE5 = float(SE5)
            if SE3 != 'NA':
                SE3 = float(SE3)
            if common_gene_name == 'GCR1':
                start, end = 412257, 412995
                print('GCR1 intron coords changed')
            
            annotated_introns[ (chromosome, start, end, strand) ] = systematic_gene_name, common_gene_name, unique_junction_reads, intron_reads5, intron_reads3, SE5, SE3, intron_reads
        infile.close()
    return annotated_introns

def load_genome(genome_fasta_filename):
    '''
    input: genome in fasta format
    output: dictionary with chromosome names as keys, where the value is the chromosomal sequence as a string
    '''
    genome_fasta_file = open(genome_fasta_filename, 'r')
    genome_dict = {}
    for line in genome_fasta_file:
        if line[0] == '>':
            current_chromosome = line.strip()[1:]
            genome_dict[current_chromosome] = ''
        else:
            genome_dict[current_chromosome] += line.strip()
    genome_fasta_file.close()
    print(genome_file + " genome file loaded")
    return genome_dict

def rev_comp(DNA):
    '''
    input: DNA string
    output: reverse complement string
    '''
    rev_comp_DNA = ''
    comp_bases = {'A':'T', 'T':'A', 'G':'C', 'C':'G'}
    for base in DNA[::-1]:
        rev_comp_DNA += comp_bases[base]
    return rev_comp_DNA

def load_splice_junctions(SJ_filename):
    '''
    input: SJ.out.tab produced by splice_junctions_from_BBMAP_and_STAR_SAM.py (file ends in '_combined_bbmap_STAR_reformatted_cigar_SJ.out.tab')
    output: sample-specific dictionary with junction_ID keys, and values a tuple containing data on each junction
    '''
    lines_to_process = 100000
    total_uniquely_mapped_gapped_reads = 0
    infile = open(SJ_filename, 'r')
    splice_junction_dict = {}
    line_num = 0
    info = infile.readline().strip('\n').split('\t')
    for line in infile:
        line_num += 1
        if line_num <= lines_to_process:
            info = line.strip('\n').split('\t')
            chromosome, start, end, strand, reads_without_mismatch_near_junction, sorted_coordinate_cigar_combos, sorted_mismatches_dict, reads_with_perfect_gapped_alignments, reads_with_gapped_alignments, unique_junction_counts, max_overhang_upstream_junction, max_overhang_perfect_matches_upstream_junction, max_overhang_downstream_junction, max_overhang_perfect_matches_downstream_junction, BBMAP_and_STAR, BBMAP_only, STAR_only, BBMAP_intron_not_in_STAR, flag_99, flag_147, flag_83, flag_163 = info
            start, end, reads_without_mismatch_near_junction, reads_with_perfect_gapped_alignments, reads_with_gapped_alignments, unique_junction_counts, max_overhang_upstream_junction, max_overhang_perfect_matches_upstream_junction, max_overhang_downstream_junction, max_overhang_perfect_matches_downstream_junction, BBMAP_and_STAR, BBMAP_only, STAR_only, BBMAP_intron_not_in_STAR, flag_99, flag_147, flag_83, flag_163 = [int(e) for e in (start, end, reads_without_mismatch_near_junction, reads_with_perfect_gapped_alignments, reads_with_gapped_alignments, unique_junction_counts, max_overhang_upstream_junction, max_overhang_perfect_matches_upstream_junction, max_overhang_downstream_junction, max_overhang_perfect_matches_downstream_junction, BBMAP_and_STAR, BBMAP_only, STAR_only, BBMAP_intron_not_in_STAR, flag_99, flag_147, flag_83, flag_163)]
            total_uniquely_mapped_gapped_reads += reads_with_gapped_alignments         
            if chromosome != 'chrMito' and chromosome != 'chrmt':
                junction_ID = (chromosome, start, end, strand)
                sorted_coordinate_cigar_combos = dict( ast.literal_eval(sorted_coordinate_cigar_combos) )
                sorted_mismatches_dict = dict( ast.literal_eval(sorted_mismatches_dict) )         
                splice_junction_dict[junction_ID] = chromosome, start, end, strand, reads_without_mismatch_near_junction, sorted_coordinate_cigar_combos, sorted_mismatches_dict, reads_with_perfect_gapped_alignments, reads_with_gapped_alignments, unique_junction_counts, max_overhang_upstream_junction, max_overhang_perfect_matches_upstream_junction, max_overhang_downstream_junction, max_overhang_perfect_matches_downstream_junction, BBMAP_and_STAR, BBMAP_only, STAR_only, BBMAP_intron_not_in_STAR, flag_99, flag_147, flag_83, flag_163
    infile.close()
    print('splice_junction_dict loaded: ', SJ_filename, ' with ', total_uniquely_mapped_gapped_reads, ' total_uniquely_mapped_gapped_reads')
    return splice_junction_dict, unique_junction_counts

## the purpose of analyzing five prime splice sites is to resolve ambiguous junctions by favoring the most likely 5'SS
def mismatches_from_five_prime_splice_site(consensus_splice_site, query_string):
    '''
    input: two strings of equal length
    output: number of mismatches
    '''
    mismatches = 0
    if query_string[0] != 'G':  ## increased penalty for not matching the first base G
        mismatches += FIVE_SS_G1_MISMATCH_PENALTY
    if query_string[1] != 'T':  ## increased penalty for not matching the second base T
        mismatches += FIVE_SS_T2_MISMATCH_PENALTY
    for idx in range(2, len(consensus_splice_site)):
        if consensus_splice_site[idx] != query_string[idx]:
            mismatches += 1
    return mismatches

def check_for_ambiguous_junction_sequence(chromosome, query_start, query_end, strand, R64_genome):
    '''
    input: intron coordinates and the reference sequence as a dictionary of chromosomes
    output: revised intron coordinates
    if there is identical sequence downstream of the three prime splice site and downstream of the five prime splice site, or upstream of the three prime splice site and upstream of the five prime splice site,
    pick the sequence that results in the closest match to the consensus splice site GTATGT
    '''
    idx = 0
    best_match_idx = 0
    coordinate_adjustment = 0
    # print query_start, query_end, R64_genome[chromosome][query_start-1: query_start + 5], R64_genome[chromosome][query_end: query_end + 6]
    ## convert 1-based coordinates to 0-based coordinates
    ## check if the first base of the five prime splice site is the same as the first base of the downstream exon, then check second, third, and so forth (for positive strand introns)
    while R64_genome[chromosome][query_end + idx] == R64_genome[chromosome][query_start -1 + idx]:
        output_outfile.write( 'checking bases downstream of ' + str(query_start) + ' and ' + str(query_end) + '\n' )
        output_outfile.write( str(idx) + ' ' + R64_genome[chromosome][query_start-1 + idx] + ' ' + R64_genome[chromosome][query_end + idx]  + '\n' )
        idx += 1
    ambiguous_downstream_bases = idx
    idx = 1
    ## check if the last base of the intron is the same as the last base of the upstream exon, then check second to last, third to last, and so forth (for positive strand introns)
    while R64_genome[chromosome][query_end - idx] == R64_genome[chromosome][query_start-1 - idx]:
        output_outfile.write('checking bases upstream of ' + str(query_start) + ' and ' + str(query_end) + '\n' )
        output_outfile.write(str(idx) + ' ' + R64_genome[chromosome][query_start-1 - idx] + ' ' +  R64_genome[chromosome][query_end - idx] + '\n' )
        idx += 1
    ambiguous_upstream_bases = idx - 1
    num_mismatches_to_consensus = []
    idx_of_num_mismatches_to_consensus = []
    ambiguous_bases_around_junction = ambiguous_downstream_bases + ambiguous_upstream_bases 
    if ambiguous_upstream_bases > 0 or ambiguous_downstream_bases > 0:
        if strand == '+':
            for idx in range(-ambiguous_upstream_bases, ambiguous_downstream_bases+1):
                possible_five_prime_splice_site = R64_genome[chromosome][query_start + idx - 1: query_start + idx + 5]
                mismatches = mismatches_from_five_prime_splice_site(CONSENSUS_FIVE_PRIME_SPLICE_SITE, possible_five_prime_splice_site)
                num_mismatches_to_consensus.append(mismatches)
                idx_of_num_mismatches_to_consensus.append(idx)
                output_outfile.write(possible_five_prime_splice_site + ' splice site score: ' + str(mismatches) + '\n' )
            best_match = min(num_mismatches_to_consensus)
            output_outfile.write( str(num_mismatches_to_consensus) + ' num_mismatches_to_consensus with index for coordinate adjustment ' + str(idx_of_num_mismatches_to_consensus) + '\n' )
            if num_mismatches_to_consensus.count(best_match) > 1:
                output_outfile.write(chromosome + ' ' + str(query_start) + ' ' + str(query_end) + ' ' + strand + ' ' + 'ambiguous junction with two or more possible five prime splice sites of equal distance from consensus!' + '\n' )
            best_match_idx = num_mismatches_to_consensus.index(best_match)
            coordinate_adjustment = idx_of_num_mismatches_to_consensus[best_match_idx]
            output_outfile.write('best_match mismatches from consensus: ' + str(best_match) + ' with coordinate adjustment of ' + str(coordinate_adjustment) + '\n' )
            new_query_start = query_start + coordinate_adjustment
            new_query_end  = query_end + coordinate_adjustment
            if coordinate_adjustment != 0:
                output_outfile.write(chromosome + ' ' + str(query_start) + ' ' + str(query_end) + ' ' + strand + ' ' +  'junction reassigned due to five prime splice sites nearby closer to consensus.' + '\n' )
                output_outfile.write(chromosome + ' ' + str(query_start) + ' ' + str(query_end) + ' ' + strand + ' ' + 'five prime splice site changed from ' + R64_genome[chromosome][query_start - 1: query_start + 5] + ' to '  + R64_genome[chromosome][new_query_start - 1: new_query_start + 5]  + '\n' )
                output_outfile.write(chromosome + ' ' + str(query_start) + ' ' + str(query_end) + ' ' + strand + ' ' +  'three prime splice site changed from ' + R64_genome[chromosome][query_end - 3: query_end] + ' to '  + R64_genome[chromosome][new_query_end - 3: new_query_end]  + '\n' )
        elif strand == '-':
            for idx in range(-ambiguous_upstream_bases, ambiguous_downstream_bases + 1):
                possible_five_prime_splice_site = rev_comp( R64_genome[chromosome][query_end + idx - 6: query_end + idx  ] )
                mismatches = mismatches_from_five_prime_splice_site(CONSENSUS_FIVE_PRIME_SPLICE_SITE, possible_five_prime_splice_site)
                num_mismatches_to_consensus.append(mismatches)
                idx_of_num_mismatches_to_consensus.append(idx)
                output_outfile.write(possible_five_prime_splice_site + ' splice site score: ' + str(mismatches) + '\n' )
            best_match = min(num_mismatches_to_consensus)
            output_outfile.write( str(num_mismatches_to_consensus) + ' num_mismatches_to_consensus with index for coordinate adjustment ' + str(idx_of_num_mismatches_to_consensus) + '\n' )            
            if num_mismatches_to_consensus.count(best_match) > 1:
                output_outfile.write(chromosome + ' ' + str(query_start) + ' ' + str(query_end) + ' ' + strand + ' ' + 'ambiguous junction with two or more possible five prime splice sites of equal distance from consensus!' + '\n' )
            best_match_idx = num_mismatches_to_consensus.index(best_match)
            coordinate_adjustment = idx_of_num_mismatches_to_consensus[best_match_idx]
            output_outfile.write('best_match mismatches from consensus: ' + str(best_match) + ' with coordinate adjustment of ' + str(coordinate_adjustment) + '\n' )
            new_query_start = query_start + coordinate_adjustment
            new_query_end = query_end + coordinate_adjustment
            coordinate_adjustment = -coordinate_adjustment
            if coordinate_adjustment != 0:
                output_outfile.write(chromosome + ' ' + str(query_start) + ' ' + str(query_end) + ' ' + strand + ' ' +  'junction reassigned due to five prime splice sites nearby closer to consensus.' + '\n' )
                output_outfile.write(chromosome + ' ' + str(query_start) + ' ' + str(query_end) + ' ' + strand + ' ' +  'five prime splice site changed from ' + rev_comp ( R64_genome[chromosome][query_end - 6: query_end] ) + ' to '  +  rev_comp ( R64_genome[chromosome][new_query_end - 6: new_query_end] )  + '\n' )
                output_outfile.write(chromosome + ' ' + str(query_start) + ' ' + str(query_end) + ' ' + strand + ' ' +  'three prime splice site changed from ' + rev_comp ( R64_genome[chromosome][query_start - 1: query_start + 2] ) + ' to '  + rev_comp ( R64_genome[chromosome][new_query_start - 1: new_query_start + 2] )  + '\n' )
        return new_query_start, new_query_end, coordinate_adjustment, ambiguous_bases_around_junction
    return query_start, query_end, coordinate_adjustment, ambiguous_bases_around_junction

def mismatches_from_branchpoint_consensus(consensus_branchpoint, query_string):
    '''
    input: two strings of equal length
    output: number of mismatches
    '''
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

def get_closest_upstream_branchpoint(chromosome, query_start, query_end, strand, R64_genome, systematic_gene_name, common_gene_name):
    '''
    input: intron coordinates and the reference sequence as a dictionary of chromosomes
    output: The distance to a sequence within the intron that is closest to the branchpoint consensus
    '''
    gene_names_to_check = 'HAC1', 'GCR1'# 'RPL40B', 'DYN2' # 'YCL002C', 'BET4', 'TAD3', 'SUS1'
    mismatches_for_best_match_to_BP_consensus = 20
    MIN_DISTANCE_BETWEEN_3SS_AND_BP = 1
    best_BP_sequence = None
    if strand == '+':
        for coord in range(query_start + 8, query_end - 10 - MIN_DISTANCE_BETWEEN_3SS_AND_BP):
            possible_BP_sequence = R64_genome[chromosome][coord:coord+7]
            mismatches = mismatches_from_branchpoint_consensus(CONSENSUS_BP, possible_BP_sequence)
            if mismatches < mismatches_for_best_match_to_BP_consensus:
                mismatches_for_best_match_to_BP_consensus = mismatches
                coord_for_best_match_to_BP_consensus = coord + 6
                BP_3SS_dist = query_end - coord_for_best_match_to_BP_consensus
                best_BP_sequence = possible_BP_sequence
                BP_5SS_dist = coord_for_best_match_to_BP_consensus - query_start
#                if len(possible_BP_sequence) < 7 or common_gene_name in gene_names_to_check:
#                    print(common_gene_name, chromosome, query_start, query_end, coord_for_best_match_to_BP_consensus, possible_BP_sequence, best_BP_sequence)

    else:
        for coord in range( query_end - 16, query_start + 3 + MIN_DISTANCE_BETWEEN_3SS_AND_BP, -1):
            possible_BP_sequence = rev_comp(R64_genome[chromosome][coord:coord+7])
     
            mismatches = mismatches_from_branchpoint_consensus(CONSENSUS_BP, possible_BP_sequence)
            if mismatches < mismatches_for_best_match_to_BP_consensus:
                mismatches_for_best_match_to_BP_consensus = mismatches
                coord_for_best_match_to_BP_consensus = coord + 2
                BP_3SS_dist = coord_for_best_match_to_BP_consensus - query_start 
                best_BP_sequence = possible_BP_sequence
                BP_5SS_dist = query_end - coord_for_best_match_to_BP_consensus
#                if len(possible_BP_sequence) < 7 or common_gene_name in gene_names_to_check:
#                    print(common_gene_name, chromosome, query_start, query_end, possible_BP_sequence, best_BP_sequence)
#       
    intron_length = query_end - query_start
    
#    if best_BP_sequence == None:
#        print(chromosome, query_start, query_end, strand,systematic_gene_name, common_gene_name )
    return [intron_length, BP_3SS_dist, BP_5SS_dist, best_BP_sequence, coord_for_best_match_to_BP_consensus]

def get_sequence_flanking_SS(chromosome, start, end, strand, R64_genome):
    '''
    input: intron coordinates and the reference sequence as a dictionary of chromosomes
    output: a list of nucleotide content metrics surround the splice sites
    '''    
    if strand == '+':
        upstream_5SS = R64_genome[chromosome][start - FLANKING_SEQ_BP - 1:start - 1]
        fiveSS_seq = R64_genome[chromosome][start - 1:start + 5]
        downstream_5SS = R64_genome[chromosome][start + 5:start + 5 + FLANKING_SEQ_BP]
        
        upstream_3SS = R64_genome[chromosome][end - FLANKING_SEQ_BP - 3:end - 3]
        threeSS_seq = R64_genome[chromosome][end - 3:end ]
        downstream_3SS = R64_genome[chromosome][end:end + FLANKING_SEQ_BP]
        seqs = [ upstream_5SS, fiveSS_seq,  downstream_5SS, upstream_3SS, threeSS_seq, downstream_3SS  ]

    if strand == '-':
        downstream_3SS = R64_genome[chromosome][start - FLANKING_SEQ_BP -1 :start -1 ] 
        threeSS_seq = R64_genome[chromosome][start -1:start + 2]
        upstream_3SS = R64_genome[chromosome][start + 2:start + 2 + FLANKING_SEQ_BP]
        
        downstream_5SS = R64_genome[chromosome][end - FLANKING_SEQ_BP - 6:end - 6]
        fiveSS_seq = R64_genome[chromosome][end - 6:end ]
        upstream_5SS = R64_genome[chromosome][end:end + FLANKING_SEQ_BP]
        seqs = [ rev_comp( e) for e in (upstream_5SS, fiveSS_seq,  downstream_5SS, upstream_3SS, threeSS_seq, downstream_3SS  ) ]

    upstream_5SS, fiveSS_seq,  downstream_5SS, upstream_3SS, threeSS_seq, downstream_3SS  =  seqs
    poly_U_count = upstream_3SS.count('T')
    poly_Y_count = poly_U_count + upstream_3SS.count('C')
    US_5SS_1 = upstream_5SS[-1]
    US_5SS_2 = upstream_5SS[-2]
    US_5SS_3 = upstream_5SS[-3]
    DS_3SS_1 = downstream_3SS[0]
    DS_3SS_2 = downstream_3SS[1]
    DS_3SS_3 = downstream_3SS[2]
    if fiveSS_seq == 'GTATGT':
        fiveSS_type = 'GTATGT'
    elif fiveSS_seq[:2] == 'GT':
        fiveSS_type = 'GT'
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
    all_flanking_seq_info = [poly_U_count, poly_Y_count, US_5SS_1, US_5SS_2, US_5SS_3, DS_3SS_1, DS_3SS_2, DS_3SS_3] + seqs + [fiveSS_type, threeSS_type]
    return all_flanking_seq_info
    
def combine_dicts(dict1, dict2):
    '''
    input: two dicts with overlapping keys with values as integers
    output: a union of the two dicts, with the values from common keys summed
    '''  
    new_dict = {}
    for e in dict1:
        if e in dict2:
            new_dict[e]  = dict1[e] + dict2[e]
        else:
            new_dict[e]  = dict1[e]
    for e in dict2:
        if e not in new_dict:
            new_dict[e] = dict2[e]
    return new_dict

def combine_clashing_junction_info(sample_specific_alternative_info, previous_sample_specific_alternative_info  ):
    '''
    input: two dicts with overlapping keys with values as integers
    output: a union of the two dicts, with the values from common keys summed
    '''  
    coordinate_adjustment, strand, annotated_and_alternative_strands_agree, reads_without_mismatch_near_junction, sorted_coordinate_cigar_combos, sorted_mismatches_dict, reads_with_perfect_gapped_alignments, reads_with_gapped_alignments, unique_junction_counts, max_overhang_upstream_junction, max_overhang_perfect_matches_upstream_junction, max_overhang_downstream_junction, max_overhang_perfect_matches_downstream_junction, BBMAP_and_STAR, BBMAP_only, STAR_only, BBMAP_intron_not_in_STAR, flag_99, flag_147, flag_83, flag_163 = sample_specific_alternative_info
    previous_coordinate_adjustment, previous_strand, previous_annotated_and_alternative_strands_agree, previous_reads_without_mismatch_near_junction, previous_sorted_coordinate_cigar_combos, previous_sorted_mismatches_dict, previous_reads_with_perfect_gapped_alignments, previous_reads_with_gapped_alignments, previous_unique_junction_counts, previous_max_overhang_upstream_junction, previous_max_overhang_perfect_matches_upstream_junction, previous_max_overhang_downstream_junction, previous_max_overhang_perfect_matches_downstream_junction, previous_BBMAP_and_STAR, previous_BBMAP_only, previous_STAR_only, previous_BBMAP_intron_not_in_STAR, previous_flag_99, previous_flag_147, previous_flag_83, previous_flag_163 = previous_sample_specific_alternative_info
    
    updated_coordinate_adjustment = previous_coordinate_adjustment + ',' + coordinate_adjustment
    updated_strand = previous_strand + ',' + strand 
    updated_annotated_and_alternative_strands_agree = annotated_and_alternative_strands_agree  
    udpated_sorted_coordinate_cigar_combos = combine_dicts(sorted_coordinate_cigar_combos, previous_sorted_coordinate_cigar_combos )
    udpated_sorted_mismatches_dict = combine_dicts( sorted_mismatches_dict, previous_sorted_mismatches_dict )
    updated_reads_with_perfect_gapped_alignments = int(previous_reads_with_perfect_gapped_alignments) + int(reads_with_perfect_gapped_alignments)
    updated_reads_with_gapped_alignments = int(previous_reads_with_gapped_alignments) + int(reads_with_gapped_alignments)
    updated_unique_junction_counts = int(previous_unique_junction_counts) + int(unique_junction_counts)
    updated_reads_without_mismatch_near_junction = int(previous_reads_without_mismatch_near_junction) + int(reads_without_mismatch_near_junction)
    updated_BBMAP_and_STAR = int(previous_BBMAP_and_STAR) + int(BBMAP_and_STAR)
    updated_BBMAP_only = int(previous_BBMAP_only) + int(BBMAP_only)
    updated_STAR_only = int(previous_STAR_only) + int(STAR_only)
    updated_BBMAP_intron_not_in_STAR = int(previous_BBMAP_intron_not_in_STAR) + int(BBMAP_intron_not_in_STAR)
    updated_flag_99 = int(previous_flag_99) + int(flag_99)
    updated_flag_147 = int(previous_flag_147) + int(flag_147)
    updated_flag_83 = int(previous_flag_83) + int(flag_83)
    updated_flag_163 = int(previous_flag_163) + int(flag_163)
    
    updated_max_overhang_upstream_junction = max( int(previous_max_overhang_upstream_junction), int(max_overhang_upstream_junction) )
    updated_max_overhang_perfect_matches_upstream_junction = max( int(previous_max_overhang_perfect_matches_upstream_junction), int(max_overhang_perfect_matches_upstream_junction) )
    updated_max_overhang_downstream_junction = max( int(previous_max_overhang_downstream_junction), int(max_overhang_downstream_junction) )
    updated_max_overhang_perfect_matches_downstream_junction = max( int(previous_max_overhang_perfect_matches_downstream_junction), int(max_overhang_perfect_matches_downstream_junction) )
    
    return updated_coordinate_adjustment, updated_strand, updated_annotated_and_alternative_strands_agree, updated_reads_without_mismatch_near_junction, udpated_sorted_coordinate_cigar_combos, udpated_sorted_mismatches_dict, updated_reads_with_perfect_gapped_alignments, updated_reads_with_gapped_alignments, updated_unique_junction_counts, updated_max_overhang_upstream_junction, updated_max_overhang_perfect_matches_upstream_junction, updated_max_overhang_downstream_junction, updated_max_overhang_perfect_matches_downstream_junction, updated_BBMAP_and_STAR, updated_BBMAP_only, updated_STAR_only, updated_BBMAP_intron_not_in_STAR, updated_flag_99, updated_flag_147, updated_flag_83, updated_flag_163

def splice_site_analysis(annotated_introns_info, junction_dict, annotated_and_alternative_outfilename):
    '''
    takes a junctions dictionary and the annotated intron stats
    writes a file, with each row a different splice junction overlapping a known intron, 
    and each column attributes of that splice junction, along with info on the annotated and unspliced reads for that junction
    '''
    gene_names_to_check = 'HAC1', 'GCR1'
    outfile = open(annotated_and_alternative_outfilename, 'w')
    ## pre-process the annotated introns
    header = 'annotated_junction_name, annotated_strand, systematic_gene_name, common_gene_name, RPG, annotated_intron_length, annotated_BP_3SS_dist, annotated_BP_5SS_dist, annotated_best_BP_sequence, annotated_coord_for_best_match_to_BP_consensus, annotated_poly_U_count, annotated_poly_Y_count, annotated_US_5SS_1, annotated_US_5SS_2, annotated_US_5SS_3, annotated_DS_3SS_1, annotated_DS_3SS_2, annotated_DS_3SS_3, annotated_upstream_5SS_seq, annotated_fiveSS_seq, annotated_downstream_5SS_seq, annotated_upstream_3SS_seq, annotated_threeSS_seq, annotated_downstream_3SS_seq, annotated_fiveSS_type, annotated_threeSS_type, junction_ID_name, chromosome, adjusted_start, adjusted_end, ambiguous_bases_around_junction, alt_5SS, alt_3SS, intron_length, BP_3SS_dist, BP_5SS_dist, best_BP_sequence, coord_for_best_match_to_BP_consensus, poly_U_count, poly_Y_count, US_5SS_1, US_5SS_2, US_5SS_3, DS_3SS_1, DS_3SS_2, DS_3SS_3, upstream_5SS, fiveSS_seq, downstream_5SS, upstream_3SS, threeSS_seq, downstream_3SS, fiveSS_type, threeSS_type, annotated_unique_junction_reads, intron_reads, total_spliced_and_unspliced_reads_associated_with_annotated_intron, annotated_SE, coordinate_adjustment, alternative_strand, SE, fraction_of_annotated_splicing, reads_without_mismatch_near_junction, coordinate_cigar_combos, mismatches_tallies, reads_with_perfect_gapped_alignments, reads_with_gapped_alignments, unique_junction_counts, max_overhang_upstream_junction, max_overhang_perfect_matches_upstream_junction, max_overhang_downstream_junction, max_overhang_perfect_matches_downstream_junction, BBMAP_and_STAR, BBMAP_only, STAR_only, BBMAP_intron_not_in_STAR, flag_99, flag_147, flag_83, flag_163'.split(', ')
    print(len(header))    
    outfile.write('\t'.join(header) + '\n')
    annotated_intron_to_nearest_alternative_junctions = {}    
    for junction_ID in junction_dict:
        chromosome, intron_start, intron_end, strand = junction_ID
        ## for the 5'SS for each splice junction, find the closest annotated five prime SS
        if MINIMUM_INTRON_SIZE <= intron_end - intron_start <= MAXIMUM_INTRON_SIZE:
            closest_five_SS = None
            closest_five_SS_dist = 100000000
            for annotated_intron_info in annotated_introns_info:
                annotated_chromosome, annotated_start, annotated_end, annotated_strand = annotated_intron_info
                if chromosome == annotated_chromosome and strand == annotated_strand:
                    if annotated_strand == '+':
                        start_diff = abs(annotated_start - intron_start)
                    else:
                        start_diff = abs(annotated_end - intron_end)
                    if start_diff < closest_five_SS_dist:
                        closest_five_SS_dist = start_diff
                        closest_five_SS = annotated_chromosome, annotated_start, annotated_end, annotated_strand 
            if closest_five_SS_dist <= MAX_BP_FOR_ALT_5SS_TO_ANNOTATED_5SS:
                annotated_chromosome, annotated_start, annotated_end, annotated_strand = closest_five_SS
                annotated_junction_ID = annotated_chromosome, annotated_start, annotated_end, annotated_strand
                ## introns must overlap at least partially with an annotated intron
                if intron_start < annotated_end and intron_end > annotated_start:
                    if annotated_junction_ID not in annotated_intron_to_nearest_alternative_junctions:
                        annotated_intron_to_nearest_alternative_junctions[annotated_junction_ID] = []
                    annotated_intron_to_nearest_alternative_junctions[annotated_junction_ID].append( junction_ID )
    ## first pass through each annotated intron exhibiting alternative splicing junction events            
    for annotated_junction_ID in annotated_intron_to_nearest_alternative_junctions:
        annotated_chromosome, annotated_start, annotated_end, annotated_strand = annotated_junction_ID
        annotated_junction_name = '_'.join( [str(e) for e in annotated_junction_ID ])
        systematic_gene_name, common_gene_name, annotated_unique_junction_reads, intron_reads5, intron_reads3, annotated_SE5, annotated_SE3, intron_reads = annotated_introns_info[ annotated_junction_ID ]        
        
#        if common_gene_name in gene_names_to_check:
#            print(annotated_introns_info[ annotated_junction_ID ])
#            print(intron_reads)
        
    
        RPG = False        
        if common_gene_name[:2] == 'RP' and common_gene_name[:3] != 'RPO':
            RPG = True
        annotated_intron_BP_info = get_closest_upstream_branchpoint(annotated_chromosome, annotated_start, annotated_end, annotated_strand, R64_genome, systematic_gene_name, common_gene_name)
        annotated_intron_all_flanking_seq_info = get_sequence_flanking_SS(annotated_chromosome, annotated_start, annotated_end, annotated_strand, R64_genome)   
        unique_junctions_after_site_adjustment = {}
        
        total_spliced_and_unspliced_reads_associated_with_annotated_intron = intron_reads
        
        for junction_ID in annotated_intron_to_nearest_alternative_junctions[annotated_junction_ID]:
            chromosome, intron_start, intron_end, strand = junction_ID      
            chromosome, start, end, strand, reads_without_mismatch_near_junction, sorted_coordinate_cigar_combos, sorted_mismatches_dict, reads_with_perfect_gapped_alignments, reads_with_gapped_alignments, unique_junction_counts, max_overhang_upstream_junction, max_overhang_perfect_matches_upstream_junction, max_overhang_downstream_junction, max_overhang_perfect_matches_downstream_junction, BBMAP_and_STAR, BBMAP_only, STAR_only, BBMAP_intron_not_in_STAR, flag_99, flag_147, flag_83, flag_163 = junction_dict[junction_ID]
            
            if common_gene_name != 'HAC1': ## don't adjust HAC1 splice junctions
            ## HAC1 is the only intron in an ORF known to splice with a different splicing mechanism - IRE1 cleavage followed by tRNA ligase
                adjusted_start, adjusted_end, coordinate_adjustment, ambiguous_bases_around_junction = check_for_ambiguous_junction_sequence(chromosome, intron_start, intron_end, annotated_strand, R64_genome)
            else:
                adjusted_start, adjusted_end, coordinate_adjustment, ambiguous_bases_around_junction = intron_start, intron_end, 0, 0
            alt_5SS, alt_3SS = False, False
            if annotated_strand == '+':
                if adjusted_start != annotated_start:
                    alt_5SS = True
                if adjusted_end != annotated_end:
                    alt_3SS = True
#                if alt_5SS and annotated_start == intron_start and adjusted_start != intron_start:
#                    print('annotated 5SS adjusted', common_gene_name, annotated_chromosome, annotated_strand, annotated_start, annotated_end)
#                    print(intron_start, intron_end, 'adjusted to', adjusted_start, adjusted_end, '\n')
#                if alt_3SS and annotated_end == intron_end and adjusted_end != intron_end:
#                    print('annotated 3SS adjusted', common_gene_name, annotated_chromosome, annotated_strand, annotated_start, annotated_end)
#                    print(intron_start, intron_end, 'adjusted to', adjusted_start, adjusted_end, '\n')
            elif annotated_strand == '-':
                if adjusted_start != annotated_start:
                    alt_3SS = True
                if adjusted_end != annotated_end:
                    alt_5SS = True
#                if alt_5SS and annotated_end == intron_end and adjusted_end != intron_end:
#                    print('annotated 5SS adjusted', common_gene_name, annotated_chromosome, annotated_strand, annotated_start, annotated_end)
#                    print(intron_start, intron_end, 'adjusted to', adjusted_start, adjusted_end, '\n')
#                if alt_3SS and annotated_start == intron_start  and adjusted_start  != intron_start :
#                    print('annotated 3SS adjusted', common_gene_name, annotated_chromosome, annotated_strand, annotated_start, annotated_end)
#                    print(intron_start, intron_end, 'adjusted to', adjusted_start, adjusted_end, '\n')
            
            BP_info = get_closest_upstream_branchpoint(annotated_chromosome, adjusted_start, adjusted_end, annotated_strand, R64_genome, systematic_gene_name, common_gene_name)
            annotated_all_flanking_seq_info = get_sequence_flanking_SS(chromosome, adjusted_start, adjusted_end, annotated_strand, R64_genome)
            junction_ID_name = '_'.join( [str(e) for e in [chromosome, adjusted_start, adjusted_end, strand ]])  
            if common_gene_name == 'GCR1':
                print(junction_ID_name)
            annotated_and_alternative_strands_agree = False
            if strand == annotated_strand:
                annotated_and_alternative_strands_agree = True
            ## common info (will not change depending on the sample data) 
            ## sample_specific info (will change depending on the sample data) for the annotated intron   
            ## for the annotated intron
            common_annotated_info = [annotated_junction_name, annotated_strand, systematic_gene_name, common_gene_name, RPG]           
            
#            if common_gene_name == 'RPL31A':
#                print(common_annotated_info, sample_specific_annotated_info)
            common_alternative_info = [junction_ID_name, chromosome, adjusted_start, adjusted_end, ambiguous_bases_around_junction,  alt_5SS, alt_3SS]
#            if junction_ID_name == 'chrII_407028_407123_-':
#                print(common_alternative_info, chromosome, start, end, strand, reads_without_mismatch_near_junction, sorted_coordinate_cigar_combos)
            coordinate_adjustment = str(coordinate_adjustment)
            sample_specific_alternative_info = [coordinate_adjustment, strand, annotated_and_alternative_strands_agree, reads_without_mismatch_near_junction, sorted_coordinate_cigar_combos, sorted_mismatches_dict,  reads_with_perfect_gapped_alignments, reads_with_gapped_alignments, unique_junction_counts, max_overhang_upstream_junction, max_overhang_perfect_matches_upstream_junction, max_overhang_downstream_junction, max_overhang_perfect_matches_downstream_junction, BBMAP_and_STAR, BBMAP_only, STAR_only, BBMAP_intron_not_in_STAR, flag_99, flag_147, flag_83, flag_163]
            
            ## [ max_overhang_upstream_junction, max_overhang_downstream_junction, ambiguous_bases_around_junction, junction_ID, strand, strands_agree, coordinate_adjustment,  reads_with_mismatch_near_junction, one_read_with_no_mismatches_near_junction, perfect_gapped_alignments, num_uniquely_mapped_reads_across_junction, unique_clusters_with_same_junction, ] +
            if junction_ID_name not in unique_junctions_after_site_adjustment:
                unique_junctions_after_site_adjustment[junction_ID_name] =  [common_annotated_info, annotated_intron_BP_info, annotated_intron_all_flanking_seq_info, common_alternative_info, BP_info, annotated_all_flanking_seq_info, sample_specific_alternative_info]
            else:
                # print('clashing junctions after site adjustment')
                previous_sample_specific_alternative_info = unique_junctions_after_site_adjustment[junction_ID_name][-1]
                updated_sample_specific_alternative_info = combine_clashing_junction_info(sample_specific_alternative_info, previous_sample_specific_alternative_info)
                unique_junctions_after_site_adjustment[junction_ID_name][-1] = updated_sample_specific_alternative_info
        ## SE will be calculated after gathering all reads (spliced and unspliced) for each intron
        
        ## update actual spliced junction counts from STAR to be the actual processed counts from STAR/BBMAP pipeline
        for junction_ID_name in unique_junctions_after_site_adjustment:
           # print(common_alternative_info[0], common_annotated_info[0] )  
            common_annotated_info, annotated_intron_BP_info, annotated_intron_all_flanking_seq_info, common_alternative_info, BP_info, annotated_all_flanking_seq_info, sample_specific_alternative_info = unique_junctions_after_site_adjustment[junction_ID_name]
            coordinate_adjustment, strand, annotated_and_alternative_strands_agree, reads_without_mismatch_near_junction, sorted_coordinate_cigar_combos, sorted_mismatches_dict,  reads_with_perfect_gapped_alignments, reads_with_gapped_alignments, unique_junction_counts, max_overhang_upstream_junction, max_overhang_perfect_matches_upstream_junction, max_overhang_downstream_junction, max_overhang_perfect_matches_downstream_junction, BBMAP_and_STAR, BBMAP_only, STAR_only, BBMAP_intron_not_in_STAR, flag_99, flag_147, flag_83, flag_163 = sample_specific_alternative_info
            annotated_junction_name, annotated_strand, systematic_gene_name, common_gene_name, RPG = common_annotated_info
            if common_gene_name == 'GCR1':
                print('GCR1 common_alternative_info', common_alternative_info[0], common_annotated_info[0])
            if common_alternative_info[0] == common_annotated_info[0]:
                
                annotated_introns_info[ annotated_junction_ID ] = systematic_gene_name, common_gene_name, unique_junction_counts, intron_reads5, intron_reads3, annotated_SE5, annotated_SE3, intron_reads
                if common_gene_name in gene_names_to_check:
                    print(common_gene_name, annotated_junction_ID, unique_junction_counts, annotated_introns_info[ annotated_junction_ID ] )
        for unique_junction_after_site_adjustments in unique_junctions_after_site_adjustment:  
            common_annotated_info, annotated_intron_BP_info, annotated_intron_all_flanking_seq_info, common_alternative_info, BP_info, annotated_all_flanking_seq_info, sample_specific_alternative_info = unique_junctions_after_site_adjustment[junction_ID_name]
            coordinate_adjustment, strand, annotated_and_alternative_strands_agree, reads_without_mismatch_near_junction, sorted_coordinate_cigar_combos, sorted_mismatches_dict,  reads_with_perfect_gapped_alignments, reads_with_gapped_alignments, unique_junction_counts, max_overhang_upstream_junction, max_overhang_perfect_matches_upstream_junction, max_overhang_downstream_junction, max_overhang_perfect_matches_downstream_junction, BBMAP_and_STAR, BBMAP_only, STAR_only, BBMAP_intron_not_in_STAR, flag_99, flag_147, flag_83, flag_163 = sample_specific_alternative_info
            junction_ID_name, chromosome, adjusted_start, adjusted_end, ambiguous_bases_around_junction,  alt_5SS, alt_3SS = common_alternative_info
            systematic_gene_name, common_gene_name, annotated_unique_junction_reads, intron_reads5, intron_reads3, annotated_SE5, annotated_SE3, intron_reads = annotated_introns_info[ annotated_junction_ID ]       
            total_spliced_and_unspliced_reads_associated_with_annotated_intron = annotated_unique_junction_reads + intron_reads
            if total_spliced_and_unspliced_reads_associated_with_annotated_intron == 0:
                annotated_SE = 'NA'
            else:
                annotated_SE = annotated_unique_junction_reads / float(total_spliced_and_unspliced_reads_associated_with_annotated_intron)    
            
            sample_specific_annotated_info = [annotated_unique_junction_reads, intron_reads, total_spliced_and_unspliced_reads_associated_with_annotated_intron, annotated_SE]
            total_spliced_and_unspliced_reads_associated_with_annotated_intron += unique_junction_counts
            
            common_annotated_info, annotated_intron_BP_info, annotated_intron_all_flanking_seq_info, common_alternative_info, BP_info, annotated_all_flanking_seq_info, sample_specific_alternative_info = unique_junctions_after_site_adjustment[unique_junction_after_site_adjustments]
            annotated_unique_junction_reads, intron_reads, total_spliced_and_unspliced_reads_associated_with_annotated_intron, annotated_SE =  sample_specific_annotated_info           
            updated_coordinate_adjustment, updated_strand, updated_annotated_and_alternative_strands_agree, updated_reads_without_mismatch_near_junction, udpated_sorted_coordinate_cigar_combos, udpated_sorted_mismatches_dict, updated_reads_with_perfect_gapped_alignments, updated_reads_with_gapped_alignments, updated_unique_junction_counts, updated_max_overhang_upstream_junction, updated_max_overhang_perfect_matches_upstream_junction, updated_max_overhang_downstream_junction, updated_max_overhang_perfect_matches_downstream_junction, updated_BBMAP_and_STAR, updated_BBMAP_only, updated_STAR_only, updated_BBMAP_intron_not_in_STAR, updated_flag_99, updated_flag_147, updated_flag_83, updated_flag_163 = sample_specific_alternative_info
            junction_ID_name, chromosome, adjusted_start, adjusted_end, ambiguous_bases_around_junction,  alt_5SS, alt_3SS  = common_alternative_info
            annotated_junction_name, annotated_strand, systematic_gene_name, common_gene_name, RPG = common_annotated_info
    
            if common_gene_name in gene_names_to_check:
#            if chromosome ==  'chrXIV' and 	adjusted_start == 545291 and adjusted_end == 545369:    
#            if common_gene_name == 'YIP3':
                print(common_gene_name, junction_ID_name, annotated_unique_junction_reads, total_spliced_and_unspliced_reads_associated_with_annotated_intron, udpated_sorted_coordinate_cigar_combos)
            ## calculate SE for alternative junction    
            if total_spliced_and_unspliced_reads_associated_with_annotated_intron == 0:
                SE = 'NA'
            else:
                SE = updated_unique_junction_counts / float(total_spliced_and_unspliced_reads_associated_with_annotated_intron)  
            ## calculate FAnS for alternative junction    
            if annotated_unique_junction_reads == 0:
                fraction_of_annotated_splicing = 'NA'
            else:
                fraction_of_annotated_splicing = updated_unique_junction_counts / float( annotated_unique_junction_reads )
                
            udpated_sorted_coordinate_cigar_combos_list = sorted(udpated_sorted_coordinate_cigar_combos.items(), key=operator.itemgetter(1), reverse = True)
            udpated_sorted_mismatches_dict_list = sorted(udpated_sorted_mismatches_dict.items(), key=operator.itemgetter(1), reverse = True)
            
            complete_sample_specific_alternative_info = [updated_coordinate_adjustment, updated_strand, SE, fraction_of_annotated_splicing, updated_reads_without_mismatch_near_junction, udpated_sorted_coordinate_cigar_combos_list, udpated_sorted_mismatches_dict_list, updated_reads_with_perfect_gapped_alignments, updated_reads_with_gapped_alignments, updated_unique_junction_counts, updated_max_overhang_upstream_junction, updated_max_overhang_perfect_matches_upstream_junction, updated_max_overhang_downstream_junction, updated_max_overhang_perfect_matches_downstream_junction, updated_BBMAP_and_STAR, updated_BBMAP_only, updated_STAR_only, updated_BBMAP_intron_not_in_STAR, updated_flag_99, updated_flag_147, updated_flag_83, updated_flag_163] 
            output = common_annotated_info + annotated_intron_BP_info + annotated_intron_all_flanking_seq_info + common_alternative_info + BP_info + annotated_all_flanking_seq_info + sample_specific_annotated_info + complete_sample_specific_alternative_info
            outfile.write( '\t'.join([str(e) for e in output]) + '\n' )
    outfile.close()

try:
    R64_genome
except:
    print('loading R64 genome')
    R64_genome = load_genome(genome_file)
    
#output_outfile = open(DIR + 'test_splice_junction_processing_output.txt', 'w')
#chromosome, query_start, query_end, strand = 'chrII_407029_407124_-'.split('_')
#query_start, query_end = int(query_start), int( query_end )
#print( check_for_ambiguous_junction_sequence(chromosome, query_start, query_end, strand, R64_genome) )
#output_outfile.close()
    
# FOR TESTING ON ONE SAMPLE
#for sample in 'upf1',: # 
#    annotated_intron_filename = DIR + sample + '_splicing_efficiency.txt'
#    annotated_introns_info = load_annotated_splicing_efficiency(annotated_intron_filename)
#    splice_junction_filename = DIR + sample + '_combined_bbmap_STAR_reformatted_cigar_SJ.out.tab'    
#    sample_junction_counts, total_uniquely_mapped_gapped_reads = load_splice_junctions(splice_junction_filename)
#    annotated_and_alternative_outfilename = DIR + sample + '_combined_bbmap_STAR_splice_junctions_annotated_and_alternative_splicing_stats_test.txt'
#    output_outfile = open(DIR + sample + '_splice_junction_processing_output_test.txt', 'w')
#    splice_site_analysis(annotated_introns_info, sample_junction_counts, annotated_and_alternative_outfilename)
#    output_outfile.close()
#    
for sample in SAMPLES: 
    output_outfile = open(DIR + sample + '_splice_junction_processing_output.txt', 'w')
    annotated_intron_filename = DIR + sample + '_splicing_efficiency.txt'
    annotated_introns_info = load_annotated_splicing_efficiency(annotated_intron_filename)
    splice_junction_filename = DIR + sample + '_combined_bbmap_STAR_reformatted_cigar_SJ.out.tab'    
    sample_junction_counts, total_uniquely_mapped_gapped_reads = load_splice_junctions(splice_junction_filename)
    annotated_and_alternative_outfilename = DIR + sample + '_combined_bbmap_STAR_splice_junctions_annotated_and_alternative_splicing_stats.txt'
    splice_site_analysis(annotated_introns_info, sample_junction_counts, annotated_and_alternative_outfilename)
    output_outfile.close()