# -*- coding: utf-8 -*-
"""
Created on Sun Feb 25 21:38:34 2018

@author: kevinroy

command line usage: python compare_splice_junctions_from_multiple_aligner_SAM.py GENOME_FASTA GFF_ANNOTATIONS ALIGNMENTS_ALIGNERS JUNCTION_STATISTICS_OUTFILENAME INTEGRATED_ALIGNMENT_SUMMARY_OUTFILENAME INTEGRATED_ALIGNMENT_SAM_OUTFILENAME  FILTERED_SAM_DIR SAMPLE_PREFIX

takes arbitrary number of SAM files containing alignments

***NOTE ENTIRETY OF ALL SAM FILES ARE LOADED INTO MEMORY***

arguments:
-each SAM file and aligner need to be separated by a comma, 
semicolon separating different alignment/aligner name combos as
SAM file_1,aligner_1;SAM file_2,aligner_2; etc.

requirements:
-SAM must have same read ID format for each aligner

-CIGAR format must be N for intron gap, M for match, S for soft-clip, 
I for insertion, D for deletion, and X for mismatch

what this script does:
-selects the optimal alignment based on user-defined criteria for penalties 
for each CIGAR operation, with the option of breaking ties for annotated 
junction alignment
-indicates which aligners agree on the selected junction for each read
-writes an "integrated" junction file with summary statistics on agreement between aligners for each junction

SRR5582776_bbmap_default.sam,BBMap_default
SRR5582776_bbmap_all_reads.sam,BBMap_tuned
SRR5582776_default_STAR_with_annotation_reformatted_cigar.sam,STAR_default_with_annotation
SRR5582776_default_STAR_no_annotation_reformatted_cigar.sam,STAR_default_without_annotation
SRR5582776_STAR_with_annotation_all_reads_reformatted_cigar.sam,STAR_noncanonical_with_annotation
SRR5582776_STAR_no_annotation_all_reads_reformatted_cigar.sam,STAR_noncanonical_without_annotation
SRR5582776_HISAT2_default_annotated_reformatted_cigar.sam,HISAT2_default_with_annotation
SRR5582776_HISAT2_all_reads_reformatted_cigar.sam,HISAT2_noncanonical_without_annotation
SRR5582776_ContextMap2_noncanonical_reformatted_cigar.sam,HISAT2_noncanonical_without_annotation
SRR5582776_ContextMap2_default_annotated_reformatted_cigar.sam,HISAT2_default_with_annotation
"""

import operator, random
import pysam

#from sys import argv

#script_placeholder, GENOME_FASTA, GFF_ANNOTATIONS, ALIGNMENTS_ALIGNERS, JUNCTION_STATISTICS_OUTFILENAME, INTEGRATED_ALIGNMENT_SUMMARY_OUTFILENAME, FILTERED_SAM_DIR, SAMPLE_PREFIX = argv
#print(GENOME_FASTA, GFF_ANNOTATIONS, ALIGNMENTS_ALIGNERS, JUNCTION_STATISTICS_OUTFILENAME, INTEGRATED_ALIGNMENT_SUMMARY_OUTFILENAME, FILTERED_SAM_DIR, SAMPLE_PREFIX)

GENOME_FASTA = '/Users/kevinroy/Google_Drive/yeast_genome_references/S288C_reference_sequence_R64-2-1_20150113_reformatted_chromosome_names.fasta'
GFF_ANNOTATIONS = '/Users/kevinroy/Google_Drive/yeast_genome_references/saccharomyces_cerevisiae_annotations.gff'
ALIGNMENTS_ALIGNERS = '/Volumes/SP_PHD/PRJNA387451/original_STAR_bam_reprocessed/COMPASS_unannotated_junctions/SRR5582776_bbmap_default.sam,BBMap_default~/Volumes/SP_PHD/PRJNA387451/original_STAR_bam_reprocessed/COMPASS_unannotated_junctions/SRR5582776_bbmap_all_reads.sam,BBMap_tuned~/Volumes/SP_PHD/PRJNA387451/original_STAR_bam_reprocessed/COMPASS_unannotated_junctions/SRR5582776_default_STAR_with_annotation_reformatted_cigar.sam,STAR_default_with_annotation~/Volumes/SP_PHD/PRJNA387451/original_STAR_bam_reprocessed/COMPASS_unannotated_junctions/SRR5582776_default_STAR_no_annotation_reformatted_cigar.sam,STAR_default_without_annotation~/Volumes/SP_PHD/PRJNA387451/original_STAR_bam_reprocessed/COMPASS_unannotated_junctions/SRR5582776_STAR_with_annotation_all_reads_reformatted_cigar.sam,STAR_noncanonical_with_annotation~/Volumes/SP_PHD/PRJNA387451/original_STAR_bam_reprocessed/COMPASS_unannotated_junctions/SRR5582776_STAR_no_annotation_all_reads_reformatted_cigar.sam,STAR_noncanonical_without_annotation~/Volumes/SP_PHD/PRJNA387451/original_STAR_bam_reprocessed/COMPASS_unannotated_junctions/SRR5582776_HISAT2_default_annotated_reformatted_cigar.sam,HISAT2_default_with_annotation~/Volumes/SP_PHD/PRJNA387451/original_STAR_bam_reprocessed/COMPASS_unannotated_junctions/SRR5582776_HISAT2_all_reads_reformatted_cigar.sam,HISAT2_noncanonical_without_annotation~/Volumes/SP_PHD/PRJNA387451/original_STAR_bam_reprocessed/COMPASS_unannotated_junctions/SRR5582776_ContextMap2_noncanonical_reformatted_cigar.sam,ContextMap2_noncanonical_without_annotation~/Volumes/SP_PHD/PRJNA387451/original_STAR_bam_reprocessed/COMPASS_unannotated_junctions/SRR5582776_ContextMap2_default_annotated_reformatted_cigar.sam,ContextMap2_default_with_annotation'
JUNCTION_STATISTICS_OUTFILENAME = '/Volumes/SP_PHD/PRJNA387451/original_STAR_bam_reprocessed/COMPASS_unannotated_junctions/SRR5582776_junction_statistics.txt'
INTEGRATED_ALIGNMENT_SUMMARY_OUTFILENAME = '/Volumes/SP_PHD/PRJNA387451/original_STAR_bam_reprocessed/COMPASS_unannotated_junctions/SRR5582776_integrated_alignments.txt'
FILTERED_SAM_DIR = '/Volumes/SP_PHD/PRJNA387451/original_STAR_bam_reprocessed/COMPASS_unannotated_junctions/'
SAMPLE_PREFIX = 'SRR5582776'

#SAMPLE_NAME = 'SRR5582776'
#GENOME_FASTA = '/Users/kevinroy/Google_Drive/yeast_genome_references/S288C_reference_sequence_R64-2-1_20150113_reformatted_chromosome_names.fasta'
#GFF_ANNOTATIONS = '/Users/kevinroy/Google_Drive/yeast_genome_references/saccharomyces_cerevisiae_annotations.gff'
#DIR = '/Volumes/MyPassportforMac_5TB/processed_sequencing_data/fast_slow_Pol_II_splicing/'
#ALIGNMENTS_ALIGNERS = DIR + SAMPLE_NAME + '_STAR_with_annotation_all_reads_reformatted_cigar.sam,STAR_with_annotation;' + DIR + SAMPLE_NAME + '_STAR_no_annotation_all_reads_reformatted_cigar.sam,STAR;' + DIR + SAMPLE_NAME + '_bbmap_all_reads.sam,BBMap;'  + DIR + SAMPLE_NAME + '_HISAT2_all_reads_reformatted_cigar.sam,HISAT2'
#JUNCTION_STATISTICS_OUTFILENAME = DIR + SAMPLE_NAME + '_junction_statistics.txt'
#INTEGRATED_ALIGNMENT_SUMMARY_OUTFILENAME = DIR + SAMPLE_NAME + '_integrated_alignments.txt'
#INTEGRATED_ALIGNMENT_SAM_OUTFILENAME = DIR + SAMPLE_NAME + '_integrated_alignments.sam'
#FILTERED_SAM_DIR = DIR
#SAMPLE_PREFIX = SAMPLE_NAME

MISMATCH_DIST_FROM_JUNCTION_DISALLOWED = 10 ## i.e.  how many bp from the junction must be mismatch free (on either side)
MINIMUM_INTRON_LENGTH = 20
MAXIMUM_INTRON_LENGTH = 2000
CIGAR_OPERATION_TO_PENALTY = {'N':1, 'X':1, 'I':1, 'D':1, 'S':1, 'M':1, '=':0  } # N is only a penalty if less than MINIMUM_INTRON_LENGTH or greater than MAXIMUM_INTRON_LENGTH
# SHORT INTRON PENALTY IS DELETION LENGTH, long introns get higher penalty
LONG_INTRON_PENALTY = 10

def get_gene_id_from_gff_annotation(string):
    '''
    input: annotation field in gff format
    output: gene_id
    '''
    if 'ID=' not in string:
        return None
    dist = len('ID=')
    ID_found = False
    for idx in range(len(string)):
        if string[idx:idx+dist] == 'ID=':
            ID_found = True
            id_start = idx+dist
            break
    if not ID_found:
        return None
    id_end = None
    for idx in range(id_start, len(string)):
        if string[idx] == ';':
            id_end = idx
            break
    if id_end == None:
        return None
    return string[id_start:id_end]

def get_parent_from_gff_annotation(string):
    '''
    input: annotation field in gff format
    output: Parent_id
    '''
    dist = len('Parent=')
    ID_found = False
    for idx in range(len(string)):
        if string[idx:idx+dist] == 'Parent=':
            ID_found = True
            id_start = idx+dist
            break
    if not ID_found:
        return None
    for idx in range(id_start, len(string)):
        if string[idx] == ';':
            id_end = idx
            break
    if ID_found:
        return string[id_start:id_end]
    else:
        return None

def load_introns(filename):
    '''
    input: gff file with intron annotations
    output: dictionary of chromosomes of start and stop coordinate tuples of strands
    '''
    introns = {}
    splice_sites = {}
    with open(filename,'r') as annotations:
        for line in annotations:
            info = line.strip().split('\t')
            chromosome, region, sequence_type, start, end, period, strand, frame, annotation = info
            if chromosome != 'chrMito':  ## exlude mitochondrial introns
                start, end = int(start), int(end)
                if sequence_type == 'intron' or sequence_type == 'five_prime_UTR_intron':
                    gene = get_gene_id_from_gff_annotation(annotation)
                    if gene == None:
                        gene = get_parent_from_gff_annotation(annotation)
                    if gene != None:
                        gene = gene.replace('_mRNA', '')
                    if chromosome not in introns:
                        introns[chromosome] = {}
                        splice_sites[chromosome] = set([])
                    if 't' not in gene:  ## exclude tRNA introns & optionally and gene != 'HAC1'
                        introns[chromosome][(start,end)] = (strand, gene)
                        splice_sites[chromosome].add(start)
                        splice_sites[chromosome].add(end)
    return introns, splice_sites

try:
    annotated_introns
except:
    print('loading annotated_introns')
    annotated_introns, splice_sites = load_introns(GFF_ANNOTATIONS)
print('annotated_introns loaded')

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
    print(GENOME_FASTA + " genome file loaded")
    return genome_dict

#
#try:
#    R64_genome
#except:
#    print('loading R64 genome')
#    R64_genome = load_genome(GENOME_FASTA)
#print('R64 genome loaded')
# 
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

def toBinary(n):
    '''
    input: binary SAM flag
    output: 0/1 representation of flag
    '''
    return ''.join(str(1 & int(n) >> i) for i in range(12)[::-1])

def process_sam(sam_fn):
    '''
    input: SAM filename (from a filtered sam file, where all read_IDs have one of the paired-end reads with a gapped alignment)
    output: dictionary of read_ID locations of read pairs, with the following tuple as a value:
    strand (string), whether or not read is mapped perfectly, and contiguously, without an intron or soft-clipping (boolean), 
    the entire line as a list (list of strings), and the flag (string in 0/1 represntation)
    '''
    read_ID_to_read_num_to_info = {}
    SAM_header = ''
    with open(sam_fn) as infile:
        print('looking for read_ID locations with gapped alignments in sam... ')
        line_count = 0
        for line in infile:
            line_count += 1
            if line[0] != "@":         
                info = line.split('\t')
                ## coordinates in SAM are 1-based!
#                info[ 5 ] = info[ 5 ].replace('=', 'M') ## BMAP uses equals sign "=" for matching bases, HISAT2 and STAR use "M", all other CIGAR operations are the same
                read_ID, flag, chromosome, coordinate, mapping_quality, cigar, equals_sign, paired_end_location, observed_template_length, sequence, quality_scores = info[:11]              
                optional_fields = info[11:]
                info = read_ID, flag, chromosome, coordinate, mapping_quality, cigar, equals_sign, paired_end_location, observed_template_length, '*', '*'
                perfect_cigar = True
                flag_as_string = flag
                if flag in ('99', '147', '83', '163'):  ## only allow properly mapped and paired reads for analysis
                    flag = toBinary(flag)
                    if flag[-8:-6] == '10':
                        read_num = '2'
                    elif flag[-8:-6] == '01':
                        read_num = '1'
                    else:
                        read_num = None
                        print('inconsistent read pair')
                        print(flag)
                        print(line)
                    ## 147 or 99
                    if flag[-8:] == '10010011' or flag[-8:] == '01100011':
                        RNA_strand = '-'
                    ## 163 or 83
                    elif flag[-8:] == '10100011' or flag[-8:] == '01010011':
                        RNA_strand = '+'
                    if 'N' in cigar or 'D' in cigar or 'I' in cigar or 'X' in cigar or 'S' in cigar or 'M' in cigar:  
                        perfect_cigar = False
                    if read_ID not in read_ID_to_read_num_to_info:
                        read_ID_to_read_num_to_info[read_ID] = {}
                    if read_num in read_ID_to_read_num_to_info[read_ID]:
                        print('read number clashing')
                        print('flag_as_string', flag_as_string)           
                    read_ID_to_read_num_to_info[read_ID][read_num] = RNA_strand, perfect_cigar, info, flag_as_string
            elif line[:3] == '@SQ':
                SAM_header += line
            if line_count % 200000 == 0:
                print('processed', line_count, 'reads')
#                break
    print('processed', line_count, 'reads')
    infile.close()
    return read_ID_to_read_num_to_info, SAM_header

def parse_cigar(cigar):
    '''
    input: CIGAR string from SAM alignment
    output: list of 2 item tuples: mapped segment length, CIGAR operation
    '''
    total_cigar_operation_bp = {'N':0, 'X':0, 'I':0, 'D':0, 'S':0, '=':0, 'M':0 }
    segment_length = ''
    parsed_cigar = []
    for char in cigar:
        if char in '0123456789':
            segment_length += char
        else:
            operation = char
            mapped_segment_length = int(segment_length)
            segment_length = ''
            parsed_cigar.append( (mapped_segment_length, operation) )
            total_cigar_operation_bp[operation] += mapped_segment_length
    return parsed_cigar, total_cigar_operation_bp

def filter_SAM_by_multiple_aligner_agreement(aligner_to_processed_alignments, FILTERED_SAM_DIR, SAMPLE_PREFIX):
    complete_agreement_outfilename = FILTERED_SAM_DIR + '/' + SAMPLE_PREFIX + '_all_aligners_agree.sam'
    complete_agreement_outfile = open(complete_agreement_outfilename, 'w')
    read_IDs = set([])
    incomplete_agreement_aligner_to_outfiles = {}
    reads_with_alignment_disagreement = 0
    reads_with_unanimous_aligner_agreement = 0
    reads_with_no_alignment_in_at_least_one_aligner = 0
    read_IDs_with_incomplete_agreement = set([])
#    read_IDs_to_process = 100
    for aligner in aligner_to_processed_alignments:
        read_ID_to_read_num_to_info, SAM_header =  aligner_to_processed_alignments[aligner]
        for read_ID in read_ID_to_read_num_to_info:
            read_IDs.add(read_ID)
        incomplete_agreement_aligner_to_outfiles[aligner] = open(FILTERED_SAM_DIR + '/' + SAMPLE_PREFIX + '_incomplete_agreement_' + aligner  + '.sam', 'w')
        incomplete_agreement_aligner_to_outfiles[aligner].write(SAM_header)
    ## process each paired-end read one at a time and compare with all aligners  
    read_pairs_processed = 0
    for read_ID in read_IDs:
        read_pairs_processed += 1
        if read_pairs_processed % 200000 == 0:
            print(read_pairs_processed, 'read_pairs_processed' )
#            break
        alignments = {}
        alignments['1'] = []
        alignments['2'] = []
        alignment_missing_in_one_aligner = False
        alignment_disagreement = False

        for aligner in aligner_to_processed_alignments:
            read_ID_to_read_num_to_info, SAM_header =  aligner_to_processed_alignments[aligner]
            if read_ID in read_ID_to_read_num_to_info:  
                for read_num in read_ID_to_read_num_to_info[read_ID]:
                    RNA_strand, perfect_cigar, info, flag_as_string = read_ID_to_read_num_to_info[read_ID][read_num]
                    read_ID, flag, chromosome, coordinate, mapping_quality, cigar, equals_sign, paired_end_location, observed_template_length, sequence, quality_scores = info[:11]
                    output = '\t'.join( info[:11] ) + '\n'
                    
                    if alignments[read_num] == []:
                        alignments[read_num] = (chromosome, coordinate, cigar)
                    else:
                        if (chromosome, coordinate, cigar) != alignments[read_num]:
                            alignment_disagreement = True
            else:
                alignment_missing_in_one_aligner = True
        if alignment_disagreement:
            reads_with_alignment_disagreement += 1
        elif alignment_missing_in_one_aligner:
            reads_with_no_alignment_in_at_least_one_aligner += 1
        else:
            complete_agreement_outfile.write( output )
            reads_with_unanimous_aligner_agreement += 1
        if alignment_missing_in_one_aligner or alignment_disagreement:
            read_IDs_with_incomplete_agreement.add(read_ID)
            for aligner in aligner_to_processed_alignments:
                read_ID_to_read_num_to_info, SAM_header =  aligner_to_processed_alignments[aligner]
                if read_ID in read_ID_to_read_num_to_info:  
                    for read_num in read_ID_to_read_num_to_info[read_ID]:
                        RNA_strand, perfect_cigar, info, flag_as_string = read_ID_to_read_num_to_info[read_ID][read_num]
                        read_ID, flag, chromosome, coordinate, mapping_quality, cigar, equals_sign, paired_end_location, observed_template_length, sequence, quality_scores = info[:11]
                        output = '\t'.join( info[:11] ) + '\n'
                        incomplete_agreement_aligner_to_outfiles[aligner].write(output)
    complete_agreement_outfile.close()
    for aligner in aligner_to_processed_alignments:
        incomplete_agreement_aligner_to_outfiles[aligner].close()
    print(SAMPLE_PREFIX + ': \n' + str(reads_with_unanimous_aligner_agreement) + 'reads_with_unanimous_aligner_agreement\n' + str(reads_with_alignment_disagreement) + 'reads_with_alignment_disagreement\n' + str(reads_with_no_alignment_in_at_least_one_aligner) + 'reads_with_no_alignment_in_at_least_one_aligner\n' )
    return read_IDs, read_IDs_with_incomplete_agreement

def find_optimal_alignments(all_aligners, aligner_to_processed_alignments, read_IDs, read_IDs_with_incomplete_agreement, INTEGRATED_ALIGNMENT_SUMMARY_OUTFILENAME):
    '''
    input: dictionaries with aligner name as keys and processed SAM as values, and an outfilename
    SAM dictionaries have the following format: read_ID_to_read_num_to_bbmap_info[read_ID][read_num] = RNA_strand, perfect_cigar, info, flag_as_string
    compares alignments for each read_ID and read_num
    chooses the gapped alignment with fewer mismatches
    output: a dictionary containing 
    '''  
    integrated_alignments = {}
    integrated_alignment_info_outfile = open(INTEGRATED_ALIGNMENT_SUMMARY_OUTFILENAME, 'w')
    header = '\t'.join('read_pair_ID, aligner, aligner_had_this_read_pair_mapped, aligner_has_best_score, alignment_score_for_read_pair, cigar_for_read_pair, intron_gap, annotated_junction, least_alignment_score, complete_agreement, ambiguous_by_same_score_with_different_cigar, total_matching_bases, total_X, total_S, total_I, total_D, total_N, total_M, total_cigar_operations'.split(', '))  + '\n'
 
 # read_pair_ID, aligner, aligner_had_this_read_pair_mapped, aligner_has_best_score, alignment_score_for_read_pair, cigar_for_read_pair, read_pair_contains_annotated_junction, least_alignment_score, complete_agreement, ambiguous_by_same_score_with_different_cigar, total_matching_bases, total_X, total_S, total_I, total_D, total_N, total_cigar_operations'  + '\n'
     
    integrated_alignment_info_outfile.write(header)       
    ## process each paired-end read one at a time and compare with all aligners  
    read_pairs_processed = 0
    for read_ID in read_IDs:
        read_pairs_processed += 1
        if read_pairs_processed % 100000 == 0:
            print(read_pairs_processed, 'read_pairs_processed' )
#            break
        aligner_to_processed_cigar = {}
        complete_agreement = True
        aligners_with_this_read_ID = []
        for aligner in all_aligners:
            read_ID_to_read_num_to_info, SAM_header =  aligner_to_processed_alignments[aligner]
            if read_ID in read_ID_to_read_num_to_info:  
                aligners_with_this_read_ID.append(aligner)
                aligner_to_processed_cigar[aligner] = {}
                for read_num in read_ID_to_read_num_to_info[read_ID]:
                    RNA_strand, perfect_cigar, info, flag_as_string = read_ID_to_read_num_to_info[read_ID][read_num]
                    read_ID, flag, chromosome, coordinate, mapping_quality, cigar, equals_sign, paired_end_location, observed_template_length, sequence, quality_scores = info[:11]          
                    alignment = (chromosome, coordinate, cigar)                   
                    coordinate = int(coordinate)
                    parsed_cigar, total_cigar_operation_bp = parse_cigar(cigar)
                    
                    ## parsing the cigar has two main functions: 
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
                    mismatch_near_junction = False
                    intron_found = False
                    annotated_junction = False
                    for cigar_idx in range(len(parsed_cigar)):
                        bp, operation = parsed_cigar[cigar_idx]    
                        if operation == '=':
                            mapped_segment_length += bp
                            coordinate += bp
                        elif operation == 'N':
                            total_gapped_bp += bp
                            mapped_segment_lengths.append(mapped_segment_length)
                            mapped_segment_length = 0
                            if parsed_cigar[cigar_idx-1][1] != '=': 
                                left_perfect_matches = 0
                            else:
                                left_perfect_matches = parsed_cigar[cigar_idx-1][0] 
                            if parsed_cigar[cigar_idx+1][1] != '=' :
                                right_perfect_matches = 0
                            else:
                                right_perfect_matches = parsed_cigar[cigar_idx+1][0]
                            if bp >= MINIMUM_INTRON_LENGTH:
                                if bp > MAXIMUM_INTRON_LENGTH:
                                    alignment_score += LONG_INTRON_PENALTY
                                intron_found = True
                                splice_sites.append(coordinate)
                                coordinate = coordinate + bp
                                splice_sites.append(coordinate - 1)
                                perfect_matches_flanking_splice_sites.append( left_perfect_matches)
                                perfect_matches_flanking_splice_sites.append( right_perfect_matches)
                            else:
                                alignment_score += bp * CIGAR_OPERATION_TO_PENALTY[operation]

                        elif operation in ('X', 'I', 'S', 'D', 'M'):
                            perfect_gapped_alignment = False
                            alignment_score += bp * CIGAR_OPERATION_TO_PENALTY[operation]
                            if ( cigar_idx-1 >= 0 and parsed_cigar[cigar_idx-1][1]  == 'N') or (cigar_idx+1 < len(parsed_cigar) and parsed_cigar[cigar_idx+1][1]  == 'N'):
                                mismatch_near_junction = True
                            elif ( cigar_idx-2 >= 0 and parsed_cigar[cigar_idx-2][1]  == 'N' and parsed_cigar[cigar_idx-1][1]  == '=' and parsed_cigar[cigar_idx-1][0]  <= MISMATCH_DIST_FROM_JUNCTION_DISALLOWED) or (cigar_idx+2 < len(parsed_cigar) and parsed_cigar[cigar_idx+2][1]  == 'N' and parsed_cigar[cigar_idx+1][1]  == '=' and parsed_cigar[cigar_idx+1][0]  <= MISMATCH_DIST_FROM_JUNCTION_DISALLOWED):
                                mismatch_near_junction = True
                        if operation in ('X', 'D', 'M'):
                            coordinate += bp
                            mapped_segment_length += bp
                        
                    mapped_segment_lengths.append(mapped_segment_length)
                    for junction_idx in range(int(len(splice_sites)/2) ) :
                        start, end = splice_sites[junction_idx*2], splice_sites[junction_idx*2 + 1]
                        intron_coords = start, end
                        if intron_coords in annotated_introns[chromosome]:
                            annotated_junction = True
                    aligner_to_processed_cigar[aligner][read_num] = alignment_score, total_cigar_operation_bp, alignment, intron_found, annotated_junction, splice_sites, perfect_matches_flanking_splice_sites, mapped_segment_lengths, perfect_gapped_alignment, mismatch_near_junction, RNA_strand, perfect_cigar, info, flag_as_string
            
        ## read pairs are processed together
        ## pick aligner with least alignment_score, tie broken by match to annotated junction
        ## if there is still a tie, and the identified splice sites differ (for e.g. due to ambiguous junction assignment), randomly pick alignment at this stage  --> keep track of how many junctions this affects
        ## identify how often ties have different alignments       
        
        complete_agreement = True
        if read_ID in read_IDs_with_incomplete_agreement:      
            complete_agreement = False
        least_alignment_score = 1000        
        aligner_to_alignment_score = {}
        aligner_to_alignments = {}
        aligner_to_combined_cigar_operation_bp = {}  
        aligner_to_intron_found = {}
        aligner_to_annotated_junction = {}
        aligner_to_total_cigar_operations = {}
        for aligner in all_aligners:
            ## default values must be set for the following lists
            aligner_to_intron_found[aligner] = False
            aligner_to_alignment_score[aligner] = 'NA'
            aligner_to_alignments[aligner] = []
            aligner_to_annotated_junction[aligner] = False
            aligner_to_combined_cigar_operation_bp[aligner] = {'N':0, 'X':0, 'I':0, 'D':0, 'S':0, '=':0, 'M':0 }
            aligner_to_total_cigar_operations[aligner] = 0
            if aligner in aligners_with_this_read_ID:
                ## process read 1 first
                aligner_to_alignment_score[aligner] = 0
                for read_num in ('1','2'):
                    if read_num not in aligner_to_processed_cigar[aligner]:
                        print(aligner, read_ID, read_num, aligner_to_processed_cigar[aligner])
                    alignment_score, total_cigar_operation_bp, alignment, intron_found, annotated_junction, splice_sites, perfect_matches_flanking_splice_sites, mapped_segment_lengths, perfect_gapped_alignment, mismatch_near_junction, RNA_strand, perfect_cigar, info, flag_as_string = aligner_to_processed_cigar[aligner][read_num]       
                    aligner_to_alignment_score[aligner] += alignment_score
                    for cigar_operation in total_cigar_operation_bp:
                        aligner_to_combined_cigar_operation_bp[aligner][cigar_operation] += total_cigar_operation_bp[cigar_operation]
                        aligner_to_total_cigar_operations[aligner] += 1
                    aligner_to_alignments[aligner] += [alignment]
                    if intron_found:
                        aligner_to_intron_found[aligner] = True
                    if annotated_junction:
                        aligner_to_annotated_junction[aligner] = True     
                if aligner_to_alignment_score[aligner] < least_alignment_score:
                    least_alignment_score = aligner_to_alignment_score[aligner]                        
        best_aligners = []
        for aligner in aligners_with_this_read_ID:
            if aligner_to_alignment_score[aligner] == least_alignment_score:                
                best_aligners.append( aligner )
        ambiguous_by_same_score_with_different_cigar = False      
#        if aligners give equal score, randomly pick one  
        # ties go to annotated aligner
        chosen_aligner = random.choice(best_aligners)
        if len(best_aligners) > 1:
            for idx in range(len(best_aligners)-1):
                aligner = best_aligners[idx]
                next_aligner = best_aligners[idx+1]
                if aligner_to_alignments[aligner] != aligner_to_alignments[next_aligner]:
                    ambiguous_by_same_score_with_different_cigar = True
                    if aligner_to_annotated_junction[aligner]:
                        chosen_aligner = aligner
                    elif aligner_to_annotated_junction[next_aligner]:
                        chosen_aligner = next_aligner
        integrated_alignments[read_ID] = chosen_aligner, least_alignment_score, complete_agreement, ambiguous_by_same_score_with_different_cigar, aligner_to_processed_cigar
        for idx in range(len(all_aligners)):
            aligner = all_aligners[idx]
            aligner_had_this_read_pair_mapped = False
            aligner_has_best_score = False
            if aligner in aligners_with_this_read_ID:
                aligner_had_this_read_pair_mapped = True
            if  aligner in best_aligners:
                aligner_has_best_score = True
            total_matching_bases = aligner_to_combined_cigar_operation_bp[aligner]['=']
            total_X = aligner_to_combined_cigar_operation_bp[aligner]['X']
            total_S = aligner_to_combined_cigar_operation_bp[aligner]['S']
            total_I = aligner_to_combined_cigar_operation_bp[aligner]['I']
            total_D = aligner_to_combined_cigar_operation_bp[aligner]['D']
            total_N = aligner_to_combined_cigar_operation_bp[aligner]['N']
            total_M = aligner_to_combined_cigar_operation_bp[aligner]['M']
            output =  read_ID, aligner, aligner_had_this_read_pair_mapped, aligner_has_best_score, aligner_to_alignment_score[aligner], aligner_to_alignments[aligner], aligner_to_intron_found[aligner], aligner_to_annotated_junction[aligner], least_alignment_score, complete_agreement, ambiguous_by_same_score_with_different_cigar, total_matching_bases, total_X, total_S, total_I, total_D, total_N, total_M, aligner_to_total_cigar_operations[aligner]
            integrated_alignment_info_outfile.write('\t'.join([str(e) for e in output]) + '\n')       
    integrated_alignment_info_outfile.close()
    return integrated_alignments, SAM_header

def process_splice_junctions_from_integrated_alignments(all_aligners, integrated_alignments, SAM_header, FILTERED_SAM_DIR, SAMPLE_PREFIX, JUNCTION_STATISTICS_OUTFILENAME):   
    intron_coords_to_read_IDs = {}    
    integrated_SAM_outfile = open(FILTERED_SAM_DIR + SAMPLE_PREFIX + '_COMPASS_integrated_alignments.sam', 'w')
    integrated_SAM_outfile.write(SAM_header)
    for read_ID in integrated_alignments:
        chosen_aligner, least_alignment_score, complete_agreement, ambiguous_by_same_score_with_different_cigar, aligner_to_processed_cigar = integrated_alignments[read_ID] 
        for read_num in aligner_to_processed_cigar[chosen_aligner]:
            alignment_score, total_cigar_operation_bp, alignment, intron_found, annotated_junction, splice_sites, perfect_matches_flanking_splice_sites, mapped_segment_lengths, perfect_gapped_alignment, mismatch_near_junction, RNA_strand, perfect_cigar, info, flag_as_string = aligner_to_processed_cigar[chosen_aligner][read_num]
            read_ID, flag, chromosome, coordinate, mapping_quality, cigar, equals_sign, paired_end_location, observed_template_length, sequence, quality_scores = info[:11]   
            integrated_SAM_outfile.write('\t'.join([str(e) for e in info[:11] ] ) + '\n')
            introns = []
            for junction_idx in range(int(len(splice_sites)/2) ) :
                start, end = splice_sites[junction_idx*2], splice_sites[junction_idx*2 + 1]
                intron_coords = start, end
                introns.append( intron_coords )
                left_perfect_matches, right_perfect_matches = perfect_matches_flanking_splice_sites[junction_idx*2], perfect_matches_flanking_splice_sites[junction_idx*2 + 1]
                left_segment, right_segment = mapped_segment_lengths[junction_idx], mapped_segment_lengths[junction_idx + 1]

                if chromosome not in intron_coords_to_read_IDs:
                    intron_coords_to_read_IDs[chromosome] = {}
                if (start, end) not in intron_coords_to_read_IDs[chromosome]:
                    intron_coords_to_read_IDs[chromosome][(start, end)] = {}
                if read_ID not in intron_coords_to_read_IDs[chromosome][(start, end)]:
                    intron_coords_to_read_IDs[chromosome][(start, end)][read_ID] = {} 
                intron_coords_to_read_IDs[chromosome][(start, end)][read_ID][read_num] =  mismatch_near_junction, left_segment, right_segment, left_perfect_matches, right_perfect_matches
## chromosome, intron_start, intron_end, strand, intron_motif, annotated, num_uniquely_mapped_reads_across_junction, num_multi_mapped_reads_across_junction, maximum_spliced_alignment_overhang,  = info
    total_candidate_introns = 0
    total_validated_introns = 0
    with open(JUNCTION_STATISTICS_OUTFILENAME, 'w') as outfile:
        header = 'chromosome, start, end, strand, reads_without_mismatch_near_junction, sorted_coordinate_cigar_combos, sorted_mismatches_dict, reads_with_perfect_gapped_alignments, reads_with_gapped_alignments, unique_junction_counts, max_overhang_upstream_junction, max_overhang_perfect_matches_upstream_junction, max_overhang_downstream_junction, max_overhang_perfect_matches_downstream_junction, flag_99, flag_147, flag_83, flag_163'.split(', ')
        aligner_counts_headers = [aligner + '_supporting_junction_counts' for aligner in  all_aligners]      
        header = '\t'.join(header + aligner_counts_headers) + '\n'
        outfile.write(header)
        for chromosome in intron_coords_to_read_IDs:
            print(chromosome)
            for intron in sorted(intron_coords_to_read_IDs[chromosome].keys()):
                (start, end) = intron

                reads_with_perfect_gapped_alignments = 0
                reads_without_mismatch_near_junction = 0
                total_candidate_introns += 1
                reads_with_gapped_alignments = 0
                max_overhang_upstream_junction = 0
                max_overhang_downstream_junction = 0
                max_overhang_perfect_matches_upstream_junction = 0
                max_overhang_perfect_matches_downstream_junction = 0

                read_IDs_counted = set([])
                unique_junction_counts = 0
                aligner_to_supporting_junction_counts = {}
                for aligner in all_aligners:
                    aligner_to_supporting_junction_counts[aligner] = 0
                flag_99, flag_147, flag_83, flag_163 = 0, 0, 0, 0
                mismatches_dict = {}
                coordinate_cigar_combos = {}
                
                for read_ID in intron_coords_to_read_IDs[chromosome][(start, end)]:
                    chosen_aligner, least_alignment_score, complete_agreement, ambiguous_by_same_score_with_different_cigar, aligner_to_processed_cigar = integrated_alignments[read_ID]
                    aligner_to_supporting_junction_counts[ chosen_aligner ] += 1
                    for read_num in intron_coords_to_read_IDs[chromosome][(start, end)][read_ID]:
                        reads_with_gapped_alignments += 1
                        mismatch_near_junction, left_segment, right_segment, left_perfect_matches, right_perfect_matches = intron_coords_to_read_IDs[chromosome][(start, end)][read_ID][read_num]
                        alignment_score, total_cigar_operation_bp, alignment, intron_found, annotated_junction, splice_sites, perfect_matches_flanking_splice_sites, mapped_segment_lengths, perfect_gapped_alignment, mismatch_near_junction, RNA_strand, perfect_cigar, info, flag_as_string = aligner_to_processed_cigar[chosen_aligner][read_num]       
                        read_ID, flag, chromosome, coordinate, mapping_quality, cigar, equals_sign, paired_end_location, observed_template_length, sequence, quality_scores = info[:11]                          
                        if read_ID not in read_IDs_counted:
                            unique_junction_counts += 1  ## prevents the read fragment counting twice if both read 1 and read 2 of the read_ID overlap the same intron
                            read_IDs_counted.add(read_ID)
                        if not mismatch_near_junction:
                            reads_without_mismatch_near_junction += 1    
                        if alignment_score not in mismatches_dict:
                            mismatches_dict[alignment_score] = 1
                        else:
                            mismatches_dict[alignment_score] += 1
                        coordinate_cigar_combo = int(coordinate), cigar
                        if coordinate_cigar_combo not in coordinate_cigar_combos:
                            coordinate_cigar_combos[coordinate_cigar_combo] = 1
                        else:
                            coordinate_cigar_combos[coordinate_cigar_combo] += 1
                        if flag_as_string == '99':
                            flag_99 += 1
                        elif flag_as_string == '147':
                            flag_147 += 1
                        elif flag_as_string == '83':
                            flag_83 += 1
                        elif flag_as_string == '163':
                            flag_163 += 1
                        if perfect_gapped_alignment:
                            reads_with_perfect_gapped_alignments += 1
                            
                        if RNA_strand == '+':
                            if left_segment > max_overhang_upstream_junction:
                                max_overhang_upstream_junction = left_segment
                            if right_segment > max_overhang_downstream_junction:
                                max_overhang_downstream_junction = right_segment
                            if left_perfect_matches > max_overhang_perfect_matches_upstream_junction:
                                max_overhang_perfect_matches_upstream_junction = left_perfect_matches
                            if right_perfect_matches > max_overhang_perfect_matches_downstream_junction:
                                max_overhang_perfect_matches_downstream_junction = right_perfect_matches
                        else:
                            if right_segment > max_overhang_upstream_junction:
                                max_overhang_upstream_junction = right_segment
                            if left_segment > max_overhang_downstream_junction:
                                max_overhang_downstream_junction = left_segment
                            if right_perfect_matches > max_overhang_perfect_matches_upstream_junction:
                                max_overhang_perfect_matches_upstream_junction = right_perfect_matches
                            if left_perfect_matches > max_overhang_perfect_matches_downstream_junction:
                                max_overhang_perfect_matches_downstream_junction = left_perfect_matches
                if reads_with_gapped_alignments > 0:
                    supporting_junction_counts = []
                    for aligner in all_aligners:
                        supporting_junction_counts.append( str(aligner_to_supporting_junction_counts[ aligner ]) )
                    sorted_coordinate_cigar_combos = sorted(coordinate_cigar_combos.items(), key=operator.itemgetter(1), reverse = True)
                    sorted_mismatches_dict = sorted(mismatches_dict.items(), key=operator.itemgetter(1), reverse = True)
                    total_validated_introns += 1
                    output = [chromosome, start, end, RNA_strand, reads_without_mismatch_near_junction, sorted_coordinate_cigar_combos, sorted_mismatches_dict, reads_with_perfect_gapped_alignments, reads_with_gapped_alignments, unique_junction_counts, max_overhang_upstream_junction, max_overhang_perfect_matches_upstream_junction, max_overhang_downstream_junction, max_overhang_perfect_matches_downstream_junction, flag_99, flag_147, flag_83, flag_163] + supporting_junction_counts 
                    outfile.write('\t'.join([str(e) for e in output]) + '\n')
    outfile.close()      
    integrated_SAM_outfile.close()          
    return total_candidate_introns, total_validated_introns

alignment_aligner_lst = ALIGNMENTS_ALIGNERS.split('~')
aligner_to_processed_alignments = {}
all_aligners = []
for alignment_aligner in alignment_aligner_lst:
    alignment_filename, aligner = alignment_aligner.split(',')
    all_aligners.append(aligner)
    aligner_to_processed_alignments[aligner] = process_sam(alignment_filename)
print('filter_SAM_by_multiple_aligner_agreement')
read_IDs, read_IDs_with_incomplete_agreement = filter_SAM_by_multiple_aligner_agreement(aligner_to_processed_alignments, FILTERED_SAM_DIR, SAMPLE_PREFIX)
print('find_optimal_alignments')
integrated_alignments, SAM_header = find_optimal_alignments(all_aligners, aligner_to_processed_alignments, read_IDs, read_IDs_with_incomplete_agreement, INTEGRATED_ALIGNMENT_SUMMARY_OUTFILENAME)
print('process_splice_junctions_from_integrated_alignments')
total_candidate_introns, total_validated_introns = process_splice_junctions_from_integrated_alignments(all_aligners, integrated_alignments, SAM_header,  FILTERED_SAM_DIR, SAMPLE_PREFIX, JUNCTION_STATISTICS_OUTFILENAME)
