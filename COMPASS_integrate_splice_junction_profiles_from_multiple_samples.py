# -*- coding: utf-8 -*-
"""
Created on Fri Feb 16 22:05:47 2018

@author: kevinroy

Integrates data tables from each sample in the dataset into a single data table
Normalizes junction counts for each sample using the total uniquely mapped reads 
(as calculated by STAR) for each sample

For the outfile each row represents a different junction, and each column attributes of that junction, 
with across-sample invariate columns, as well as sample-specific columns
"""

import ast

PARENT_DIR= '/Volumes/SPxDrive/' # '/Users/kevinroy/Google_Drive/' # 
DIR=PARENT_DIR + 'splice_junction_analysis/'
SAMPLES='upf1_prp18 BY_upf1del_prp18del_2 BY_upf1del_prp18del_3 upf1 upf1_8.14.15 upf1_8.17.15'.split(' ') 

## normalization info
infilename = PARENT_DIR + 'splice_junction_analysis/upf1_prp18_first_pass_mapping_statistics.tsv'
with open(infilename, 'r') as infile:
    header = infile.readline()
    sample_name_to_uniquely_mapped_reads = {}
    for line in infile:
        info = line.strip().split('\t')
        sample_name, uniquely_mapped_reads, uniquely_mapped_gapped_reads = info[0], info[7], info[-1]
        sample_name_to_uniquely_mapped_reads[sample_name] = int(uniquely_mapped_reads)
    print(sample_name_to_uniquely_mapped_reads)
infile.close()

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
    
genome_file = PARENT_DIR + 'yeast_genome_references/saccharomyces_cerevisiae_sequence.fasta'

try:
    R64_genome
except:
    print('loading R64 genome')
    R64_genome = load_genome(genome_file)

def poly_U_PWMS(upstream_3SS_seq):
    '''
    takes the sequences upstream of the 3'SS
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
    print(U_score, 'U_score', upstream_3SS_seq)
    return U_score

poly_U_PWMS('T'*10)

junction_ID_to_common_info = {}
junction_ID_to_sample_to_info = {}
gene_name_to_junction_ID = {}

annotated_introns_to_sample_specific_info = {}

sample_unique_info = 'annotated_unique_junction_reads, intron_reads, total_spliced_and_unspliced_reads_associated_with_annotated_intron, annotated_SE, alternative_strand, SE, fraction_of_annotated_splicing, reads_without_mismatch_near_junction, coordinate_adjustment, top_coordinate_cigar_combo, mismatches_tallies, reads_with_perfect_gapped_alignments, reads_with_gapped_alignments, unique_junction_counts, max_overhang_upstream_junction, max_overhang_perfect_matches_upstream_junction, max_overhang_downstream_junction, max_overhang_perfect_matches_downstream_junction, BBMAP_and_STAR, BBMAP_only, STAR_only, BBMAP_intron_not_in_STAR, flag_99, flag_147, flag_83, flag_163, unique_junction_reads_rpm, annotated_unique_junction_reads_rpm, intron_reads_rpm, total_spliced_and_unspliced_reads_associated_with_annotated_intron_rpm'.split(', ')

NORMALIZATION_FACTOR = 100000000
MOST_FREQUENT_COORDINATE_CIGAR_COMBOS_TO_SHOW = 4  ## report only this many of the unique CIGARs for each junction, otherwise abundant junctions contain large numbers of unique CIGARs of low abundance

header_sample_unique_info = []
for sample_name in SAMPLES:
    total_mapped_reads = sample_name_to_uniquely_mapped_reads[sample_name]
    annotated_introns_to_sample_specific_info[sample_name] = {}
    header_sample_unique_info += [ (sample_name +':'+e) for e in sample_unique_info ]
    annotated_and_alternative_infilename = DIR + sample_name + '_combined_bbmap_STAR_splice_junctions_annotated_and_alternative_splicing_stats.txt'
    with open(annotated_and_alternative_infilename, 'r') as infile:
        header = infile.readline()
        for line in infile:
            info = line.strip().split('\t')
            
            ## check that the number of columns in each row matches the expected number of 80
            if len(info) != 80:
                print(len(info))
            else:

                annotated_junction_name, annotated_strand, systematic_gene_name, common_gene_name, RPG, annotated_intron_length, annotated_BP_3SS_dist, annotated_BP_5SS_dist, annotated_best_BP_sequence, annotated_coord_for_best_match_to_BP_consensus, annotated_poly_U_count, annotated_poly_Y_count, annotated_US_5SS_1, annotated_US_5SS_2, annotated_US_5SS_3, annotated_DS_3SS_1, annotated_DS_3SS_2, annotated_DS_3SS_3, annotated_upstream_5SS_seq, annotated_fiveSS_seq, annotated_downstream_5SS_seq, annotated_upstream_3SS_seq, annotated_threeSS_seq, annotated_downstream_3SS_seq, annotated_fiveSS_type, annotated_threeSS_type, junction_ID_name, chromosome, adjusted_start, adjusted_end, max_ambiguous_bases, alt_5SS, alt_3SS, intron_length, BP_3SS_dist, BP_5SS_dist, best_BP_sequence, coord_for_best_match_to_BP_consensus, poly_U_count, poly_Y_count, US_5SS_1, US_5SS_2, US_5SS_3, DS_3SS_1, DS_3SS_2, DS_3SS_3, upstream_5SS, fiveSS_seq, downstream_5SS, upstream_3SS, threeSS_seq, downstream_3SS, fiveSS_type, threeSS_type, annotated_unique_junction_reads, intron_reads, total_spliced_and_unspliced_reads_associated_with_annotated_intron, annotated_SE, coordinate_adjustment, alternative_strand, SE, fraction_of_annotated_splicing, reads_without_mismatch_near_junction, coordinate_cigar_combos, mismatches_tallies, reads_with_perfect_gapped_alignments, reads_with_gapped_alignments, unique_junction_counts, max_overhang_upstream_junction, max_overhang_perfect_matches_upstream_junction, max_overhang_downstream_junction, max_overhang_perfect_matches_downstream_junction, BBMAP_and_STAR, BBMAP_only, STAR_only, BBMAP_intron_not_in_STAR, flag_99, flag_147, flag_83, flag_163 = info
                ## read in coordinate cigar combos as a list of tuples                
                coordinate_cigar_combos = ( ast.literal_eval(coordinate_cigar_combos) )
                coordinate_cigar_combos = coordinate_cigar_combos[:MOST_FREQUENT_COORDINATE_CIGAR_COMBOS_TO_SHOW]
                
                if junction_ID_name not in junction_ID_to_sample_to_info:     
                    junction_ID_to_sample_to_info[junction_ID_name] = {}
                if systematic_gene_name not in gene_name_to_junction_ID:
                    gene_name_to_junction_ID[systematic_gene_name] = []
                if junction_ID_name not in gene_name_to_junction_ID[systematic_gene_name]:
                    gene_name_to_junction_ID[systematic_gene_name].append(junction_ID_name)
                
                ## normalize the read counts based on total mapped reads in each sample to reads per 100 million reads (rpm)
                # print( total_mapped_reads, unique_junction_counts, annotated_unique_junction_reads, intron_reads)
                unique_junction_reads_rpm, annotated_unique_junction_reads_rpm, intron_reads_rpm, total_spliced_and_unspliced_reads_associated_with_annotated_intron_rpm = float(unique_junction_counts) / float(total_mapped_reads) * NORMALIZATION_FACTOR, float(annotated_unique_junction_reads) / float(total_mapped_reads) * NORMALIZATION_FACTOR, float(intron_reads) / float(total_mapped_reads) * NORMALIZATION_FACTOR, float(total_spliced_and_unspliced_reads_associated_with_annotated_intron) / float(total_mapped_reads) * NORMALIZATION_FACTOR
                junction_ID_to_sample_to_info[junction_ID_name][sample_name] = annotated_unique_junction_reads, intron_reads, total_spliced_and_unspliced_reads_associated_with_annotated_intron, annotated_SE, alternative_strand, SE, fraction_of_annotated_splicing, reads_without_mismatch_near_junction, coordinate_adjustment, coordinate_cigar_combos, mismatches_tallies, reads_with_perfect_gapped_alignments, reads_with_gapped_alignments, unique_junction_counts, max_overhang_upstream_junction, max_overhang_perfect_matches_upstream_junction, max_overhang_downstream_junction, max_overhang_perfect_matches_downstream_junction, BBMAP_and_STAR, BBMAP_only, STAR_only, BBMAP_intron_not_in_STAR, flag_99, flag_147, flag_83, flag_163, unique_junction_reads_rpm, annotated_unique_junction_reads_rpm, intron_reads_rpm, total_spliced_and_unspliced_reads_associated_with_annotated_intron_rpm
                common_info = annotated_junction_name, annotated_strand, systematic_gene_name, common_gene_name, RPG, annotated_intron_length, annotated_BP_3SS_dist, annotated_BP_5SS_dist, annotated_best_BP_sequence, annotated_coord_for_best_match_to_BP_consensus, annotated_poly_U_count, annotated_poly_Y_count, annotated_US_5SS_1, annotated_US_5SS_2, annotated_US_5SS_3, annotated_DS_3SS_1, annotated_DS_3SS_2, annotated_DS_3SS_3, annotated_upstream_5SS_seq, annotated_fiveSS_seq, annotated_downstream_5SS_seq, annotated_upstream_3SS_seq, annotated_threeSS_seq, annotated_downstream_3SS_seq, annotated_fiveSS_type, annotated_threeSS_type, junction_ID_name, chromosome, adjusted_start, adjusted_end, max_ambiguous_bases, alt_5SS, alt_3SS, intron_length, BP_3SS_dist, BP_5SS_dist, best_BP_sequence, coord_for_best_match_to_BP_consensus, poly_U_count, poly_Y_count, US_5SS_1, US_5SS_2, US_5SS_3, DS_3SS_1, DS_3SS_2, DS_3SS_3, upstream_5SS, fiveSS_seq, downstream_5SS, upstream_3SS, threeSS_seq, downstream_3SS, fiveSS_type, threeSS_type
                if junction_ID_name not in junction_ID_to_common_info:
                    junction_ID_to_common_info[junction_ID_name] = common_info
                else:
                    if common_info != junction_ID_to_common_info[junction_ID_name]:
                        print('inconsistent junction info\n', common_info, '\n', junction_ID_to_common_info[junction_ID_name], '\n\n' )
                
                if annotated_junction_name not in annotated_introns_to_sample_specific_info[sample_name]:
                    ## ensure an annotated intron junction entry for each sample
                    ##, just in case this particular sample happens to have zero reads for the annotated junction (e.g. GCR1 for example)
                    annotated_introns_to_sample_specific_info[sample_name][annotated_junction_name] = [annotated_unique_junction_reads, intron_reads, total_spliced_and_unspliced_reads_associated_with_annotated_intron, annotated_SE, alternative_strand, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, annotated_unique_junction_reads_rpm, intron_reads_rpm, total_spliced_and_unspliced_reads_associated_with_annotated_intron_rpm]
    infile.close()
 
outfilename = DIR + 'combined_prp18_splice_junction_stats.txt'
outfile = open(outfilename , 'w')
header_common_info = 'annotated_chromosome, annotated_start, annotated_end, annotated_strand, systematic_gene_name, common_gene_name, RPG, annotated_intron_length, annotated_BP_3SS_dist, annotated_BP_5SS_dist, annotated_best_BP_sequence, annotated_coord_for_best_match_to_BP_consensus, annotated_U_score, alternative_U_score, annotated_poly_U_count, annotated_poly_Y_count, annotated_US_5SS_1, annotated_US_5SS_2, annotated_US_5SS_3, annotated_DS_3SS_1, annotated_DS_3SS_2, annotated_DS_3SS_3, annotated_upstream_5SS_seq, annotated_fiveSS_seq, annotated_downstream_5SS_seq, annotated_upstream_3SS_seq, annotated_threeSS_seq, annotated_downstream_3SS_seq, annotated_fiveSS_type, annotated_threeSS_type, junction_ID_name, chromosome, adjusted_start, adjusted_end, annotated_to_alternative_5SS_dist, annotated_to_alternative_3SS_dist, overhang_perfect_matches_upstream_junction, overhang_perfect_matches_downstream_junction, overhang_upstream_junction, overhang_downstream_junction, polyA_only_upstream_overhang, polyA_only_downstream_overhang, homopolymer_only_upstream_overhang, homopolymer_only_downstream_overhang, max_ambiguous_bases, alt_5SS, alt_3SS, intron_length, BP_3SS_dist, BP_5SS_dist, best_BP_sequence, coord_for_best_match_to_BP_consensus, poly_U_count, poly_Y_count, US_5SS_1, US_5SS_2, US_5SS_3, DS_3SS_1, DS_3SS_2, DS_3SS_3, upstream_5SS, fiveSS_seq, downstream_5SS, upstream_3SS, threeSS_seq, downstream_3SS, fiveSS_type, threeSS_type, total_reads_without_mismatch_near_junction, total_reads_with_perfect_gapped_alignments, total_gapped_reads, total_BBMAP_and_STAR, total_BBMAP_only, total_STAR_only, total_BBMAP_intron_not_in_STAR'.split(', ')
header = header_common_info + header_sample_unique_info
print(len(header))
header = '\t'.join(header) + '\n'
outfile.write(header)

for systematic_gene_name in gene_name_to_junction_ID:
    for junction_ID_name in gene_name_to_junction_ID[systematic_gene_name]:
        
        ## aggregate data across samples for these attributes:
        combined_sample_unique_info = []
        overhang_perfect_matches_upstream_junction = 0
        overhang_upstream_junction = 0
        overhang_downstream_junction = 0
        overhang_perfect_matches_downstream_junction = 0
        total_reads_without_mismatch_near_junction = 0
        total_reads_with_perfect_gapped_alignments = 0
        total_gapped_reads = 0
        total_BBMAP_and_STAR, total_BBMAP_only, total_STAR_only, total_BBMAP_intron_not_in_STAR = 0, 0, 0, 0    
        
        common_info = junction_ID_to_common_info[junction_ID_name]
        annotated_junction_name, annotated_strand, systematic_gene_name, common_gene_name, RPG, annotated_intron_length, annotated_BP_3SS_dist, annotated_BP_5SS_dist, annotated_best_BP_sequence, annotated_coord_for_best_match_to_BP_consensus, annotated_poly_U_count, annotated_poly_Y_count, annotated_US_5SS_1, annotated_US_5SS_2, annotated_US_5SS_3, annotated_DS_3SS_1, annotated_DS_3SS_2, annotated_DS_3SS_3, annotated_upstream_5SS_seq, annotated_fiveSS_seq, annotated_downstream_5SS_seq, annotated_upstream_3SS_seq, annotated_threeSS_seq, annotated_downstream_3SS_seq, annotated_fiveSS_type, annotated_threeSS_type, junction_ID_name, chromosome, adjusted_start, adjusted_end, max_ambiguous_bases, alt_5SS, alt_3SS, intron_length, BP_3SS_dist, BP_5SS_dist, best_BP_sequence, coord_for_best_match_to_BP_consensus, poly_U_count, poly_Y_count, US_5SS_1, US_5SS_2, US_5SS_3, DS_3SS_1, DS_3SS_2, DS_3SS_3, upstream_5SS, fiveSS_seq, downstream_5SS, upstream_3SS, threeSS_seq, downstream_3SS, fiveSS_type, threeSS_type = common_info
        annotated_chromosome, annotated_start, annotated_end, annotated_strand = annotated_junction_name.split('_')
        annotated_start, annotated_end = int(annotated_start), int(annotated_end)
        
        for sample_name in SAMPLES:
            if sample_name in junction_ID_to_sample_to_info[junction_ID_name]: ## check if the sample has any reads mapping to the junction, otherwise introduce zeros with the else statement below
                sample_unique_info = junction_ID_to_sample_to_info[junction_ID_name][sample_name]
                annotated_unique_junction_reads, intron_reads, total_spliced_and_unspliced_reads_associated_with_annotated_intron, annotated_SE, alternative_strand, SE, fraction_of_annotated_splicing, reads_without_mismatch_near_junction, coordinate_adjustment, coordinate_cigar_combos, mismatches_tallies, reads_with_perfect_gapped_alignments, reads_with_gapped_alignments, unique_junction_counts, max_overhang_upstream_junction, max_overhang_perfect_matches_upstream_junction, max_overhang_downstream_junction, max_overhang_perfect_matches_downstream_junction, BBMAP_and_STAR, BBMAP_only, STAR_only, BBMAP_intron_not_in_STAR, flag_99, flag_147, flag_83, flag_163, unique_junction_reads_rpm, annotated_unique_junction_reads_rpm, intron_reads_rpm, total_spliced_and_unspliced_reads_associated_with_annotated_intron_rpm = sample_unique_info
                overhang_perfect_matches_upstream_junction = max(overhang_perfect_matches_upstream_junction, int(max_overhang_perfect_matches_upstream_junction) )  ## max_overhang_perfect_matches_upstream_junction is the sample-specific value
                overhang_upstream_junction = max(overhang_upstream_junction, int(max_overhang_upstream_junction) )
                overhang_perfect_matches_downstream_junction = max(overhang_perfect_matches_downstream_junction, int(max_overhang_perfect_matches_downstream_junction) )
                overhang_downstream_junction = max(overhang_downstream_junction, int(max_overhang_downstream_junction) )   
                total_reads_without_mismatch_near_junction += int(reads_without_mismatch_near_junction)
                total_reads_with_perfect_gapped_alignments += int(reads_with_perfect_gapped_alignments)
                total_BBMAP_and_STAR += int(BBMAP_and_STAR)
                total_BBMAP_only += int(BBMAP_only)
                total_STAR_only += int(STAR_only)
                total_BBMAP_intron_not_in_STAR += int(BBMAP_intron_not_in_STAR)
                total_gapped_reads += int(reads_with_gapped_alignments)
                combined_sample_unique_info += [str(e) for e in sample_unique_info]
            else:
                if annotated_junction_name not in annotated_introns_to_sample_specific_info[sample_name]:
                    print(common_gene_name, sample_name)
                    combined_sample_unique_info += [str(0) for e in range( len(sample_unique_info) ) ]
                else:
                    combined_sample_unique_info += [str(e) for e in annotated_introns_to_sample_specific_info[sample_name][annotated_junction_name] ]

        annotated_U_score = poly_U_PWMS(annotated_upstream_3SS_seq)
        alternative_U_score = poly_U_PWMS(upstream_3SS)  # annotated_U_score, alternative_U_score, 
        ## polyA_only_upstream_overhang, polyA_only_downstream_overhang, homopolymer_only_upstream_overhang, homopolymer_only_downstream_overhang
        polyA_only_upstream_overhang = False
        polyA_only_downstream_overhang = False  
        homopolymer_only_upstream_overhang = False
        homopolymer_only_downstream_overhang = False 
        adjusted_start = int(adjusted_start)
        adjusted_end = int(adjusted_end)
        if annotated_strand == '+':
            annotated_to_alternative_5SS_dist = adjusted_start - annotated_start
            annotated_to_alternative_3SS_dist = adjusted_end - annotated_end
        else:
            annotated_to_alternative_3SS_dist = annotated_start - adjusted_start
            annotated_to_alternative_5SS_dist = annotated_end - adjusted_end
        upstream_overhang_seq = R64_genome[chromosome][adjusted_start - overhang_perfect_matches_upstream_junction - 1: adjusted_start - 1]
        downstream_overhang_seq = R64_genome[chromosome][adjusted_end: adjusted_end + overhang_perfect_matches_downstream_junction]
        if annotated_strand == '-':    
            upstream_overhang_seq =  rev_comp(downstream_overhang_seq)
            downstream_overhang_seq = rev_comp(upstream_overhang_seq)
        if upstream_overhang_seq == len(upstream_overhang_seq) * 'A':
            polyA_only_upstream_overhang = True
        if downstream_overhang_seq == len(downstream_overhang_seq) * 'A':
            polyA_only_downstream_overhang = True
        if upstream_overhang_seq == len(upstream_overhang_seq) * 'A' or upstream_overhang_seq == len(upstream_overhang_seq) * 'T' or upstream_overhang_seq == len(upstream_overhang_seq) * 'G' or upstream_overhang_seq == len(upstream_overhang_seq) * 'C':
            homopolymer_only_upstream_overhang = True
        if downstream_overhang_seq == len(downstream_overhang_seq) * 'A' or downstream_overhang_seq == len(downstream_overhang_seq) * 'T' or downstream_overhang_seq == len(downstream_overhang_seq) * 'G' or downstream_overhang_seq == len(downstream_overhang_seq) * 'C':
            homopolymer_only_downstream_overhang = True  
        ## Data table stratified into common info first, followed by sample specific info
        ## Common info contains the annotated intron info first (overlapping annotated junction as determined by closest 5'SS), followed by the alternative intron info
        updated_annotated_common_info = annotated_chromosome, annotated_start, annotated_end, annotated_strand, systematic_gene_name, common_gene_name, RPG, annotated_intron_length, annotated_BP_3SS_dist, annotated_BP_5SS_dist, annotated_best_BP_sequence, annotated_coord_for_best_match_to_BP_consensus, annotated_U_score, alternative_U_score, annotated_poly_U_count, annotated_poly_Y_count, annotated_US_5SS_1, annotated_US_5SS_2, annotated_US_5SS_3, annotated_DS_3SS_1, annotated_DS_3SS_2, annotated_DS_3SS_3, annotated_upstream_5SS_seq, annotated_fiveSS_seq, annotated_downstream_5SS_seq, annotated_upstream_3SS_seq, annotated_threeSS_seq, annotated_downstream_3SS_seq, annotated_fiveSS_type, annotated_threeSS_type
        ## note that each row represents a different junction overlapping, or identical to, to an annotated intron junction
        updated_junction_common_info = junction_ID_name, chromosome, adjusted_start, adjusted_end, annotated_to_alternative_5SS_dist, annotated_to_alternative_3SS_dist, overhang_perfect_matches_upstream_junction, overhang_perfect_matches_downstream_junction, overhang_upstream_junction, overhang_downstream_junction, polyA_only_upstream_overhang, polyA_only_downstream_overhang, homopolymer_only_upstream_overhang, homopolymer_only_downstream_overhang, max_ambiguous_bases, alt_5SS, alt_3SS, intron_length, BP_3SS_dist, BP_5SS_dist, best_BP_sequence, coord_for_best_match_to_BP_consensus, poly_U_count, poly_Y_count, US_5SS_1, US_5SS_2, US_5SS_3, DS_3SS_1, DS_3SS_2, DS_3SS_3, upstream_5SS, fiveSS_seq, downstream_5SS, upstream_3SS, threeSS_seq, downstream_3SS, fiveSS_type, threeSS_type, total_reads_without_mismatch_near_junction, total_reads_with_perfect_gapped_alignments, total_gapped_reads, total_BBMAP_and_STAR, total_BBMAP_only, total_STAR_only, total_BBMAP_intron_not_in_STAR
        output = [str(e) for e in updated_annotated_common_info] + [str(e) for e in updated_junction_common_info] + combined_sample_unique_info
        output = '\t'.join(output) + '\n'
        outfile.write(output)
outfile.close()
                
            
