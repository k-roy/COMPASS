# -*- coding: utf-8 -*-
"""
Created on Fri Feb 16 22:05:47 2018

@author: kevinroy

takes all junctions aggregated by the script integrate_splice_junction_profiles_from_multiple_samples.py
writes the sequences encompassing the BP-3'SS, and downstream bases of 5, 10, 15, and 20 to fasta file

After running this script, need to utilize the Vienna RNA package and generate RNA folds for each BP-3'SS sequence:

DIR=/Volumes/SPxDrive/prp18_RNA_folding_50nt
mkdir $DIR
cd $DIR
rnafold -p --constraint --batch --infile=annotated_BP_3SS.fasta --outfile=annotated_BP_3SS_folds.txt

rnafold -p --constraint --batch --infile=all_BP_3SS.fasta --outfile=all_BP_3SS_folds.txt

rnafold -p --constraint --batch --infile=alt_BP_3SS.fasta --outfile=alt_BP_3SS_folds.txt

for filename in *dp.ps;
do prefix=${filename%_dp.ps};
echo $prefix;
/usr/local/share/ViennaRNA/bin/mountain.pl $filename > $prefix\_mountain.txt;
done
"""

PARENT_DIR= '/Volumes/SPxDrive/' # '/Users/kevinroy/Google_Drive/' # 
DIR=PARENT_DIR + 'splice_junction_analysis/'
genome_file = PARENT_DIR + 'yeast_genome_references/saccharomyces_cerevisiae_sequence.fasta'
 
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
    input: = DNA string
    output: reverse complement string
    '''
    rev_comp_DNA = ''
    comp_bases = {'A':'T', 'T':'A', 'G':'C', 'C':'G'}
    for base in DNA[::-1]:
        rev_comp_DNA += comp_bases[base]
    return rev_comp_DNA

try:
    R64_genome
except:
    print('loading R64 genome')
    R64_genome = load_genome(genome_file)

infilename =  DIR + 'combined_averaged_prp18_splice_junction_stats.txt'

FOLD_DIR = PARENT_DIR + 'prp18_RNA_folding_50nt/'
outfilename = FOLD_DIR + 'annotated_BP_3SS.fasta'
annotated_BP_3SS_outfile = open(outfilename, 'w')
outfilename = FOLD_DIR + 'alt_BP_3SS.fasta'
alt_BP_3SS_outfile = open(outfilename, 'w')

outfilename = FOLD_DIR + 'all_BP_3SS.fasta'
all_BP_3SS_outfile = open(outfilename, 'w')

with open(infilename, 'r') as infile:
    header = infile.readline()
    for line in infile:
        info = line.strip().split()
        annotated_chromosome, annotated_start, annotated_end, annotated_strand, systematic_gene_name, common_gene_name, RPG, annotated_intron_length, annotated_BP_3SS_dist, annotated_BP_5SS_dist, annotated_best_BP_sequence, annotated_coord_for_best_match_to_BP_consensus, annotated_U_score, alternative_U_score, annotated_poly_U_count, annotated_poly_Y_count, annotated_US_5SS_1, annotated_US_5SS_2, annotated_US_5SS_3, annotated_DS_3SS_1, annotated_DS_3SS_2, annotated_DS_3SS_3, annotated_upstream_5SS_seq, annotated_fiveSS_seq, annotated_downstream_5SS_seq, annotated_upstream_3SS_seq, annotated_threeSS_seq, annotated_downstream_3SS_seq, annotated_fiveSS_type, annotated_threeSS_type, junction_ID_name, chromosome, adjusted_start, adjusted_end, annotated_to_alternative_5SS_dist, annotated_to_alternative_3SS_dist, overhang_perfect_matches_upstream_junction, overhang_perfect_matches_downstream_junction, overhang_upstream_junction, overhang_downstream_junction, polyA_only_upstream_overhang, polyA_only_downstream_overhang, homopolymer_only_upstream_overhang, homopolymer_only_downstream_overhang, max_ambiguous_bases, alt_5SS, alt_3SS, intron_length, BP_3SS_dist, BP_5SS_dist, best_BP_sequence, coord_for_best_match_to_BP_consensus, poly_U_count, poly_Y_count, US_5SS_1, US_5SS_2, US_5SS_3, DS_3SS_1, DS_3SS_2, DS_3SS_3, upstream_5SS, fiveSS_seq, downstream_5SS, upstream_3SS, threeSS_seq, downstream_3SS, fiveSS_type, threeSS_type, total_reads_without_mismatch_near_junction, total_reads_with_perfect_gapped_alignments, total_gapped_reads, total_BBMAP_and_STAR, total_BBMAP_only, total_STAR_only, total_BBMAP_intron_not_in_STAR = info[:75]
       #  chromosome, adjusted_start, adjusted_end, annotated_strand, annotated_coord_for_best_match_to_BP_consensus
        annotated_to_alternative_5SS_dist, annotated_to_alternative_3SS_dist, annotated_coord_for_best_match_to_BP_consensus, coord_for_best_match_to_BP_consensus, adjusted_start, adjusted_end = [int(e) for e in [annotated_to_alternative_5SS_dist, annotated_to_alternative_3SS_dist, annotated_coord_for_best_match_to_BP_consensus, coord_for_best_match_to_BP_consensus, adjusted_start, adjusted_end]]
        for DS_bp_to_add in 50, 55, 60, 65:
            if annotated_strand == '+':
                BP_3SS_length = adjusted_end - (coord_for_best_match_to_BP_consensus - 6)
                seq_to_fold =   R64_genome[chromosome][coord_for_best_match_to_BP_consensus -6:adjusted_end + DS_bp_to_add]  
            else:
                BP_3SS_length = (coord_for_best_match_to_BP_consensus + 5) - adjusted_start + 1
                seq_to_fold =   rev_comp( R64_genome[chromosome][adjusted_start - 1 - DS_bp_to_add:coord_for_best_match_to_BP_consensus + 5 ] )
                if seq_to_fold[3:7] != 'TAAC':
                    print(common_gene_name, junction_ID_name, annotated_best_BP_sequence, seq_to_fold, '\n')
            ## force bases just downstream of BP to be unpaired, per Meyer et al. Mol Cell, 2012 (PMID: 21925391; DOI: 10.1016/j.molcel.2011.07.030)
            unpaired_bases = min(15, BP_3SS_length)
            length_unrestricted_bases =  BP_3SS_length + DS_bp_to_add - unpaired_bases
            constraint_string = 'x'*unpaired_bases + '.'*length_unrestricted_bases
            seq_name = '>' + '_'.join([junction_ID_name, common_gene_name, 'annotated_BP:' + str(annotated_coord_for_best_match_to_BP_consensus), 'BP:' + str(coord_for_best_match_to_BP_consensus), 'DS_bp:' + str(DS_bp_to_add)])   

            if annotated_to_alternative_5SS_dist == 0 and annotated_to_alternative_3SS_dist == 0:
                annotated_BP_3SS_outfile.write(seq_name + '\n' + seq_to_fold + '\n' +  constraint_string + '\n' )
            else:
                alt_BP_3SS_outfile.write(seq_name + '\n' + seq_to_fold + '\n' +  constraint_string + '\n' )
            all_BP_3SS_outfile.write(seq_name + '\n' + seq_to_fold + '\n' +  constraint_string + '\n' )


all_BP_3SS_outfile.close()
annotated_BP_3SS_outfile.close()
alt_BP_3SS_outfile.close()


