import operator
import pysam
import pandas as pd

def get_ambiguous_junctions(chrom, start, stop, genome_fasta):
    '''
    takes intron coordinates
    checks for ambiguous junction based on nt upstream and 
    downstream of exon-intron and intron-exon junctions
    returns list of ambiguous junctions
    '''
    junctions = []
    idx = 1
    while start-idx > 0 and genome_fasta.fetch(chrom, start-idx, start-idx+1).upper() \
        == genome_fasta.fetch(chrom, stop-idx+1, stop-idx+2).upper():
        ambiguous_intron = (chrom, start-idx, stop-idx)
        junctions.append(ambiguous_intron)
        idx += 1
    idx = 0
    while stop-idx+1 < genome_fasta.get_reference_length(chrom) and \
        genome_fasta.fetch(chrom, start-idx, start-idx+1).upper() == \
            genome_fasta.fetch(chrom, stop-idx+1, stop-idx+2).upper():
        ambiguous_intron = (chrom, start-idx+1, stop-idx+1)
        junctions.append(ambiguous_intron)
        idx -= 1
    return junctions

def get_ambiguous_junctions_in_annotated_introns(introns_file, genome_fasta):
    '''
    takes a file of introns in gff format with 1-based coords
    returns annotated_intron_df, ambiguous_junction_to_annotated_junction, ambiguous_annotated_junctions, annotated_introns, junction_to_intron_type
    '''
    annotated_intron_df = pd.read_csv(introns_file, sep = '\t', names=['chrom','start', 'stop', 'strand'])
    # adjust 1-based introns to 0-based python coords
    annotated_intron_df['adjusted_start'] = annotated_intron_df['start'] + 1 # - 1
    # The HISAT2 splice sites file already has the start coord 0-based
    annotated_intron_df['adjusted_stop'] = annotated_intron_df['stop'] - 1
    ambiguous_junction_to_annotated_junction = {}
    ambiguous_annotated_junctions = set([])
    annotated_introns = set([])
    junction_to_intron_type = {}
    for index, row in annotated_intron_df.iterrows():
        chrom = row['chrom']
        start = row['adjusted_start']
        stop = row['adjusted_stop']
        annotated_intron = (chrom, start, stop)
        ambiguous_junctions = get_ambiguous_junctions(chrom, start, stop, genome_fasta)
        junction_to_intron_type[annotated_intron] = 'intron'
        annotated_introns.add(annotated_intron)
        if ambiguous_junctions != []:
            ambiguous_annotated_junctions.add(annotated_intron)
            for ambiguous_junction in ambiguous_junctions:
                ambiguous_junction_to_annotated_junction[ambiguous_junction] \
                    = annotated_intron
                junction_to_intron_type[ambiguous_junction] = 'intron'
    print('number of annotated introns:', len(annotated_intron_df.index))
    print('number of annotated introns with identical nt upstream or downstream of junctions:', len(ambiguous_annotated_junctions))
    print('total number of intron coords mapping to annotated junctions:', len(ambiguous_junction_to_annotated_junction))
    return annotated_intron_df, ambiguous_junction_to_annotated_junction, ambiguous_annotated_junctions, annotated_introns, junction_to_intron_type

def parse_cigar(cigar, minimum_intron_length):
    '''
    input: CIGAR string from SAM alignment
    output: list of 2 item tuples: mapped segment length, CIGAR operation
    '''
    total_cigar_operation_bp = {'N':0, 'X':0, 'I':0, 'D':0, 'S':0, 'H':0, '=':0, 'M':0 }
    segment_length = ''
    cigartuples = []
    for char in cigar:
        if char in '0123456789':
            segment_length += char
        else:
            operation = char
            mapped_segment_length = int(segment_length)
            # convert D to N for BBMap, as it does not use N for introns
            if mapped_segment_length >= minimum_intron_length and operation == 'D':
                operation = 'N'
            elif mapped_segment_length < minimum_intron_length and operation == 'N':
                operation = 'D'
            segment_length = ''
            cigartuples.append((mapped_segment_length, operation))
            total_cigar_operation_bp[operation] += mapped_segment_length
    return cigartuples, total_cigar_operation_bp

def hamming(a, b):
    return len([i for i in filter(lambda x: x[0] != x[1], zip(a, b))])
    
def get_five_SS_score(query, consensus_5SS, penalties_5SS):
    return sum(penalties_5SS[i] for i in range(len(query)) if query[i] not in consensus_5SS[i])

def get_three_SS_score(query, consensus_3SS, penalties_3SS):
    return sum(penalties_3SS[i] for i in range(len(query)) if query[i] not in consensus_3SS[i])

def rc(seq):
    return seq.translate(str.maketrans("ACTG", "TGAC"))[::-1]

def count_frequency(lst, num_items_to_report):
    counts = {}
    for i in lst: counts[i] = counts.get(i, 0) + 1
    sort_counts= sorted(counts.items(), key=operator.itemgetter(1))
    return sort_counts[:num_items_to_report]

def adjust_ambiguous_junctions(chrom, start, stop, RNA_strand, genome_fasta, \
    annotated_intron_df, annotated_introns, consensus_5SS, penalties_5SS, \
        consensus_3SS, penalties_3SS):
    '''
    This function takes intron coordinates and
    checks for ambiguous junction based on nt upstream and downstream of 
    exon-intron and intron-exon junctions.
    if ambiguous, adjusts based on closest match to 5SS and 3SS consensus. 
    The 3SS motif adds a smaller scoring component as it is less conserved in 
    S. cerevisiae.  This can be changed for other organisms.
    # max penalty score for 5SS = 10, max penalty score for 3SS = 6
    returns adjusted coordinates and other info on the intron: 
    five_SS, three_SS, RNA_strand, num_amb_junctions, annotated_junction, ann_5SS, ann_3SS, canonical_5SS, canonical_3SS, intron_size, intron_coords_adjusted
    '''
    junctions = [(chrom, start, stop)]
    idx = 1
    # don't consider junctions with 10 identical nt upstream or downstream 
    # as these are almost certainly artifactual gapped alignments
    while idx < 10 and start-idx > 0 and genome_fasta.fetch(chrom, start-idx, start-idx+1).upper() \
        == genome_fasta.fetch(chrom, stop-idx+1, stop-idx+2).upper():
        #print(chrom, start, stop, idx, genome_fasta.fetch(chrom, start-idx, start-idx+1), 
        # genome_fasta.fetch(chrom, stop-idx+1, stop-idx+2))
        ambiguous_intron = (chrom, start-idx, stop-idx)
        junctions.append(ambiguous_intron)
        idx += 1
    idx = 0
    while idx > -10 and stop-idx+1 < genome_fasta.get_reference_length(chrom) and \
        genome_fasta.fetch(chrom, start-idx, start-idx+1).upper() == \
            genome_fasta.fetch(chrom, stop-idx+1, stop-idx+2).upper():
        #print(chrom, start, stop, idx, genome_fasta.fetch(chrom, start-idx, start-idx+1), 
        # genome_fasta.fetch(chrom, stop-idx+1, stop-idx+2))
        ambiguous_intron = (chrom, start-idx+1, stop-idx+1)
        junctions.append(ambiguous_intron)
        idx -= 1
    potential_5SS = []
    potential_3SS = []
    SS_scores = []
    junctions = sorted(junctions) 
    # This is critical to ensure that different ambiguous junction inputs 
    # yield identical adjusted junctions, 
    # e.g. in cases where there are two potential amb junctions with identical scores.
    for junction in junctions:
        chrom, amb_start, amb_stop = junction
        if RNA_strand == '+':
            fiveSS = genome_fasta.fetch(chrom, amb_start, amb_start+6).upper()
            threeSS = genome_fasta.fetch(chrom, amb_stop-2, amb_stop+1).upper()
        if RNA_strand == '-':
            fiveSS = rc(genome_fasta.fetch(chrom, amb_stop-5, amb_stop+1)).upper()
            threeSS = rc(genome_fasta.fetch(chrom, amb_start, amb_start+3)).upper()
        potential_5SS.append(fiveSS)
        potential_3SS.append(threeSS)
        five_SS_score = get_five_SS_score(fiveSS, consensus_5SS, penalties_5SS)
        three_SS_score = get_three_SS_score(threeSS, consensus_3SS, penalties_3SS)
        SS_scores.append(five_SS_score + three_SS_score)
    most_likely_intron = junctions[SS_scores.index(min(SS_scores))] 
    chrom, adj_start, adj_stop = most_likely_intron
    intron_coords_adjusted = (adj_start == start)
    if RNA_strand == '+':
        ann_5SS = ((annotated_intron_df['strand'] == RNA_strand) & (annotated_intron_df['adjusted_start'] == adj_start)).any()
        ann_3SS = ((annotated_intron_df['strand'] == RNA_strand) & (annotated_intron_df['adjusted_stop'] == adj_stop)).any()
    else:
        ann_3SS = ((annotated_intron_df['strand'] == RNA_strand) & (annotated_intron_df['adjusted_start'] == adj_start)).any()
        ann_5SS = ((annotated_intron_df['strand'] == RNA_strand) & (annotated_intron_df['adjusted_stop'] == adj_stop)).any()
    five_SS = potential_5SS[SS_scores.index(min(SS_scores))]
    three_SS = potential_3SS[SS_scores.index(min(SS_scores))]
    canonical_5SS = (five_SS[:2] == 'GT')
    canonical_3SS = (three_SS[1:] == 'AG')
    intron_size = abs(stop - start)
    annotated_junction = (most_likely_intron in annotated_introns)
    num_amb_junctions = len(potential_5SS)
    # These key variables returned by this script are: chrom, adj_start, adj_stop, annotated_junction. 
    # Because analyzing the motifs is required for the adjustment process, this info is returned back by this function and saved in a dictionary.
    # The dictionary
    # The rest of the returned values in the list could be generated later in the pipeline, but it is convenient to already process this here. 
    return [chrom, adj_start, adj_stop, five_SS, three_SS, RNA_strand, num_amb_junctions, annotated_junction, ann_5SS, ann_3SS, canonical_5SS, canonical_3SS, intron_size, intron_coords_adjusted]

def process_alignment(read, quality_scores_of_current_read_num, junctions_to_ambiguous_junctions, genome_fasta, \
    annotated_intron_df, annotated_introns, consensus_5SS, penalties_5SS, consensus_3SS, \
        penalties_3SS, MISMATCH_DIST_FROM_JUNCTION_DISALLOWED, minimum_intron_length):
    '''
    takes a read (a pysam AlignedSegment object), junctions_to_ambiguous_junctions (a growing dictionary containing output of
    adjust_ambiguous_junctions for each encountered junction), and two global variables

    returns dictionary of analysis on the read alignment, and an updated dictionary junctions_to_ambiguous_junctions
    '''
    common_intron_info = 'chrom, adj_start, adj_stop, five_SS, three_SS, RNA_strand, num_amb_junctions, annotated_junction, ann_5SS, ann_3SS, canonical_5SS, canonical_3SS, intron_size'.split(', ')

    sample_specific_intron_info = 'intron_coords_adjusted, mismatch_near_junction, US_perfect_matches, DS_perfect_matches, five_SS_Q_score, three_SS_Q_score'.split(', ')

    adjusted_intron_keys = common_intron_info + sample_specific_intron_info

    keys = 'flag, chrom, coord, cigar, NH, adjusted_introns, perfect_gapped_alignment, splice_sites, alignment_score, RNA_strand'.split(', ')          
    
    if type(read) != pysam.libcalignedsegment.AlignedSegment or read.cigartuples == None:
        # values = read.flag, read.reference_name, read.reference_start, read.cigarstring.replace('N', 'D'), read.get_tag('NH'), adjusted_introns, alignment_score
        values = [None] * 5 + [[], False, [], 1000]  # alignment score is set to 1000 for unmapped reads
        return dict(zip(keys, values)), junctions_to_ambiguous_junctions
    
    quality_scores = read.query_qualities
    if quality_scores == None: # MAGIC_BLAST does not report quality scores for the reads
        quality_scores = quality_scores_of_current_read_num
        if quality_scores == []:
            quality_scores = [0]*150
    cigartuples, total_cigar_operation_bp = parse_cigar(read.cigarstring, minimum_intron_length)
    coordinate = read.reference_start
    chrom = read.reference_name
    RNA_strand = '+' if ((read.is_read1 and read.is_reverse) or (read.is_read2 and not read.is_reverse)) else '-'
    # Parsing the cigar serves three main purposes: 
    # (1) assess whether alignment has a gap within the allowable range for an intron
    # (2) calculate the putative intron boundaries (i.e. 5' and 3' splice sites)
    # (3) gather metrics on the alignment, including the total mismatches, 
    # lengths of left and right mapped segments around the gap,
    # and the length of each segment matching perfectly on the left and right of the gap.
    mapped_segment_length = 0
    total_gapped_bp = 0
    alignment_score = 0
    splice_sites = []  # a single read may have more than one splice site
    quality_scores_at_splice_sites = []
    perfect_matches_flanking_splice_sites = []
    mapped_segment_lengths = []
    mismatch_near_junction = []
    five_SS_Q_score = 0
    three_SS_Q_score = 0
    for cigar_idx in range(len(cigartuples)):  
        bp, operation = cigartuples[cigar_idx]
        if operation in ('=', 'X', 'D', 'M'):
            mapped_segment_length += bp
            coordinate += bp
        elif operation == 'N':
            total_gapped_bp += bp
            mapped_segment_lengths.append(mapped_segment_length)
            mapped_segment_length = 0
            # the cigar operation flanking N should always be =, 
            # but need to check this just in case
            if cigartuples[cigar_idx-1][1] != '=':
                left_perfect_matches = 0
            else:
                left_perfect_matches = cigartuples[cigar_idx-1][0]
            if cigartuples[cigar_idx+1][1] != '=' :
                right_perfect_matches = 0
            else:
                right_perfect_matches = cigartuples[cigar_idx+1][0]   
            if bp >= minimum_intron_length:
                splice_sites.append((coordinate, coordinate + bp - 1))
                if RNA_strand == '+':
                    US_perfect_matches = left_perfect_matches
                    DS_perfect_matches = right_perfect_matches
                    five_SS_Q_score = quality_scores[mapped_segment_length - 1]
                    three_SS_Q_score = quality_scores[mapped_segment_length]
                else:
                    DS_perfect_matches = left_perfect_matches
                    US_perfect_matches = right_perfect_matches
                    three_SS_Q_score = quality_scores[mapped_segment_length - 1]
                    five_SS_Q_score = quality_scores[mapped_segment_length]
                quality_scores_at_splice_sites.append((five_SS_Q_score, three_SS_Q_score))
                coordinate = coordinate + bp
                perfect_matches_flanking_splice_sites.append((US_perfect_matches, DS_perfect_matches))
                mismatch_near_junction.append((cigar_idx >= 2 and left_perfect_matches < MISMATCH_DIST_FROM_JUNCTION_DISALLOWED) or (cigar_idx+2 < len(cigartuples) and right_perfect_matches < MISMATCH_DIST_FROM_JUNCTION_DISALLOWED))
            else:
                alignment_score += bp
        if operation in ('X', 'I', 'S', 'H', 'D', 'M'):
            alignment_score += bp
    mapped_segment_lengths.append(mapped_segment_length)
    
    adjusted_introns = []
    for idx in range(len(splice_sites)):
        start, stop = splice_sites[idx]
        US_perfect_matches, DS_perfect_matches = perfect_matches_flanking_splice_sites[idx]
        five_SS_Q_score, three_SS_Q_score = quality_scores_at_splice_sites[idx]
        intron = chrom, start, stop
        if intron not in junctions_to_ambiguous_junctions:
            junctions_to_ambiguous_junctions[intron] = adjust_ambiguous_junctions(chrom, start, stop, RNA_strand, genome_fasta, annotated_intron_df, annotated_introns, consensus_5SS, penalties_5SS, consensus_3SS, penalties_3SS)
            # adjust_ambiguous_junctions returns: chrom, amb_start, amb_stop, five_SS, three_SS, RNA_strand, num_amb_junctions, annotated_junction, ann_5SS, ann_3SS, canonical_5SS, canonical_3SS, intron_size, intron_coords_adjusted, 
            # additional info added: mismatch_near_junction, US_perfect_matches, DS_perfect_matches
            # values in adjusted_introns list: 
            # chrom, adj_start, adj_stop, five_SS, three_SS, RNA_strand, num_amb_junctions, annotated_junction, ann_5SS, ann_3SS, canonical_5SS, canonical_3SS, intron_size, intron_coords_adjusted, mismatch_near_junction, US_perfect_matches, DS_perfect_matches
        adjusted_intron_values = junctions_to_ambiguous_junctions[intron] + [mismatch_near_junction[idx], US_perfect_matches, DS_perfect_matches, five_SS_Q_score, three_SS_Q_score]
        adjusted_introns.append(dict(zip(adjusted_intron_keys, adjusted_intron_values)))
    perfect_gapped_alignment = (splice_sites != []) & (alignment_score == 0)
    values = read.flag, read.reference_name, read.reference_start, read.cigarstring.replace('D', 'N').replace('M', 'X'), read.get_tag('NH'), adjusted_introns, perfect_gapped_alignment, splice_sites, alignment_score, RNA_strand
    return dict(zip(keys, values)), junctions_to_ambiguous_junctions
