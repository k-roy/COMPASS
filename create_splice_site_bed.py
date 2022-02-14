import pandas as pd
from sys import argv
COMPASS_JUNCTIONS_DIR = argv[1]
sample_name = argv[2]

splice_junction_counts_infilename = COMPASS_JUNCTIONS_DIR + sample_name + '_COMPASS_splice_junctions_with_sequence_info.tsv'

junction_df = pd.read_csv(splice_junction_counts_infilename, sep = '\t')
junction_df.columns

junction_df['adj_start_plus_one'] = junction_df['adj_start'] + 1
junction_df['adj_stop_plus_one'] = junction_df['adj_stop'] + 1

plus_strand_fiveSS_splice_site_bed_outfilename = COMPASS_JUNCTIONS_DIR + sample_name + '_plus_strand_fiveSS_splice_site.bed'
plus_strand_threeSS_splice_site_bed_outfilename = COMPASS_JUNCTIONS_DIR + sample_name + '_plus_strand_threeSS_splice_site.bed'
minus_strand_fiveSS_splice_site_bed_outfilename = COMPASS_JUNCTIONS_DIR + sample_name + '_minus_strand_fiveSS_splice_site.bed'
minus_strand_threeSS_splice_site_bed_outfilename = COMPASS_JUNCTIONS_DIR + sample_name + '_minus_strand_threeSS_splice_site.bed'

junction_df[junction_df['RNA_strand'] == '+'][['chrom', 'adj_start', 'adj_start_plus_one']].to_csv(plus_strand_fiveSS_splice_site_bed_outfilename, sep='\t', index=False, header=False)
junction_df[junction_df['RNA_strand'] == '+'][['chrom', 'adj_stop', 'adj_stop_plus_one']].to_csv(plus_strand_threeSS_splice_site_bed_outfilename, sep='\t', index=False, header=False)
junction_df[junction_df['RNA_strand'] == '-'][['chrom', 'adj_stop', 'adj_stop_plus_one']].to_csv(minus_strand_fiveSS_splice_site_bed_outfilename, sep='\t', index=False, header=False)
junction_df[junction_df['RNA_strand'] == '-'][['chrom', 'adj_start', 'adj_start_plus_one']].to_csv(minus_strand_threeSS_splice_site_bed_outfilename, sep='\t', index=False, header=False)