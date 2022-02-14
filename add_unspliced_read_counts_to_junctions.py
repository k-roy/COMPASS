import pandas as pd
from sys import argv
COMPASS_JUNCTIONS_DIR = argv[1]
sample_name = argv[2]

splice_junction_counts_infilename = COMPASS_JUNCTIONS_DIR + sample_name + '_COMPASS_splice_junctions_with_sequence_info.tsv'
splice_junction_counts_outfilename = COMPASS_JUNCTIONS_DIR + sample_name + '_COMPASS_splice_junctions_with_sequence_info_and_unspliced_signal.tsv'

junction_df = pd.read_csv(splice_junction_counts_infilename, sep = '\t')
junction_df.columns

junction_df['adj_start_plus_one'] = junction_df['adj_start'] + 1
junction_df['adj_stop_plus_one'] = junction_df['adj_stop'] + 1

plus_strand_fiveSS_splice_site_coverage_bed_outfilename = COMPASS_JUNCTIONS_DIR + sample_name + '_plus_strand_fiveSS_splice_site_coverage.bed'
plus_strand_threeSS_splice_site_coverage_bed_outfilename = COMPASS_JUNCTIONS_DIR + sample_name + '_plus_strand_threeSS_splice_site_coverage.bed'
minus_strand_fiveSS_splice_site_coverage_bed_outfilename = COMPASS_JUNCTIONS_DIR + sample_name + '_minus_strand_fiveSS_splice_site_coverage.bed'
minus_strand_threeSS_splice_site_coverage_bed_outfilename = COMPASS_JUNCTIONS_DIR + sample_name + '_minus_strand_threeSS_splice_site_coverage.bed'

plus_strand_fiveSS_splice_site_coverage_bed = pd.read_csv(plus_strand_fiveSS_splice_site_coverage_bed_outfilename, sep = '\t', names = ['chrom', 'adj_start_plus_one', 'fiveSS_coverage'])
plus_strand_fiveSS_splice_site_coverage_bed['RNA_strand'] = '+'
plus_strand_fiveSS_splice_site_coverage_bed

plus_strand_threeSS_splice_site_coverage_bed = pd.read_csv(plus_strand_threeSS_splice_site_coverage_bed_outfilename, sep = '\t', names = ['chrom', 'adj_stop_plus_one', 'threeSS_coverage'])
plus_strand_threeSS_splice_site_coverage_bed['RNA_strand'] = '+'
plus_strand_threeSS_splice_site_coverage_bed

minus_strand_fiveSS_splice_site_coverage_bed = pd.read_csv(minus_strand_fiveSS_splice_site_coverage_bed_outfilename, sep = '\t', names = ['chrom', 'adj_stop_plus_one', 'fiveSS_coverage'])
minus_strand_fiveSS_splice_site_coverage_bed['RNA_strand'] = '-'
minus_strand_fiveSS_splice_site_coverage_bed

minus_strand_threeSS_splice_site_coverage_bed = pd.read_csv(minus_strand_threeSS_splice_site_coverage_bed_outfilename, sep = '\t', names = ['chrom', 'adj_start_plus_one', 'threeSS_coverage'])
minus_strand_threeSS_splice_site_coverage_bed['RNA_strand'] = '-'
minus_strand_threeSS_splice_site_coverage_bed

junction_df_plus_strand = junction_df[junction_df['RNA_strand'] == '+']
junction_df_minus_strand = junction_df[junction_df['RNA_strand'] == '-']

junction_df_plus_strand = junction_df_plus_strand.merge(plus_strand_fiveSS_splice_site_coverage_bed, how='left')
junction_df_plus_strand.fiveSS_coverage.isna().sum()
junction_df_minus_strand = junction_df_minus_strand.merge(minus_strand_fiveSS_splice_site_coverage_bed, how='left')
junction_df_minus_strand.fiveSS_coverage.isna().sum()

junction_df_plus_strand = junction_df_plus_strand.merge(plus_strand_threeSS_splice_site_coverage_bed, how='left')
junction_df_plus_strand.threeSS_coverage.isna().sum()
junction_df_minus_strand = junction_df_minus_strand.merge(minus_strand_threeSS_splice_site_coverage_bed, how='left')
junction_df_minus_strand.threeSS_coverage.isna().sum()

junction_df_combined = pd.concat([junction_df_plus_strand, junction_df_minus_strand])
junction_df_combined['fiveSS_coverage'] = junction_df_combined['fiveSS_coverage'].fillna(0)
junction_df_combined['threeSS_coverage'] = junction_df_combined['threeSS_coverage'].fillna(0)

junction_df_combined['mean_5SS_3SS_coverage'] = (junction_df_combined['fiveSS_coverage'] + junction_df_combined['threeSS_coverage'])/2

junction_df_combined['splicing_efficiency'] = junction_df_combined['COMPASS_counts'] / (junction_df_combined['COMPASS_counts'] + junction_df_combined['mean_5SS_3SS_coverage'])

junction_df_combined.to_csv(splice_junction_counts_outfilename, sep='\t', index=False)