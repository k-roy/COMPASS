#install.packages("tidyverse")
#install.packages("RSvgDevice")
library(tidyverse)
library(RSvgDevice)
library(cowplot)

DIR <- '/Users/kevinroy/Google_Drive/splice_analysis_results/' # '/Volumes/SPxDrive/splice_junction_analysis/' #  
FIGURE_DIR <- paste(DIR, 'figures/', sep = '') 
tbl <- read_tsv(paste(DIR,'combined_prp18_splice_junction_stats.txt', sep = '') )

PSEUDO_COUNT <- 1
PSEUDO_FRACTION <- 0.000001

tbl$upf1_prp18_mean_SE <- tbl %>% select('upf1_prp18:SE', 'BY_upf1del_prp18del_2:SE', 'BY_upf1del_prp18del_3:SE') %>% apply(1, mean, na.rm = TRUE)
tbl$upf1_prp18_sd_SE <- tbl %>% select('upf1_prp18:SE', 'BY_upf1del_prp18del_2:SE', 'BY_upf1del_prp18del_3:SE') %>% apply(1, sd, na.rm = TRUE)
tbl$upf1_mean_SE <- tbl %>% select('upf1:SE', 'upf1_8.17.15:SE', 'upf1_8.14.15:SE') %>% apply(1, mean, na.rm = TRUE)
tbl$upf1_sd_SE <- tbl %>% select('upf1:SE', 'upf1_8.17.15:SE', 'upf1_8.14.15:SE') %>% apply(1, sd, na.rm = TRUE)
tbl$upf1_prp18_mean_reads <- tbl %>% select('upf1_prp18:unique_junction_reads_rpm', 'BY_upf1del_prp18del_2:unique_junction_reads_rpm', 'BY_upf1del_prp18del_3:unique_junction_reads_rpm') %>% apply(1, mean, na.rm = TRUE)
tbl$upf1_mean_reads <- tbl %>% select('upf1:unique_junction_reads_rpm', 'upf1_8.17.15:unique_junction_reads_rpm', 'upf1_8.14.15:unique_junction_reads_rpm') %>% apply(1, mean, na.rm = TRUE)
tbl$upf1_prp18_mean_FAnS <- tbl %>% select('upf1_prp18:fraction_of_annotated_splicing', 'BY_upf1del_prp18del_2:fraction_of_annotated_splicing', 'BY_upf1del_prp18del_3:fraction_of_annotated_splicing') %>% apply(1, mean, na.rm = TRUE)
tbl$upf1_mean_FAnS <- tbl %>% select('upf1:fraction_of_annotated_splicing', 'upf1_8.17.15:fraction_of_annotated_splicing', 'upf1_8.14.15:fraction_of_annotated_splicing') %>% apply(1, mean, na.rm = TRUE)
tbl$upf1_prp18_annotated_mean_reads <- tbl %>% select('upf1_prp18:annotated_unique_junction_reads_rpm', 'BY_upf1del_prp18del_2:annotated_unique_junction_reads_rpm', 'BY_upf1del_prp18del_3:annotated_unique_junction_reads_rpm') %>% apply(1, mean, na.rm = TRUE)
tbl$upf1_annotated_mean_reads <- tbl %>% select('upf1:annotated_unique_junction_reads_rpm', 'upf1_8.17.15:annotated_unique_junction_reads_rpm', 'upf1_8.14.15:annotated_unique_junction_reads_rpm') %>% apply(1, mean, na.rm = TRUE)
tbl$upf1_prp18_mean_intron_reads <- tbl %>% select('upf1_prp18:intron_reads_rpm', 'BY_upf1del_prp18del_2:intron_reads_rpm', 'BY_upf1del_prp18del_3:intron_reads_rpm') %>% apply(1, mean, na.rm = TRUE)
tbl$upf1_mean_intron_reads <- tbl %>% select('upf1:intron_reads_rpm', 'upf1_8.17.15:intron_reads_rpm', 'upf1_8.14.15:intron_reads_rpm') %>% apply(1, mean, na.rm = TRUE)

tbl$upf1_annotated_mean_SE <- tbl %>% select('upf1:annotated_SE', 'upf1_8.17.15:annotated_SE', 'upf1_8.14.15:annotated_SE') %>% apply(1, mean, na.rm = TRUE)
tbl$upf1_annotated_prp18_mean_SE <- tbl %>% select('upf1_prp18:annotated_SE', 'BY_upf1del_prp18del_2:annotated_SE', 'BY_upf1del_prp18del_3:annotated_SE') %>% apply(1, mean, na.rm = TRUE)

tbl <- tbl %>% 
  mutate(upf1_prp18_annotated_mean_reads_with_pseudo_count = upf1_prp18_annotated_mean_reads + PSEUDO_COUNT,
         upf1_prp18_mean_reads_with_pseudo_count = upf1_prp18_mean_reads + PSEUDO_COUNT,
         upf1_prp18_mean_intron_reads_with_pseudo_count = upf1_prp18_mean_intron_reads + PSEUDO_COUNT,
         upf1_prp18_mean_SE_with_pseudo_count = upf1_prp18_mean_reads_with_pseudo_count / (upf1_prp18_mean_reads_with_pseudo_count + upf1_prp18_mean_intron_reads_with_pseudo_count),
         upf1_prp18_mean_FAnS_with_pseudo_count = upf1_prp18_mean_reads_with_pseudo_count / upf1_prp18_annotated_mean_reads_with_pseudo_count,
         
         upf1_annotated_mean_reads_with_pseudo_count = upf1_annotated_mean_reads + PSEUDO_COUNT,
         upf1_mean_reads_with_pseudo_count = upf1_mean_reads + PSEUDO_COUNT,
         upf1_mean_intron_reads_with_pseudo_count = upf1_mean_intron_reads+ PSEUDO_COUNT,
         upf1_mean_SE_with_pseudo_count = upf1_mean_reads_with_pseudo_count / (upf1_mean_reads_with_pseudo_count + upf1_mean_intron_reads_with_pseudo_count),
         upf1_mean_FAnS_with_pseudo_count = upf1_mean_reads_with_pseudo_count / upf1_annotated_mean_reads_with_pseudo_count,
         
         alt_minus_ann_U_score = alternative_U_score - annotated_U_score,
         
         prp18_effect_SE = upf1_prp18_mean_SE / upf1_mean_SE,
         prp18_effect_FAnS = upf1_prp18_mean_FAnS / upf1_mean_FAnS,
         prp18_effect_SE_with_pseudo_counts = upf1_prp18_mean_SE_with_pseudo_count / upf1_mean_SE_with_pseudo_count,
         prp18_effect_FAnS_with_pseudo_counts = upf1_prp18_mean_FAnS_with_pseudo_count / upf1_mean_FAnS_with_pseudo_count,
         fraction_of_gapped_reads_without_junction_proximal_mismatches = total_reads_without_mismatch_near_junction / total_gapped_reads,
         annotated = ifelse(alt_5SS == FALSE & alt_3SS == FALSE, TRUE, FALSE),
         threeSS_type_motif = ifelse(threeSS_type %in% c('CAG', 'TAG'), 'YAG', ifelse(threeSS_type %in% c('AC', 'other'), 'other', threeSS_type)),
         threeSS_type_motif_factor = factor(threeSS_type_motif, levels=c('YAG', 'AAG', 'GAG', 'BG', 'HAT', 'other')),
         threeSS_type_factor = factor(threeSS_type, levels=c('CAG', 'TAG','AAG', 'GAG', 'BG', 'HAT', 'AC', 'other')),
         fiveSS_type_factor = factor(fiveSS_type, levels=c('GTATGT', 'GT', 'non-GT')),
         alt_5SS_alt_3SS = paste(alt_5SS,alt_3SS), 
         annotated_junction = paste(annotated_chromosome, annotated_start, annotated_end, annotated_strand, sep = '_' ),
         threeSS_type_annotated = ifelse( annotated & threeSS_type_motif == 'YAG',  'annotated_YAG',

                                          ifelse( alt_3SS & threeSS_type_motif == 'YAG',  'alternative_YAG',
                                                  ifelse( alt_5SS & !alt_3SS & threeSS_type_motif == 'YAG',  'alternative_5SS_YAG',
                                                  threeSS_type_motif
                                          ) ) ),
         threeSS_type_annotated_factor =  factor(threeSS_type_annotated, levels=c('annotated_YAG', 'alternative_YAG', 'AAG', 'GAG', 'BG', 'HAT', 'other', 'alternative_5SS_YAG')),
         upf1_prp18_mean_FAnS_with_pseudo_fraction = upf1_prp18_mean_FAnS_with_pseudo_count + PSEUDO_FRACTION,
         upf1_mean_FAnS_with_pseudo_fraction =  upf1_mean_FAnS_with_pseudo_count + PSEUDO_FRACTION,
         upf1_prp18_mean_SE_with_pseudo_fraction =  upf1_prp18_mean_SE + PSEUDO_FRACTION,
         upf1_mean_SE_with_pseudo_fraction =  upf1_mean_SE + PSEUDO_FRACTION,
         log2_BP_3SS_dist = log(BP_3SS_dist, 2), 
         log2_BP_5SS_dist = log(BP_5SS_dist, 2), 
         log2_annotated_BP_3SS_dist = log(annotated_BP_3SS_dist, 2),
         log10_prp18_effect_SE_with_pseudo_counts = log(prp18_effect_SE_with_pseudo_counts,10),
         log10_upf1_mean_SE = log(upf1_mean_SE_with_pseudo_fraction, 10), 
         log10_upf1_prp18_mean_SE = log(upf1_prp18_mean_SE_with_pseudo_fraction, 10), 
         log10_prp18_effect_FAnS_with_pseudo_counts = log(prp18_effect_FAnS_with_pseudo_counts, 10),
         log10_upf1_mean_FAnS_with_pseudo_fraction = log(upf1_mean_FAnS_with_pseudo_fraction, 10),
         log10_upf1_prp18_mean_FAnS_with_pseudo_fraction = log(upf1_prp18_mean_FAnS_with_pseudo_fraction, 10) )
    

summary(tbl$threeSS_type_motif)
tbl$threeSS_type_motif

filtered_tbl <- tbl %>% 
  ## parameters below guided by analysis in Fig S3
  filter( fraction_of_gapped_reads_without_junction_proximal_mismatches > 0.5, overhang_perfect_matches_upstream_junction > 10 & overhang_perfect_matches_downstream_junction > 40 ) %>%
  arrange(desc( prp18_effect_SE_with_pseudo_counts ) ) %>%
  ## filter below guided by manual inspection of systematic sequencing errors at annotated splice junctions with homopolymers downstream of the 3'SS
  ## this filters out alt 5'SS which are very close to the annotated, and likely caused by sequencing error
  filter( ! ( (0 < annotated_to_alternative_5SS_dist &  annotated_to_alternative_5SS_dist < 4 & 0 < annotated_to_alternative_3SS_dist & annotated_to_alternative_3SS_dist < 4) & ( substr(downstream_3SS,1,4 ) == 'AAAA' | substr(downstream_3SS,1,4 ) == 'TTTT' | substr(downstream_3SS,1,4 ) == 'CCCC' | substr(downstream_3SS,1,4 ) == 'GGGG'  ) ) ) %>%
  filter( ! ( (0 > annotated_to_alternative_5SS_dist &  annotated_to_alternative_5SS_dist >-4 & 0 > annotated_to_alternative_3SS_dist & annotated_to_alternative_3SS_dist >-4) & ( substr(downstream_3SS,1,4 ) == 'AAAA' | substr(downstream_3SS,1,4 ) == 'TTTT' | substr(downstream_3SS,1,4 ) == 'CCCC' | substr(downstream_3SS,1,4 ) == 'GGGG'  ) ) ) %>%
  ## later, an anti-join with the full table can be used to assign 'passing_junction_filter = FALSE'
  mutate(passing_junction_filter = TRUE) 

filtered_tbl %>% filter(common_gene_name == 'PHO85') %>% select(annotated_junction, junction_ID_name, total_gapped_reads, BP_5SS_dist, BP_3SS_dist, fiveSS_seq, annotated_to_alternative_3SS_dist, threeSS_seq, annotated,  upf1_mean_reads, upf1_prp18_mean_reads, upf1_mean_FAnS, upf1_prp18_mean_FAnS, upf1_mean_SE, upf1_prp18_mean_SE) %>% arrange(BP_3SS_dist )
filtered_tbl %>% filter(annotated)
filtered_tbl %>% filter(BP_3SS_dist < 11, total_gapped_reads > 5 ) %>% select(common_gene_name, junction_ID_name, fiveSS_seq,  upstream_3SS, threeSS_seq, downstream_3SS, annotated_to_alternative_5SS_dist, annotated_to_alternative_3SS_dist, total_gapped_reads, BP_3SS_dist) %>% arrange(total_gapped_reads ) %>% print(n=100) 
                                                                    #  , upf1_mean_SE, upf1_prp18_mean_SE, upf1_prp18_mean_FAnS,
                                                                   #   BP_5SS_dist, BP_3SS_dist, annotated,  upf1_mean_reads, upf1_prp18_mean_reads, upf1_mean_FAnS, upf1_prp18_mean_FAnS, upf1_mean_SE, upf1_prp18_mean_SE) 
filtered_tbl %>% filter(annotated & fraction_of_gapped_reads_without_junction_proximal_mismatches < .95) %>% select(common_gene_name, junction_ID_name, upstream_5SS, downstream_3SS, fraction_of_gapped_reads_without_junction_proximal_mismatches, total_gapped_reads, BP_5SS_dist, BP_3SS_dist, fiveSS_seq, annotated_to_alternative_3SS_dist, threeSS_seq, annotated,  upf1_mean_reads, upf1_prp18_mean_reads, upf1_mean_FAnS, upf1_prp18_mean_FAnS, upf1_mean_SE, upf1_prp18_mean_SE) %>% arrange(BP_3SS_dist )
filtered_tbl_min_6_reads %>% filter(threeSS_seq == 'TAG' & annotated_to_alternative_3SS_dist < -40 & !alt_5SS) %>% select(common_gene_name, annotated_to_alternative_3SS_dist,upf1_prp18_mean_FAnS, upf1_mean_FAnS,  threeSS_seq, 
                                                                 BP_3SS_dist, total_gapped_reads, upf1_mean_SE, upf1_prp18_mean_SE, upf1_prp18_mean_FAnS, junction_ID_name, upstream_5SS, downstream_3SS, fraction_of_gapped_reads_without_junction_proximal_mismatches, total_gapped_reads, BP_5SS_dist, BP_3SS_dist, fiveSS_seq,  annotated,  upf1_mean_reads, upf1_prp18_mean_reads, upf1_mean_FAnS, upf1_prp18_mean_FAnS) %>% arrange(BP_3SS_dist )

filtered_tbl_min_6_reads %>% filter(annotated & threeSS_seq == 'AAG') %>% select(common_gene_name, annotated_to_alternative_3SS_dist, threeSS_seq, 
                                                                 BP_3SS_dist, total_gapped_reads,  upf1_prp18_mean_FAnS, upf1_mean_FAnS, upf1_prp18_mean_SE, junction_ID_name, upstream_5SS, downstream_3SS, fraction_of_gapped_reads_without_junction_proximal_mismatches, total_gapped_reads, BP_5SS_dist, BP_3SS_dist, fiveSS_seq,  annotated,  upf1_mean_reads, upf1_prp18_mean_reads, upf1_mean_FAnS, upf1_prp18_mean_FAnS) %>% arrange(BP_3SS_dist )

filtered_tbl %>% filter( is.na(threeSS_type_annotated_factor) )  %>% select(common_gene_name, annotated_to_alternative_3SS_dist, threeSS_seq, 
                                                                            BP_3SS_dist, total_gapped_reads, upf1_mean_SE, upf1_prp18_mean_SE, upf1_prp18_mean_FAnS, junction_ID_name, upstream_5SS, downstream_3SS, fraction_of_gapped_reads_without_junction_proximal_mismatches, total_gapped_reads, BP_5SS_dist, BP_3SS_dist, fiveSS_seq,  annotated,  upf1_mean_reads, upf1_prp18_mean_reads, upf1_mean_FAnS, upf1_prp18_mean_FAnS) %>% arrange(BP_3SS_dist )

filtered_tbl_min_6_reads %>% filter(threeSS_type_motif == 'HAT' |  threeSS_type_motif == 'other') %>% select(common_gene_name, annotated_to_alternative_3SS_dist, threeSS_seq,BP_3SS_dist, total_gapped_reads,  upf1_prp18_mean_FAnS, upf1_mean_FAnS, upf1_prp18_mean_SE, junction_ID_name, upstream_5SS, downstream_3SS, fraction_of_gapped_reads_without_junction_proximal_mismatches, total_gapped_reads, BP_5SS_dist, BP_3SS_dist, fiveSS_seq,  annotated,  upf1_mean_reads, upf1_prp18_mean_reads, upf1_mean_FAnS, upf1_prp18_mean_FAnS) %>% arrange(desc(upf1_prp18_mean_FAnS ))
filtered_tbl_min_6_reads$threeSS_type_motif

BP_3SS_nt_acessibility_tbl <- read_tsv(paste(DIR,'BP_3SS_nt_effective_distances_accessibility_with_common_gene_name.txt', sep = '') )
filtered_tbl <- filtered_tbl %>% left_join(BP_3SS_nt_acessibility_tbl, by = c("junction_ID_name",  "common_gene_name", "threeSS_seq"))


all_trint_randomized <- read_tsv(paste(DIR,'all_trinucleotides_between_bp_and_annotated_3SS_plus_5nt_downstream_randomized.txt', sep = '') )
all_trint_randomized <- all_trint_randomized %>% mutate(threeSS_type_factor = factor(motif, levels=c('CAG', 'TAG','AAG', 'GAG', 'BG', 'HAT', 'AC', 'other')) )
all_trint_randomized <- all_trint_randomized %>% mutate(annotated_junction = paste(annotated_chromosome, annotated_start, annotated_end, annotated_strand, sep = '_' ) )

all_trint <- read_tsv(paste(DIR,'all_trinucleotides_between_bp_and_annotated_3SS_plus_50nt_downstream.txt', sep = '') )
all_trint <- all_trint %>% mutate(threeSS_type_factor = factor(motif, levels=c('CAG', 'TAG','AAG', 'GAG', 'BG', 'HAT', 'AC', 'other')) )
all_trint$threeSS_type_factor
all_trint <- all_trint %>% mutate(threeSS_type_annotated = ifelse( annotated_to_alternative_3SS_dist == 0 & seq %in% c('CAG', 'TAG'),  'annotated_YAG',
                                                                   
                                                                   ifelse( annotated_to_alternative_3SS_dist != 0 & seq %in% c('CAG', 'TAG'),  'alternative_YAG',
                                                                           motif
                                                                           ) ) ) 
all_trint <- all_trint %>% mutate(threeSS_type_annotated_factor =  factor(threeSS_type_annotated, levels=c('annotated_YAG', 'alternative_YAG', 'AAG', 'GAG', 'BG', 'HAT', 'AC', 'other')) )

all_trint <- all_trint %>% mutate(annotated_junction = paste(annotated_chromosome, annotated_start, annotated_end, annotated_strand, sep = '_' ) )
#all_trint <- all_trint %>% left_join(filtered_tbl %>% select(annotated_junction, RPG))
all_trint$threeSS_type_annotated
all_trint$annotated_to_alternative_3SS_dist

filtered_tbl <- filtered_tbl %>%
  mutate(median_effective_BP_3SS_dist = median(effective_BP_3SS_dist), 
         mean_effective_BP_3SS_dist = mean(effective_BP_3SS_dist), 
         log2_effective_BP_3SS_dist = log(effective_BP_3SS_dist, 2), 
         RNA_structure_present = ifelse( (BP_3SS_dist - effective_BP_3SS_dist) > 5, TRUE,  FALSE),
         structured_bp =  BP_3SS_dist - effective_BP_3SS_dist) 

filtered_tbl_min_6_reads <- filtered_tbl %>% 
  filter( total_gapped_reads > 5) 

filtered_tbl_min_6_reads_min_one_thousandth_FAnS <- filtered_tbl_min_6_reads %>% 
  filter(upf1_prp18_mean_FAnS_with_pseudo_fraction >= 0.001 | upf1_mean_FAnS_with_pseudo_fraction >= 0.001)

filtered_tbl_min_6_reads_with_no_structure <- filtered_tbl_min_6_reads %>% filter(RNA_structure_present == FALSE)
filtered_tbl_min_6_reads_with_structure <- filtered_tbl_min_6_reads %>% filter(RNA_structure_present == TRUE)

annotated_only <- filtered_tbl %>% filter(alt_3SS == FALSE & alt_5SS == FALSE) 
alternative_only <- filtered_tbl %>% filter(alt_5SS == TRUE | alt_3SS == TRUE) 
alternative_3SS_only <- filtered_tbl %>% filter(alt_3SS == TRUE & alt_5SS == FALSE) 
alternative_5SS_only <- filtered_tbl %>% filter(alt_5SS == TRUE & alt_3SS == FALSE)
alternative_5SS_3SS_only <- filtered_tbl %>% filter(alt_5SS == TRUE & alt_3SS == TRUE) 

alternative_3SS_only$threeSS_type
alternative_3SS_only %>% filter(total_gapped_reads > 5) %>% group_by(threeSS_type) %>% count()

filtered_tbl %>% filter(annotated_to_alternative_3SS_dist > 0, threeSS_seq == 'GAG', log10_prp18_effect_FAnS_with_pseudo_counts > 1) %>% select(common_gene_name, threeSS_seq, annotated_to_alternative_5SS_dist, annotated_to_alternative_3SS_dist, upf1_prp18_mean_FAnS_with_pseudo_fraction, upf1_mean_FAnS_with_pseudo_fraction, log10_prp18_effect_FAnS_with_pseudo_counts) %>%
  arrange(desc(upf1_prp18_mean_FAnS_with_pseudo_fraction)) %>% print(n=20)

filtered_tbl %>% filter(common_gene_name == 'MTR2') %>% select(common_gene_name, threeSS_seq, upstream_5SS, annotated_best_BP_sequence, fiveSS_seq, upstream_3SS, annotated_to_alternative_5SS_dist, annotated_to_alternative_3SS_dist, upf1_prp18_mean_SE_with_pseudo_fraction, upf1_mean_SE_with_pseudo_fraction, log10_prp18_effect_FAnS_with_pseudo_counts) %>%
  arrange(desc(upf1_prp18_mean_SE_with_pseudo_fraction)) %>% print(n=20)

## should also sum the total FAnS for each annotated intron  (i.e. generate Splicing "Error" Frequency metrics)
alt_junctions_per_intron <- filtered_tbl_min_6_reads  %>% group_by(annotated_junction, systematic_gene_name, common_gene_name, RPG,  alt_5SS_alt_3SS ) %>% tally()

alt_junctions_per_intron_min_one_thousandth_FAnS <- filtered_tbl_min_6_reads %>% filter(upf1_mean_FAnS >= 0.001 | upf1_prp18_mean_FAnS >= 0.001 ) %>% group_by(annotated_junction, systematic_gene_name, common_gene_name, RPG,  alt_5SS_alt_3SS ) %>% tally()

alt_junctions_per_intron <- alt_junctions_per_intron %>% spread(alt_5SS_alt_3SS, n, fill = 0) %>% 
  mutate('total number of alternative 3SS junctions' = sum(`FALSE TRUE` )) %>%
  mutate('total number of alternative 5SS junctions' = sum(`TRUE FALSE`)) %>%
  mutate('total number of alternative 5SS & 3SS junctions' = sum(`TRUE TRUE` )) %>%
  rename('annotated junction identified (1=Yes, 0=No)' = `FALSE FALSE`, 'alt. 3SS' = `FALSE TRUE`,'alt. 5SS' = `TRUE FALSE`, 'alt. 5SS & alt. 3SS' = `TRUE TRUE`   ) 

alt_junctions_per_intron_min_one_thousandth_FAnS <- alt_junctions_per_intron_min_one_thousandth_FAnS %>% spread(alt_5SS_alt_3SS, n, fill = 0) %>% 
  mutate('total number of alternative 3SS junctions min_one_thousandth_FAnS' = sum(`FALSE TRUE` )) %>%
  mutate('total number of alternative 5SS junctions min_one_thousandth_FAnS' = sum(`TRUE FALSE`)) %>%
  mutate('total number of alternative 5SS & 3SS junctions min_one_thousandth_FAnS' = sum(`TRUE TRUE` )) %>%
  rename('alt. 3SS min_one_thousandth_FAnS' = `FALSE TRUE`,'alt. 5SS min_one_thousandth_FAnS' = `TRUE FALSE`, 'alt. 5SS & alt. 3SS min_one_thousandth_FAnS' = `TRUE TRUE`   ) 

combined_alt_tally <- left_join(alt_junctions_per_intron, alt_junctions_per_intron_min_one_thousandth_FAnS)

tallied_upf1_mean_FAnS <- alternative_3SS_only %>% group_by(annotated_junction) %>% tally(upf1_mean_FAnS)  %>% rename( total_upf1_FAnS = n)
tallied_upf1_prp18_mean_FAnS <- alternative_3SS_only %>% group_by(annotated_junction) %>% tally(upf1_prp18_mean_FAnS)  %>% rename( total_upf1_prp18_FAnS = n)

combined_alt_tally <- combined_alt_tally %>% left_join(tallied_upf1_mean_FAnS ) %>% left_join(tallied_upf1_prp18_mean_FAnS)

combined_alt_tally %>% arrange(desc(total_upf1_prp18_FAnS)) %>% write_tsv(paste(DIR,'number_of_alt_junctions_per_annotated_intron.txt', sep = ''))

# total_upf1_prp18_FAnS == 0 |
combined_alt_tally <- combined_alt_tally %>% mutate( total_upf1_prp18_FAnS = ifelse( is.na(total_upf1_prp18_FAnS), 0, total_upf1_prp18_FAnS ) )
combined_alt_tally <- combined_alt_tally %>% mutate( total_upf1_FAnS = ifelse( is.na(total_upf1_FAnS), 0, total_upf1_FAnS ) )

combined_alt_tally %>% mutate(junction_ID_name =annotated_junction)%>% left_join(filtered_tbl %>% select(junction_ID_name,  upf1_mean_SE) ) %>% arrange(desc(total_upf1_prp18_FAnS)) %>% select(common_gene_name, total_upf1_FAnS, total_upf1_prp18_FAnS, upf1_mean_SE) %>% print(n=60)

filtered_tbl$upf1_mean_SE
combined_alt_tally_gathered <- combined_alt_tally %>% gather( key = strain, value = total_FAnS, total_upf1_prp18_FAnS,  total_upf1_FAnS  )
combined_alt_tally_gathered <- combined_alt_tally_gathered %>% mutate(RPG_label = ifelse(RPG, paste('RPG',strain), paste('non-RPG',strain) ) )
combined_alt_tally_gathered$RPG_label
combined_alt_tally_gathered %>% filter(common_gene_name == 'NYV1') %>% select(total_FAnS, RPG_label)
# value = measurement,
# Sepal.Length, Sepal.Width, Petal.Length, Petal.Width)
combined_alt_tally$total_FAnS
p1 <- ggplot(combined_alt_tally, aes(total_upf1_prp18_FAnS, colour = c(RPG) ) ) + geom_density() + scale_x_log10()  # + facet_grid(strain ~ .)
p2 <- ggplot(combined_alt_tally, aes(total_upf1_FAnS, colour = c(RPG) ) ) + geom_density() + scale_x_log10()  # + facet_grid(strain ~ .)
ggplot(combined_alt_tally_gathered, aes(total_FAnS, colour = c(RPG_label) ) ) + geom_density( adjust = 1 ) + scale_x_log10()
filename <- (paste(FIGURE_DIR,'alt_3SS_distributions_of_total_FAnS_per_annotated_intron.eps', sep = '') )
ggsave(filename, width = 8, height = 4)

plot_grid(p1, p2, nrow =2 )

sum( combined_alt_tally$total_upf1_prp18_FAnS )
sum( combined_alt_tally$total_upf1_FAnS )

dev.off()
combined_alt_tally_RPG <- combined_alt_tally %>% filter(RPG)
combined_alt_tally_non_RPG <- combined_alt_tally %>% filter(!RPG)

max_value_index <-which.max(density(combined_alt_tally_RPG$total_upf1_prp18_FAnS)$y)
density(combined_alt_tally_RPG$total_upf1_prp18_FAnS)$x[max_value_index]

max_value_index <-which.max(density(combined_alt_tally_RPG$total_upf1_FAnS)$y)
density(combined_alt_tally_RPG$total_upf1_FAnS)$x[max_value_index]

combined_alt_tally$total_upf1_prp18_FAnS

max_value_index <-which.max(density(combined_alt_tally_non_RPG$total_upf1_prp18_FAnS)$y)
density(combined_alt_tally_non_RPG$total_upf1_prp18_FAnS)$x[max_value_index]

max_value_index <-which.max(density(combined_alt_tally_non_RPG$total_upf1_FAnS)$y)
density(combined_alt_tally_non_RPG$total_upf1_FAnS)$x[max_value_index]

filename <- (paste(FIGURE_DIR,'upf1_prp18_total_FAnS.eps', sep = '') )
ggsave(filename, width = 4, height = 4)


## check features of junctions removed by filter
anti_filtered_tbl <- anti_join(tbl, filtered_tbl)
anti_filtered_tbl <- anti_filtered_tbl %>% 
  mutate(passing_junction_filter = FALSE) 

filtered_and_nonfiltered <- bind_rows(filtered_tbl, anti_filtered_tbl )

filtered_tbl$structured_bp
## WRITE OUTFILES
filtered_tbl %>% filter(common_gene_name == 'HAC1') %>% write_tsv(paste(DIR,'HAC1_filtered_junctions.txt', sep = '') ) 
tbl %>% write_tsv(paste(DIR,'unfiltered_junctions.txt', sep = '') ) 

filtered_tbl %>% write_tsv(paste(DIR,'filtered_junctions.txt', sep = '') ) 

filtered_tbl %>% 
  select(common_gene_name, junction_ID_name, threeSS_seq, annotated, alt_5SS, alt_3SS, accessibility, BP_3SS_dist, effective_BP_3SS_dist, 
         prp18_effect_SE, prp18_effect_FAnS,  prp18_effect_SE_with_pseudo_counts, prp18_effect_FAnS_with_pseudo_counts, 
         total_gapped_reads, upf1_prp18_mean_SE, upf1_mean_SE,upf1_prp18_mean_SE_with_pseudo_count, upf1_mean_SE_with_pseudo_count,
         upf1_prp18_mean_FAnS, upf1_mean_FAnS, upf1_prp18_mean_FAnS_with_pseudo_count, upf1_mean_FAnS_with_pseudo_count,
         BP_3SS_seq_with_5_bp_DS, BP_3SS_fold_with_5_bp_DS) %>% 
  arrange(effective_BP_3SS_dist) %>% 
  write_csv(paste(DIR,'prp18_SE_FAnS_summary.csv', sep = '') )


sum(as.logical(tbl$alt_5SS) )
sum(as.logical(tbl$alt_3SS) )
summary( as.factor(tbl$fiveSS_type) )
summary( as.factor(tbl$threeSS_type) )
tbl %>% group_by(alt_5SS_alt_3SS) %>% tally()
tbl %>% group_by(threeSS_type) %>% tally()
tbl %>% group_by(fiveSS_type) %>% tally()


sum(as.logical(filtered_tbl$alt_5SS) )
sum(as.logical(filtered_tbl$alt_3SS) )
summary( as.factor(filtered_tbl$fiveSS_type) )
summary( as.factor(filtered_tbl$threeSS_type) )

sum(as.logical(filtered_tbl_min_6_reads$alt_5SS) )
sum(as.logical(filtered_tbl_min_6_reads$alt_3SS) )
summary( as.factor(filtered_tbl_min_6_reads$fiveSS_type) )
summary( as.factor(filtered_tbl_min_6_reads$threeSS_type) )

filtered_tbl_min_6_reads_min_one_thousandth_FAnS %>% 
  filter(alt_5SS_alt_3SS != "FALSE FALSE") %>%
  count()

## stats for supp tables
sum(as.logical(filtered_tbl_min_6_reads_min_one_thousandth_FAnS$alt_5SS) )
sum(as.logical(filtered_tbl_min_6_reads_min_one_thousandth_FAnS$alt_3SS) )
summary( as.factor(filtered_tbl_min_6_reads_min_one_thousandth_FAnS$fiveSS_type) )
summary( as.factor(filtered_tbl_min_6_reads_min_one_thousandth_FAnS$threeSS_type) )
tbl %>% count()
filtered_tbl %>% count()
anti_filtered_tbl %>% count()
filtered_tbl_min_6_reads %>% count()
filtered_tbl_min_6_reads_min_one_thousandth_FAnS %>% count()

tbl %>% filter(coord_for_best_match_to_BP_consensus == annotated_coord_for_best_match_to_BP_consensus) %>% select(junction_ID_name)%>% count()
tbl %>% filter(coord_for_best_match_to_BP_consensus != annotated_coord_for_best_match_to_BP_consensus) %>% select(junction_ID_name)%>% count()
filtered_tbl %>% filter(coord_for_best_match_to_BP_consensus == annotated_coord_for_best_match_to_BP_consensus) %>% select(junction_ID_name) %>% count()
filtered_tbl %>% filter(coord_for_best_match_to_BP_consensus != annotated_coord_for_best_match_to_BP_consensus) %>% select(junction_ID_name)%>% count()
filtered_tbl_min_6_reads %>% filter(coord_for_best_match_to_BP_consensus == annotated_coord_for_best_match_to_BP_consensus) %>% select(junction_ID_name) %>% count()
filtered_tbl_min_6_reads %>% filter(coord_for_best_match_to_BP_consensus != annotated_coord_for_best_match_to_BP_consensus) %>% select(junction_ID_name) %>% count()
filtered_tbl_min_6_reads_min_one_thousandth_FAnS %>% filter(coord_for_best_match_to_BP_consensus == annotated_coord_for_best_match_to_BP_consensus) %>% select(junction_ID_name) %>% count()
filtered_tbl_min_6_reads_min_one_thousandth_FAnS %>% filter(coord_for_best_match_to_BP_consensus != annotated_coord_for_best_match_to_BP_consensus) %>% select(junction_ID_name) %>% count()
## distribution of effective vs linear distance for annotated and alternative
summary(filtered_tbl_min_6_reads$RNA_structure_present)
greater_45 <- annotated_only %>% filter(BP_3SS_dist > 45)
summary(greater_45$RNA_structure_present)
summary(annotated_only$RNA_structure_present)
summary(alternative_3SS_only$RNA_structure_present)
summary(alternative_5SS_only$RNA_structure_present)
summary(alternative_5SS_3SS_only$RNA_structure_present)
filtered_tbl_min_6_reads %>% group_by(alt_5SS_alt_3SS) %>% count()
filtered_tbl_min_6_reads_min_one_thousandth_FAnS %>% group_by(alt_5SS_alt_3SS) %>% count()


params_to_write <- filtered_tbl_min_6_reads %>%
  select(junction_ID_name, effective_BP_3SS_dist, accessibility, systematic_gene_name, common_gene_name,  chromosome, adjusted_start, adjusted_end, annotated_strand, annotated_to_alternative_5SS_dist, annotated_to_alternative_3SS_dist,  best_BP_sequence, BP_5SS_dist , BP_3SS_dist ,  upstream_5SS, fiveSS_seq, downstream_5SS, upstream_3SS,  threeSS_seq, downstream_3SS, annotated_upstream_5SS_seq, annotated_fiveSS_seq, annotated_downstream_5SS_seq, annotated_upstream_3SS_seq, annotated_threeSS_seq, annotated_downstream_3SS_seq, annotated_poly_U_count, annotated_poly_Y_count, prp18_effect_FAnS_with_pseudo_counts, prp18_effect_SE_with_pseudo_counts, total_gapped_reads, total_reads_without_mismatch_near_junction, total_reads_with_perfect_gapped_alignments, overhang_perfect_matches_upstream_junction, overhang_perfect_matches_downstream_junction, max_ambiguous_bases, total_BBMAP_and_STAR, total_BBMAP_only, total_STAR_only, total_BBMAP_intron_not_in_STAR, upf1_prp18_mean_reads, upf1_mean_reads, upf1_prp18_annotated_mean_reads, upf1_annotated_mean_reads, upf1_prp18_mean_intron_reads, upf1_mean_intron_reads,  upf1_prp18_mean_FAnS_with_pseudo_count, upf1_mean_FAnS_with_pseudo_count, upf1_prp18_mean_FAnS, upf1_mean_FAnS, fiveSS_type, threeSS_type, alt_5SS, alt_3SS) %>% 
  filter(upf1_prp18_mean_FAnS > 0.001 | upf1_mean_FAnS > 0.001) %>%
  arrange(desc(upf1_prp18_mean_FAnS))

params_to_write %>% filter(alt_3SS == TRUE) %>% write_tsv(paste(DIR,'combined_upf1_prp18_splicing_alt_3SS_filtered.tsv', sep = ''))
params_to_write %>% write_tsv(paste(DIR,'combined_upf1_prp18_splicing_all_filtered.tsv', sep = ''))
params_to_write %>% filter(threeSS_type == 'non-G') %>% write_tsv(paste(DIR,'combined_upf1_prp18_splicing_non-G_3SS_filtered.tsv', sep = ''))
params_to_write %>% filter(threeSS_type == 'YAG') %>% write_tsv(paste(DIR,'combined_upf1_prp18_splicing_YAG_3SS_filtered.tsv', sep = ''))
params_to_write %>% filter(threeSS_type == 'RAG')%>% write_tsv(paste(DIR,'combined_upf1_prp18_splicing_RAG_3SS_filtered.tsv', sep = ''))
params_to_write %>% filter(threeSS_type == 'BG') %>%  write_tsv(paste(DIR,'combined_upf1_prp18_splicing_BG_3SS_filtered.tsv', sep = ''))
params_to_write %>% filter(fiveSS_type == 'GTATGT') %>%  write_tsv(paste(DIR,'combined_upf1_prp18_splicing_GTATGT_5SS_filtered.tsv', sep = ''))
params_to_write %>% filter(fiveSS_type == 'GT') %>%  write_tsv(paste(DIR,'combined_upf1_prp18_splicing_GT_5SS_filtered.tsv', sep = ''))
params_to_write %>% filter(fiveSS_type == 'non-GT') %>%  write_tsv(paste(DIR,'combined_upf1_prp18_splicing_non-GT_5SS_filtered.tsv', sep = ''))

## horizontal lines drawn for a minimum of 6 reads above line
## Fig S3A
p1 <- ggplot(tbl, aes(fraction_of_gapped_reads_without_junction_proximal_mismatches, total_gapped_reads) ) + 
  geom_jitter(aes(colour = annotated), width = 0.01, height = 0.03, size = .5, alpha = 0.5)  + 
  scale_y_log10() +  geom_hline(yintercept = 5, linetype = 'dotted') 
my_svg(paste(FIGURE_DIR, 'fraction_of_gapped_reads_without_junction_proximal_mismatches.svg'), 8, 4) 
p1
dev.off()

## these annotated junctions have a high level of junction proximal mismatches
## both have polyA stretch downstream 3SS
## label these points in the figure
tbl %>% filter( fraction_of_gapped_reads_without_junction_proximal_mismatches < .9, annotated == TRUE) %>% select(common_gene_name, upstream_5SS, downstream_3SS, total_gapped_reads, total_reads_with_perfect_gapped_alignments)
tbl %>% filter( overhang_perfect_matches_upstream_junction < 70, annotated == TRUE) %>% select(common_gene_name, upstream_5SS, downstream_3SS, total_gapped_reads, total_reads_with_perfect_gapped_alignments)
tbl %>% filter( overhang_perfect_matches_downstream_junction < 85, annotated == TRUE) %>% select(common_gene_name, upstream_5SS, downstream_3SS, total_gapped_reads, total_reads_with_perfect_gapped_alignments, overhang_perfect_matches_downstream_junction)

## Fig S3B
p2 <- ggplot(tbl, aes(overhang_perfect_matches_upstream_junction, total_gapped_reads) ) + 
  geom_jitter(aes(colour = annotated), width = 0.001,height = 0.03, size = .5, alpha = 0.5)  + 
  scale_y_log10() +   scale_x_continuous(breaks = seq(  0, 105, 5))  +  geom_hline(yintercept = 5, linetype = 'dotted') 
my_svg(paste(FIGURE_DIR, 'overhang_perfect_matches_upstream_junction.eps'), 8, 4) 
p2
dev.off()

## Fig S3C
p3 <- ggplot(tbl, aes(overhang_perfect_matches_downstream_junction, total_gapped_reads) ) + 
  geom_jitter(aes(colour = annotated), width = 0.001,height = 0.03, size = 1, alpha = 0.5)  + 
  scale_y_log10() +   scale_x_continuous(breaks = seq(  0, 105, 5)) +  geom_hline(yintercept = 5, linetype = 'dotted') 
my_svg(paste(FIGURE_DIR, 'overhang_perfect_matches_downstream_junction.svg'), 8, 4) 
p3
dev.off()

## BEGIN Fig S4a
p1 <- ggplot(tbl, aes(overhang_perfect_matches_downstream_junction)) + 
  geom_density() + facet_grid(alt_5SS_alt_3SS ~ ., scales = "free_y") +
  scale_x_continuous(breaks = seq(  0, 105, 10)) + geom_vline(xintercept = 40, linetype = 'dotted') 

p2 <- ggplot(tbl, aes(overhang_perfect_matches_upstream_junction)) + 
  geom_density() + facet_grid(alt_5SS_alt_3SS ~ ., scales = "free_y") +
  scale_x_continuous(breaks = seq(  0, 105, 10)) + geom_vline(xintercept = 10, linetype = 'dotted') 

p3 <- ggplot(tbl, aes(overhang_perfect_matches_downstream_junction)) + 
  geom_density() + facet_grid(threeSS_type ~ ., scales = "free_y") +
  scale_x_continuous(breaks = seq(  0, 105, 10)) + geom_vline(xintercept = 40, linetype = 'dotted')

p4 <- ggplot(tbl, aes(overhang_perfect_matches_upstream_junction)) + 
  geom_density() + facet_grid(threeSS_type ~ ., scales = "free_y") +
  scale_x_continuous(breaks = seq(  0, 105, 10)) + geom_vline(xintercept = 10, linetype = 'dotted')

p5 <- ggplot(tbl, aes(overhang_perfect_matches_downstream_junction)) + 
  geom_density() + facet_grid(fiveSS_type ~ ., scales = "free_y") +
  scale_x_continuous(breaks = seq(  0, 105, 10)) + geom_vline(xintercept = 40, linetype = 'dotted')

p6 <- ggplot(tbl, aes(overhang_perfect_matches_upstream_junction)) + 
  geom_density() + facet_grid(fiveSS_type ~ ., scales = "free_y") +
  scale_x_continuous(breaks = seq(  0, 105, 10)) + geom_vline(xintercept = 10, linetype = 'dotted')

library(ggpubr)
ggarrange(p2, p1, p4, p3, p6, p5, nrow = 3, ncol = 2)
filename <- (paste(FIGURE_DIR,'all_density_plots_overhang_perfect_matches_flanking_junction_v2.eps', sep = '') )
ggsave(filename, width = 12, height = 16)

## below is actual order of plots used for Fig S4
ggarrange(p2, p4, p6, p1, p3, p5, nrow = 2, ncol = 3)
filename <- (paste(FIGURE_DIR,'all_density_plots_overhang_perfect_matches_flanking_junction.eps', sep = '') )
ggsave(filename, width = 16, height = 12)

## BEGIN Fig S4b


filtered_and_nonfiltered <- filtered_and_nonfiltered %>%
  mutate(read_max_50 = ifelse(total_gapped_reads > 50, 51,
                              ifelse(total_gapped_reads <= 50, total_gapped_reads, NA) ))
summary( filtered_and_nonfiltered$read_max_50 )

ggplot(data = filtered_and_nonfiltered, aes(read_max_50) )+ geom_histogram(binwidth = 1,  aes(fill = passing_junction_filter) ) + 
  scale_x_continuous(limits = c(0, 52) ) # + scale_y_log10()
filename <- (paste(FIGURE_DIR,'all_filtered_junctions_total_gapped_reads_uneven_bins.eps', sep = '') )
ggsave(filename, width = 8, height = 6)

## binned version below ofplot above used for actual figure
ggplot(data = filtered_and_nonfiltered, aes(read_max_50) )+ 
  stat_bin(breaks = c(0,1,2,3,4,5,6,7,8,9,10,20,30,40,50,60), color = 'black', aes(fill = passing_junction_filter) ) + 
  scale_x_continuous(limits = c(0, 61), breaks = seq(  0, 50, 10) )  
# scale_y_continuous(limits = c(0, 650), breaks = seq(  0, 650, 100) )
filename <- (paste(FIGURE_DIR,'all_filtered_junctions_total_gapped_reads_histogram.eps', sep = '') )
ggsave(filename, width = 8, height = 6)

## BEGIN Fig S4c
## plot how many alt 3SS junctions observed by at different levels of SE or FAnS
filtered_tbl <- filtered_tbl %>% mutate(max_FAnS = ifelse(upf1_prp18_mean_FAnS > upf1_mean_FAnS, upf1_prp18_mean_FAnS, upf1_mean_FAnS) )
filtered_tbl <- filtered_tbl %>% arrange(desc( max_FAnS ))
filtered_tbl <- filtered_tbl %>% group_by(alt_5SS_alt_3SS) %>%  mutate(ordered_num = row_number() )
filtered_tbl$ordered_num
filtered_tbl %>% filter(annotated == TRUE)
ggplot(filtered_tbl %>% filter(!annotated), aes(colour = alt_5SS_alt_3SS)) + geom_line(aes(ordered_num, max_FAnS) ) + scale_y_log10() + scale_x_continuous(breaks = seq(0,1500, 100) )
filename <- (paste(FIGURE_DIR,'number_of_junctions_by_FAnS_by_alt5SS_3SS.eps', sep = '') )
ggsave(filename, width = 8, height = 6)

filtered_tbl_min_6_reads <- filtered_tbl_min_6_reads %>% mutate(max_FAnS = ifelse(upf1_prp18_mean_FAnS > upf1_mean_FAnS, upf1_prp18_mean_FAnS, upf1_mean_FAnS) )
filtered_tbl_min_6_reads <- filtered_tbl_min_6_reads %>% arrange(desc( max_FAnS ))
filtered_tbl_min_6_reads <- filtered_tbl_min_6_reads %>% group_by(alt_5SS_alt_3SS) %>%  mutate(ordered_num = row_number() )
filtered_tbl_min_6_reads$ordered_num
ggplot(filtered_tbl_min_6_reads %>% filter(!annotated), aes(colour = alt_5SS_alt_3SS)) + geom_line(aes(ordered_num, max_FAnS) ) + scale_y_log10() + scale_x_continuous(breaks = seq(0,1500, 100) )
filename <- (paste(FIGURE_DIR,'number_of_junctions_by_FAnS_by_alt5SS_3SS_6_read_min.eps', sep = '') )
ggsave(filename, width = 8, height = 6)

## END Fig S4

## BEGIN Fig S4c
## how many alt junctions per annotated introns
## Fig S5
p1 <- ggplot(alt_junctions_per_intron, aes(`total number of alternative 3SS junctions`) ) +
  geom_histogram(binwidth = 1, color = 'black', aes(fill = RPG)) +
  scale_x_continuous(breaks = seq(0,20,1))

p2 <- ggplot(alt_junctions_per_intron_min_one_thousandth_FAnS, aes(`total number of alternative 3SS junctions min_one_thousandth_FAnS`) ) +
  geom_histogram(binwidth = 1, color = 'black', aes(fill = RPG)) +
  scale_x_continuous(breaks = seq(0,20,1))

p3 <- ggplot(alt_junctions_per_intron, aes(`total number of alternative 5SS junctions`) ) +
  geom_histogram(binwidth = 1, color = 'black', aes(fill = RPG)) +
  scale_x_continuous(breaks = seq(0,20,1))

p4 <- ggplot(alt_junctions_per_intron_min_one_thousandth_FAnS, aes(`total number of alternative 5SS junctions min_one_thousandth_FAnS`) ) +
  geom_histogram(binwidth = 1, color = 'black', aes(fill = RPG)) +
  scale_x_continuous(breaks = seq(0,20,1))

filename <- paste(FIGURE_DIR, 'num_alt_junctions_per_annotated_intron_min_6_reads.svg')
ggarrange(p1, p2, p3, p4)
ggsave(filename, width = 12, height = 6)

## END Fig S5
