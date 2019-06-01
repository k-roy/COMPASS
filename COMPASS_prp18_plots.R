library(GGally)
library(ggplot2)
library(ggpubr)
require(gridExtra)
library(cowplot)
library(svglite)
library(sfsmisc)
library(RColorBrewer)
library(ggsignif)
library(ggrepel)
theme_pubr()

HEIGHT <- 5
WIDTH <- 5
ADJUST_VALUE <- 1
SIZE <- 2
BOX_WIDTH <- .4
NUDGE_X <- -.4
w <- 6
h <- 4

ggplotRegression <- function (fit, additional_label) {
  
  require(ggplot2)
  
  ggplot(fit$model, aes_string(x = names(fit$model)[2], y = names(fit$model)[1])) + 
    geom_point() +
    stat_smooth(method = "lm", col = "red") +
    labs(title = paste(additional_label, 
                       #"Int=",signif(fit$coef[[1]],3 ), "R2=",signif(summary(fit)$adj.r.squared, 3),
                       "slope=",signif(fit$coef[[2]], 3),
                       "p=",signif(summary(fit)$coef[2,4], 3)))
}

DIR <- '/Users/kevinroy/Google_Drive/splice_analysis_results/' # '/Volumes/SPxDrive/splice_junction_analysis/' #  
FIGURE_DIR <- paste(DIR, 'figures/', sep = '') 

# Fig 1b
ggplot(annotated_only, aes(upf1_mean_SE, upf1_prp18_mean_SE) ) + geom_point( aes(colour = factor(RPG)), size = 2, alpha = .5 )  + 
  geom_abline(intercept = 0, slope = 1, linetype = 'dotted', alpha = .5) +
  geom_text_repel(aes(label=ifelse( upf1_prp18_mean_SE >=  upf1_mean_SE & upf1_prp18_mean_SE < .8,    as.character(common_gene_name),'')), nudge_y = .1 ) +
  geom_text_repel(aes(label=ifelse(  (upf1_prp18_mean_SE - upf1_mean_SE) < -.5,    as.character(common_gene_name), '' )), nudge_y = -.1) +
  geom_text_repel(aes(label=ifelse(  RPG & upf1_prp18_mean_SE < .5,    as.character(common_gene_name), '' )), nudge_y = -.1) 
ggsave(filename = paste(FIGURE_DIR, 'prp18_SE_annotated_only_RPGs.svg'), device = 'svg', width = 6, height = 4) 

# Fig 1c
ggplot(filtered_tbl, aes(upf1_mean_FAnS, upf1_prp18_mean_FAnS) ) + geom_point( aes(colour = factor(RPG)), size = 2, alpha = .5)   + 
  geom_abline(intercept = 0, slope = 1, linetype = 'dotted', alpha = .5) + facet_grid(. ~ alt_5SS_alt_3SS) +
  # geom_text_repel(aes(label=ifelse( upf1_prp18_mean_FAnS > 1,    as.character(paste(common_gene_name, threeSS_seq, sep = ' ')),'')), nudge_x = -1.5 ) +
  scale_y_log10() + scale_x_log10() 
ggsave(filename = paste(FIGURE_DIR, 'prp18_FAnS_annotated_only_RPGs_no_label.svg'), device = 'svg', width = 16, height = 4) 

## other way of plotting the above (shown for SE)
ggscatter(filtered_tbl, y = "upf1_prp18_mean_SE", x = "upf1_mean_SE", 
          size = 2, alpha = .2, 
          facet.by = c("alt_5SS_alt_3SS"), nrow = 1, color = c("RPG") ) + 
  geom_abline(slope = 1, linetype = 'dotted') +
  scale_x_log10(limits = c(-Inf, 1) ) + scale_y_log10(limits = c(-Inf, 1))
filename <- (paste(FIGURE_DIR,'SE_alternative_5SS_only_by_5SS_type_without_pseudocounts.svg', sep = '') )
ggsave(filename, width = 12, height = 4)

# Fig 1d
my_comparisons <- list( c(TRUE, FALSE) )
ggboxplot(annotated_only, x = "RPG", y = "log10_prp18_effect_SE_with_pseudo_counts",
          color = "RPG", width = .4, facet.by = c("alt_5SS_alt_3SS"), nrow = 1) +
  geom_jitter(width = 0.1, size = 1, alpha = 0.2)  + 
  geom_abline(intercept = 0, slope = 0, linetype = 'dashed', alpha = .3)  + 
  geom_text_repel(aes(label=ifelse( log10_prp18_effect_SE_with_pseudo_counts > 0 | log10_prp18_effect_SE_with_pseudo_counts < -.2,    as.character(paste(common_gene_name, threeSS_seq, sep = ' ')),'')), nudge_x = .5 ) +
  # geom_text(label = filtered_tbl$common_gene_name, nudge_x  = NUDGE_X, check_overlap = TRUE, size = 5) +
  stat_compare_means(comparisons = my_comparisons) # Add pairwise comparisons p-value
aggsave(filename = paste(FIGURE_DIR, 'prp18_effect_SE_faceted_RPGs.svg'), device = 'svg', width = 3, height = 6) 

# Fig 1e
my_comparisons <- list( c(TRUE, FALSE) )
summary(filtered_tbl$RPG)
ggboxplot(filtered_tbl, x = "RPG", y = "log10_prp18_effect_FAnS_with_pseudo_counts",
          color = "RPG", width = .4, facet.by = c("alt_5SS_alt_3SS"), nrow = 1) +
  geom_jitter(width = 0.1, size = 1, alpha = 0.1) + 
  geom_abline(intercept = 0, slope = 0, linetype = 'dashed', alpha = .5)  + 
  # geom_text(label = filtered_tbl$common_gene_name, nudge_x  = NUDGE_X, check_overlap = TRUE, size = 5) +
  stat_compare_means(comparisons = my_comparisons) # Add pairwise comparisons p-value
ggsave(filename = paste(DIR, 'prp18_effect_FAnS_faceted_RPGs.svg'), device = 'svg', width = 12, height = 6) 

## Fig 2A

genes_to_highlight <- c('PHO85', 'MUD1', 'SEC14', 'SPT14', 'MAF1', 'UBC12', 'BET4', 'YCL002C', 'NYV1', 'RAD14', 'RPL30')
ggscatter(alternative_3SS_only %>% filter(threeSS_type != 'other',threeSS_type != 'AC'), y = "upf1_prp18_mean_FAnS_with_pseudo_fraction", x = "upf1_mean_FAnS_with_pseudo_fraction", 
          size = 2, alpha = .3, 
          facet.by = c("threeSS_type_factor"), color = c("fiveSS_type") ) + 
  geom_text_repel(aes(label=ifelse(common_gene_name %in% genes_to_highlight & total_gapped_reads > 5 & annotated_to_alternative_3SS_dist < 0 & upf1_prp18_mean_FAnS_with_pseudo_fraction >= 0.1, as.character(paste(common_gene_name, threeSS_seq, sep = ' ')),'')) , nudge_y = 1) +
  geom_abline(slope = 1, linetype = 'dotted') +
  scale_x_log10(limits = c(-Inf, 10) ) + scale_y_log10(limits = c(-Inf, 10)) #
# geom_segment(aes(x=1e-03,xend=1e+01,y=1e-03,yend=1e-03))
# geom_hline( yintercept = 1e-03, linetype = 'dotted') + geom_vline(xintercept = 1e-03, linetype = 'dotted')
filename <- (paste(FIGURE_DIR,'all_introns_5SS_by_3SS_type_AF.svg', sep = '') )
ggsave(filename, width = 9, height = 6)

## Supp fig based on above
ggscatter(filtered_tbl %>% filter(threeSS_type != 'other'), y = "upf1_prp18_mean_FAnS_with_pseudo_fraction", x = "upf1_mean_FAnS_with_pseudo_fraction", 
          size = 2, alpha = .3, 
          facet.by = c( "fiveSS_type_factor","threeSS_type_factor"), color = c("alt_5SS_alt_3SS") ) + 
 # geom_text_repel(aes(label=ifelse(common_gene_name %in% genes_to_highlight & total_gapped_reads > 5 & annotated_to_alternative_3SS_dist < 0 & upf1_prp18_mean_FAnS_with_pseudo_fraction >= 0.1, as.character(paste(common_gene_name, threeSS_seq, sep = ' ')),'')) , nudge_y = 1) +
  geom_abline(slope = 1, linetype = 'dotted') +
  scale_x_log10(limits = c(-Inf, 10) ) + scale_y_log10(limits = c(-Inf, 10)) #
# geom_segment(aes(x=1e-03,xend=1e+01,y=1e-03,yend=1e-03))
# geom_hline( yintercept = 1e-03, linetype = 'dotted') + geom_vline(xintercept = 1e-03, linetype = 'dotted')
filename <- (paste(FIGURE_DIR,'all_introns_5SS_by_3SS_type_AF.svg', sep = '') )
ggsave(filename, width = 8, height = 8)

## Figure 2B
my_comparisons <- list( c("AAG", "GAG"), c("AAG", "BG"),c("AAG", "HAT") )
genes_to_highlight <- c('PHO85', 'MUD1', 'SEC14', 'SPT14', 'MAF1', 'UBC12', 'BET4', 'YCL002C', 'NYV1', 'RAD14')
alternative_3SS_only$threeSS_type_motif_factor
p2 <- ggboxplot(alternative_3SS_only %>% filter(threeSS_type_motif_factor != 'other', threeSS_type_motif_factor != 'AC'), x = "threeSS_type_motif_factor", y = "log10_prp18_effect_FAnS_with_pseudo_counts",
                color = "threeSS_type_motif_factor",  width = BOX_WIDTH)+ 
  geom_hline(yintercept = 0, linetype = 'dotted') + 
  geom_jitter(width = 0.15, size = 1, alpha = .15) +
  geom_text_repel(aes(label=ifelse(common_gene_name %in% genes_to_highlight & total_gapped_reads > 5 & annotated_to_alternative_3SS_dist < 0 & upf1_prp18_mean_FAnS_with_pseudo_fraction >= 0.1, as.character(common_gene_name),'')) , nudge_x = .5, size = 5, force = 1, fontface = 'italic') #+
 # stat_compare_means(comparisons = my_comparisons) # Add pairwise comparisons p-value stat_compare_means(label.y = 3.5)     # Add global p-value
p2 # , nudge_x = .2, size = 4 
ggsave(paste(FIGURE_DIR,'all_alternative_junctions_threeSS_type_by_alt_fraction_FAnS_prp18_effect.svg', sep = '') , width = 6, height = 6)

my_comparisons <- list( c("AAG", "GAG"), c("AAG", "BG"),c("AAG", "HAT") )
genes_to_highlight <- c('PHO85', 'MUD1', 'SEC14', 'SPT14', 'MAF1', 'UBC12', 'BET4', 'YCL002C', 'NYV1', 'RAD14')
alternative_3SS_only$threeSS_type_motif_factor
p2 <- ggboxplot(alternative_3SS_only %>% filter(threeSS_type_motif_factor != 'other', threeSS_type_motif_factor != 'AC'), x = "threeSS_type_motif_factor", y = "log10_prp18_effect_SE_with_pseudo_counts",
                color = "threeSS_type_motif_factor",  width = BOX_WIDTH)+ 
  geom_hline(yintercept = 0, linetype = 'dotted') + 
  geom_jitter(width = 0.15, size = 1, alpha = .15) +
  geom_text_repel(aes(label=ifelse(common_gene_name %in% genes_to_highlight & total_gapped_reads > 5 & annotated_to_alternative_3SS_dist < 0 & upf1_prp18_mean_FAnS_with_pseudo_fraction >= 0.1, as.character(common_gene_name),'')) , nudge_x = .5, size = 5, force = 1, fontface = 'italic') #+
# stat_compare_means(comparisons = my_comparisons) # Add pairwise comparisons p-value stat_compare_means(label.y = 3.5)     # Add global p-value
p2 # , nudge_x = .2, size = 4 
ggsave(paste(FIGURE_DIR,'all_alternative_junctions_threeSS_type_by_alt_fraction_FAnS_prp18_effect.svg', sep = '') , width = 6, height = 6)

## START Branchpoint analysis
alternative_3SS_only %>% filter(common_gene_name %in% genes_to_highlight, annotated_best_BP_sequence != 'TACTAAC') %>% 
  select(common_gene_name, threeSS_seq, annotated_to_alternative_3SS_dist, coord_for_best_match_to_BP_consensus, best_BP_sequence, annotated_coord_for_best_match_to_BP_consensus, annotated_best_BP_sequence ) %>%
  print(n=100)
alternative_3SS_only %>% filter(coord_for_best_match_to_BP_consensus == annotated_coord_for_best_match_to_BP_consensus)
ggscatter(alternative_3SS_only %>% filter(threeSS_type != 'other',threeSS_type != 'AC'), y = "upf1_prp18_mean_FAnS_with_pseudo_fraction", x = "upf1_mean_FAnS_with_pseudo_fraction", 
          size = 2, alpha = .3, 
          facet.by = c("annotated_best_BP_sequence"), color = c("threeSS_type_factor") ) + 
  geom_text_repel(aes(label=ifelse(annotated_best_BP_sequence != 'TACTAAC' & total_gapped_reads > 5& upf1_prp18_mean_FAnS_with_pseudo_fraction >= 0.1, as.character(paste(common_gene_name, threeSS_seq, sep = ' ')),'')) , nudge_y = 1) +
  geom_abline(slope = 1, linetype = 'dotted') +
  scale_x_log10(limits = c(-Inf, 10) ) + scale_y_log10(limits = c(-Inf, 10)) 

annotated_only <- annotated_only %>% mutate(BP_category = if_else(substring(annotated_best_BP_sequence, 1, 7) == 'TACTAAC', annotated_best_BP_sequence,
                                                                  'non-ACTAAC') )
my_comparisons <- list( c('TACTAAC', 'non-ACTAAC') )

ggboxplot(annotated_only, x = "canonical_BP", y = "log10_prp18_effect_SE_with_pseudo_counts",
          color = "RPG", width = .4, facet.by = c("RPG"), nrow = 1) +
  geom_jitter(width = 0.1, size = 1, alpha = 0.2)  + 
  geom_abline(intercept = 0, slope = 0, linetype = 'dashed', alpha = .3)  + 
  geom_text_repel(aes(label=ifelse( log10_prp18_effect_SE_with_pseudo_counts > 0 | log10_prp18_effect_SE_with_pseudo_counts < -.4,    as.character(paste(common_gene_name, threeSS_seq, sep = ' ')),'')), nudge_x = .5 ) +
  # geom_text(label = filtered_tbl$common_gene_name, nudge_x  = NUDGE_X, check_overlap = TRUE, size = 5) +
  stat_compare_means(comparisons = my_comparisons) # Add pairwise comparisons p-value

ggplot(annotated_only, aes( x = RPG, y = log10_prp18_effect_SE_with_pseudo_counts) ) + geom_boxplot()  +
  geom_jitter(width = 0.1, size = 1, alpha = 0.2)  + facet_grid(~ BP_category) + ggtitle('Canonical Branchpoint TACTAAC') +
  stat_compare_means(comparisons = my_comparisons)

alternative_3SS_only <- alternative_3SS_only %>% mutate(canonical_BP = best_BP_sequence == 'TACTAAC')

my_comparisons <- list( c("CAG", "TAG"), c("CAG", "AAG"),c("CAG", "GAG"), c("TAG", "AAG") , c("TAG", "GAG"), c("GAG", "AAG") )
my_comparisons <- list( c(TRUE, FALSE) )
alternative_3SS_only$threeSS_type
ggboxplot(alternative_3SS_only, x = "canonical_BP", y = "log10_prp18_effect_FAnS_with_pseudo_counts",
          color = "threeSS_type_factor", nrow = 1) +
  geom_jitter(width = 0.2, size = 1, alpha = 0.1) + 
  geom_abline(intercept = 0, slope = 0, linetype = 'dashed', alpha = .5)  + 
  facet_grid(~ threeSS_type_factor) +
  geom_text_repel(aes(label=ifelse(total_gapped_reads > 5 & upf1_prp18_mean_FAnS_with_pseudo_fraction >= 0.3, as.character(paste(common_gene_name, threeSS_seq, BP_3SS_dist, sep = ' ')),'')) , nudge_x = .5, size = 3, force = 1) +
  stat_compare_means(comparisons = my_comparisons) # Add pairwise comparisons p-value
## END Branchpoint analysis

my_comparisons <- list( c("CAG", "TAG"), c("CAG", "AAG"),c("CAG", "GAG"), c("TAG", "AAG") , c("TAG", "GAG"), c("GAG", "AAG") )
ggboxplot(alternative_3SS_only %>% group_by(threeSS_type), x = "threeSS_type_factor", y = "log10_prp18_effect_FAnS_with_pseudo_counts",
          color = "threeSS_type", nrow = 1) +
  geom_jitter(width = 0.2, size = 1, alpha = 0.1) + 
  geom_abline(intercept = 0, slope = 0, linetype = 'dashed', alpha = .5)  + 
  geom_text_repel(aes(label=ifelse(total_gapped_reads > 5 & upf1_prp18_mean_FAnS_with_pseudo_fraction >= 0.1, as.character(paste(common_gene_name, threeSS_seq, BP_3SS_dist, sep = ' ')),'')) , nudge_x = .5, size = 3, force = 1) +
  stat_compare_means(comparisons = my_comparisons) # Add pairwise comparisons p-value
ggsave(paste(FIGURE_DIR,'all_alternative_junctions_threeSS_type_by_alt_fraction_FAnS_prp18_effect_labeled.svg', sep = '') , width = 6, height = 6)

## BEGIN Fig S7
p1 <- ggscatter(alternative_only %>% filter(total_gapped_reads > 5), y = "upf1_prp18_mean_SE_with_pseudo_fraction", x = "upf1_mean_SE_with_pseudo_fraction", 
                color = "threeSS_seq", fill ="threeSS_type_factor", size = 1, alpha = .2, stroke = .5,
                facet.by = "threeSS_seq", axes = FALSE) +
  geom_abline(slope = 1, linetype = 'dotted') +
  # geom_text_repel(aes(label=ifelse( (threeSS_type %in% c("BG", "non-G") & upf1_prp18_mean_FAnS_with_pseudo_fraction >= 0.01) | ( upf1_prp18_mean_FAnS_with_pseudo_fraction > .9 & !(common_gene_name  %in% c('GCR1')) &  !(threeSS_type  %in% c('YAG', 'RAG') ) ),    as.character(paste(common_gene_name, threeSS_seq, sep = ' ')),'')), nudge_y = 1 ) +
  # geom_text_repel(aes(label=ifelse(upf1_prp18_mean_FAnS_with_pseudo_fraction >= 1,    as.character(paste(common_gene_name, threeSS_seq, sep = ' ')),'')), nudge_y = 1 ) +
  
  scale_x_log10(limits = c(-Inf, 10) ) + scale_y_log10(limits = c(-Inf, 10))  

p1
ggsave(paste(FIGURE_DIR, 'upf1_prp18_vs_upf1_FAnS_by_3SS_type.svg'), width = 12, height = 12)

p2 <- ggscatter(alternative_only, y = "upf1_prp18_mean_SE_with_pseudo_fraction", x = "upf1_mean_SE_with_pseudo_fraction", 
                color = "threeSS_type_factor",fill ="threeSS_type_factor", size = 1, alpha = .2, stroke = .5,
                facet.by = "threeSS_type_factor", axes = FALSE, nrow = 2) + 
  geom_abline(slope = 1, linetype = 'dotted') +
  scale_x_log10(limits = c(-Inf, 10) ) + scale_y_log10(limits = c(-Inf, 10)) 
p2
plot_grid(p2, p1)
ggsave(paste(FIGURE_DIR, 'upf1_prp18_vs_upf1_AS_SE_by_3SS_type.svg'), width = 9, height = 6)
## END Fig S7

# BEGIN ACCESSIBILITY FIGURES additional figures for Supp
# %>% filter(alt_5SS_alt_3SS == 'FALSE TRUE')
my_comparisons <- list( c("CAG", "TAG"), c("CAG", "AAG"),c("CAG", "GAG"), c("TAG", "AAG") , c("TAG", "GAG"), c("GAG", "AAG") )

## these plots below suggest overall there is generally equal accessibility for sites that are used, regardless of 3'SS motif
ggboxplot(filtered_tbl, x = "threeSS_type_factor", y = "accessibility",
          color = "threeSS_type_factor", facet.by = c("RPG"), nrow = 1) +
  geom_jitter(width = 0.2, size = 1, alpha = 0.1) + 
  geom_abline(intercept = 0, slope = 0, linetype = 'dashed', alpha = .5)  + 
  # geom_text(label = filtered_tbl$common_gene_name, nudge_x  = NUDGE_X, check_overlap = TRUE, size = 5) +
  stat_compare_means(comparisons = my_comparisons) # Add pairwise comparisons p-value

all_trint %>% filter(seq == 'GAG', alt_BP_3SS_dist > 10  ) %>% select(annotated_junction, linear_distance, alt_BP_3SS_dist, annotated_to_alternative_3SS_dist, accessibility, annotated_start, annotated_end) %>% arrange(annotated_to_alternative_3SS_dist) %>% print(n=100)
all_trint$annotated_to_alternative_3SS_dist
all_trint <- all_trint %>% filter( alt_BP_3SS_dist > 15)
all_trint_medians <- all_trint %>% group_by(seq) %>% summarize(the_medians = median(as.numeric(accessibility), na.rm = TRUE)) %>% arrange(the_medians) 
all_trint_medians %>% print(n=200)
all_trint_medians_factor <- factor(all_trint_medians$seq, levels = all_trint_medians$seq)
all_trint <- all_trint %>% mutate( seq_factor = factor(seq, levels = all_trint_medians_factor)  )


all_trint_randomized %>% print(n=200)
all_trint_randomized <- all_trint_randomized %>% filter( alt_BP_3SS_dist > 15)
all_trint_randomized %>% filter(seq == 'GAG', alt_BP_3SS_dist > 10  ) %>% select(annotated_junction, linear_distance, alt_BP_3SS_dist, annotated_to_alternative_3SS_dist, accessibility, annotated_start, annotated_end) %>% arrange(annotated_to_alternative_3SS_dist) %>% print(n=100)
all_trint_medians <- all_trint_randomized %>% group_by(seq) %>% summarize(the_medians = median(as.numeric(accessibility), na.rm = TRUE)) %>% arrange(the_medians) 
all_trint_medians_factor <- factor(all_trint_medians$seq, levels = all_trint_medians$seq)
all_trint_randomized <- all_trint_randomized %>% mutate( seq_factor = factor(seq, levels = all_trint_medians_factor)  )

my_comparisons <- list( c("AAG", "GAG"), c("AAG", "BG"),c("AAG", "HAT"), c("TAG", "AAG") , c("CAG", "AAG") )
ggboxplot(all_trint_randomized, x = "seq_factor", y = "accessibility",
          color = "threeSS_type_factor", nrow = 1) +
  geom_jitter(width = 0.2, size = 1, alpha = 0.1) + 
  geom_abline(intercept = 0, slope = 0, linetype = 'dashed', alpha = .5)  + 
  # geom_text(label = filtered_tbl$common_gene_name, nudge_x  = NUDGE_X, check_overlap = TRUE, size = 5) +
  stat_compare_means(comparisons = my_comparisons) + # Add pairwise comparisons p-value
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(filename = paste(FIGURE_DIR, 'accessibility_by_threeSS_type.svg'), device = 'svg', width = 6, height = 6) 

my_comparisons <- list( c("AAG", "GAG"), c("AAG", "BG"),c("AAG", "HAT"), c("TAG", "AAG") , c("CAG", "AAG") )
ggboxplot(all_trint, x = "seq_factor", y = "accessibility",
          color = "threeSS_type_factor", nrow = 1) +
  geom_jitter(width = 0.2, size = 1, alpha = 0.1) + 
  geom_abline(intercept = 0, slope = 0, linetype = 'dashed', alpha = .5)  + 
  # geom_text(label = filtered_tbl$common_gene_name, nudge_x  = NUDGE_X, check_overlap = TRUE, size = 5) +
  stat_compare_means(comparisons = my_comparisons) + # Add pairwise comparisons p-value
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
# END ACCESSIBILITY FIGURES additional figures for Supp

# BEGIN Fig S6
my_comparisons <- list( c(TRUE, FALSE) )
ggboxplot(filtered_tbl, x = "RPG", y = "upf1_prp18_mean_reads",
          color = "RPG", width = .4, facet.by = c("alt_5SS_alt_3SS"), nrow = 1) +
  geom_jitter(width = 0.1, size = 1, alpha = 0.1) + 
  geom_abline(intercept = 0, slope = 0, linetype = 'dashed', alpha = .5)  + scale_y_log10() +
  # geom_text(label = filtered_tbl$common_gene_name, nudge_x  = NUDGE_X, check_overlap = TRUE, size = 5) +
  stat_compare_means(comparisons = my_comparisons) # Add pairwise comparisons p-value
ggsave(filename = paste(DIR, 'upf1_prp18_mean_reads_faceted_RPG_vs_nonRPG_ann.svg'), device = 'svg', width = 12, height = 6) 

my_comparisons <- list( c(TRUE, FALSE) )
ggboxplot(filtered_tbl, x = "RPG", y = "upf1_mean_reads",
          color = "RPG", width = .4, facet.by = c("alt_5SS_alt_3SS"), nrow = 1) +
  geom_jitter(width = 0.1, size = 1, alpha = 0.1) + 
  geom_abline(intercept = 0, slope = 0, linetype = 'dashed', alpha = .5)  + scale_y_log10() +
  # geom_text(label = filtered_tbl$common_gene_name, nudge_x  = NUDGE_X, check_overlap = TRUE, size = 5) +
  stat_compare_means(comparisons = my_comparisons) # Add pairwise comparisons p-value
ggsave(filename = paste(DIR, 'upf1_mean_reads_faceted_RPG_vs_nonRPG_ann.svg'), device = 'svg', width = 12, height = 6) 

my_comparisons <- list( c(TRUE, FALSE) )
ggboxplot(filtered_tbl, x = "RPG", y = "upf1_prp18_mean_intron_reads",
          color = "RPG", width = .4, facet.by = c("alt_5SS_alt_3SS"), nrow = 1) +
  geom_jitter(width = 0.1, size = 1, alpha = 0.1) + 
  geom_abline(intercept = 0, slope = 0, linetype = 'dashed', alpha = .5)  + scale_y_log10() +
  # geom_text(label = filtered_tbl$common_gene_name, nudge_x  = NUDGE_X, check_overlap = TRUE, size = 5) +
  stat_compare_means(comparisons = my_comparisons) # Add pairwise comparisons p-value
ggsave(filename = paste(DIR, 'upf1_prp18_mean_intron_reads_faceted_RPG_vs_nonRPG_ann.svg'), device = 'svg', width = 12, height = 6) 

my_comparisons <- list( c(TRUE, FALSE) )
ggboxplot(filtered_tbl, x = "RPG", y = "upf1_mean_intron_reads",
          color = "RPG", width = .4, facet.by = c("alt_5SS_alt_3SS"), nrow = 1) +
  geom_jitter(width = 0.1, size = 1, alpha = 0.1) + 
  geom_abline(intercept = 0, slope = 0, linetype = 'dashed', alpha = .5)  + scale_y_log10() +
  # geom_text(label = filtered_tbl$common_gene_name, nudge_x  = NUDGE_X, check_overlap = TRUE, size = 5) +
  stat_compare_means(comparisons = my_comparisons) # Add pairwise comparisons p-value
ggsave(filename = paste(DIR, 'upf1_mean_intron_reads_faceted_RPG_vs_nonRPG_ann.svg'), device = 'svg', width = 12, height = 6) 

ggboxplot(filtered_tbl %>% filter(annotated == TRUE), x = "RPG", y = "log10_prp18_effect_SE_with_pseudo_counts",
          color = "RPG", width = .4, facet.by = c("threeSS_type"), nrow = 1) +
  geom_jitter(width = 0.1, size = 1, alpha = 0.2)  + 
  geom_abline(intercept = 0, slope = 0, linetype = 'dashed', alpha = .3)  + 
 # geom_text_repel(aes(label=ifelse( log10_prp18_effect_SE_with_pseudo_counts > 0 | log10_prp18_effect_SE_with_pseudo_counts < -.5,    as.character(paste(common_gene_name, threeSS_seq, sep = ' ')),'')), nudge_x = .5 ) +
  # geom_text(label = filtered_tbl$common_gene_name, nudge_x  = NUDGE_X, check_overlap = TRUE, size = 5) +
  stat_compare_means(comparisons = my_comparisons) # Add pairwise comparisons p-value
ggsave(filename = paste(FIGURE_DIR, 'prp18_effect_SE_faceted_RPGs.svg'), device = 'svg', width = 3, height = 6) 
# END Fig S6

my_comparisons <- list( c(TRUE, FALSE) )
ggboxplot(alternative_3SS_only %>% filter(threeSS_type_factor != 'other', threeSS_type_factor != 'AC' ), x = "RPG", y = "alternative_U_score",
          color = "threeSS_type_factor", facet.by =  "threeSS_type_factor", nrow = 1) +
  geom_jitter(width = 0.2, size = 1, alpha = 0.1) + 
  geom_abline(intercept = 0, slope = 0, linetype = 'dashed', alpha = .5)  + 
  # geom_text(label = filtered_tbl$common_gene_name, nudge_x  = NUDGE_X, check_overlap = TRUE, size = 5) +
  stat_compare_means(comparisons = my_comparisons) # Add pairwise comparisons p-value + 
ggsave(filename = paste(FIGURE_DIR, 'alternative_U_score_by_threeSS_type_upstream_polyU.svg'), device = 'svg', width = 8, height = 4) 

my_comparisons <- list(  c("AAG", "BG"),c("AAG", "HAT"))
ggboxplot(annotated_only, x = "RPG", y = "alternative_U_score", color = "RPG", 
          nrow = 1) +
  geom_jitter(width = 0.2, size = 1, alpha = 0.1) + 
  geom_abline(intercept = 0, slope = 0, linetype = 'dashed', alpha = .5)  + 
  # geom_text(label = filtered_tbl$common_gene_name, nudge_x  = NUDGE_X, check_overlap = TRUE, size = 5) +
  stat_compare_means(comparisons = my_comparisons) # Add pairwise comparisons p-value
ggsave(filename = paste(FIGURE_DIR, 'alternative_U_score_by_threeSS_type_upstream_polyU.svg'), device = 'svg', width = 2, height = 4) 

alternative_3SS_only$threeSS_type_motif_factor
ggboxplot(alternative_3SS_only, x = "threeSS_type_motif_factor", y = "alt_minus_ann_U_score",
          color = "RPG", nrow = 1, facet.by = c( "RPG")) +
  geom_jitter(width = 0.2, size = 1, alpha = 0.1) + 
  geom_abline(intercept = 0, slope = 0, linetype = 'dashed', alpha = .5)  + 
  # geom_text(label = alternative_3SS_only$common_gene_name, nudge_x  = NUDGE_X, check_overlap = TRUE, size = 5) +
  stat_compare_means(comparisons = my_comparisons) # Add pairwise comparisons p-value
ggsave(filename = paste(FIGURE_DIR, 'alt_minus_ann_U_score_by_threeSS_type_by_RPG.svg'), device = 'svg', width = 8, height = 4) 


## need to revise this code below for # of alt. junctions by ann. U score
df <- combined_alt_tally %>% left_join(annotated_only %>% select(annotated_junction, annotated_U_score), by = c("annotated_junction"))
df <- df %>% mutate( alt_j = ifelse(`total number of alternative 3SS junctions min_one_thousandth_FAnS` > 2, 3, `total number of alternative 3SS junctions min_one_thousandth_FAnS` ))
p1 <- ggplot( df, aes(x = alt_j, y = annotated_U_score )) + geom_boxplot(aes(colour = as.factor(alt_j) )) + geom_jitter()
p2 <- ggplot( df, aes(x = alt_j )) + geom_density(aes(colour = as.factor(alt_j) )) 

df %>% group_by(RPG) %>% tally(`total number of alternative 3SS junctions min_one_thousandth_FAnS` )
df %>% group_by(`total number of alternative 5SS junctions min_one_thousandth_FAnS`) %>%  count()
## For some reason, getting NA's on some of these
df$alt_j <- as.factor(df$alt_j)

p1
filename <- (paste(FIGURE_DIR,'impact_of_annotated_U_score_on_SE.svg', sep = '') )
ggsave(filename, width = 6, height = 6)
p2
filename <- (paste(FIGURE_DIR,'impact_of_annotated_U_score_on_alt_junctions_identified.svg', sep = '') )
ggsave(filename, width = 8, height = 6)

my_comparisons <- list( c(2, 1), c(2, 0) )
ggboxplot(df, x = "alt_j", y = "annotated_U_score",
          color = "alt_j", nrow = 1) +
  geom_jitter(width = 0.2, size = 1, alpha = 0.1) + 
  geom_abline(intercept = 0, slope = 0, linetype = 'dashed', alpha = .5)  + 
  # geom_text(label = alternative_3SS_only$common_gene_name, nudge_x  = NUDGE_X, check_overlap = TRUE, size = 5) +
  stat_compare_means(comparisons = my_comparisons) # Add pairwise comparisons p-value
ggsave(filename = paste(FIGURE_DIR, 'alt_minus_ann_U_score_by_threeSS_type_by_RPG.svg'), device = 'svg', width = 8, height = 4) 


my_comparisons <- list( c(TRUE, FALSE) )
ggboxplot(filtered_tbl, x = "threeSS_type", y = "alternative_U_score",
           color = "threeSS_type", width = .4, facet.by = c("alt_5SS_alt_3SS"), nrow = 1) +
  geom_jitter(width = 0.1, size = 1, alpha = 0.1) + 
  geom_abline(intercept = 0, slope = 0, linetype = 'dashed', alpha = .5)  + scale_y_log10() +
  # geom_text(label = filtered_tbl$common_gene_name, nudge_x  = NUDGE_X, check_overlap = TRUE, size = 5) +
  # Add pairwise comparisons p-value
  stat_compare_means(comparisons = my_comparisons)  + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(filename = paste(DIR, 'upf1_prp18_mean_reads_faceted_RPG_vs_nonRPG_ann.svg'), device = 'svg', width = 12, height = 6) 

filtered_tbl %>% filter(upf1_prp18_mean_FAnS > 1) %>% select(common_gene_name, threeSS_seq, annotated_to_alternative_3SS_dist, BP_3SS_dist,  total_gapped_reads)
# %>% group_by(common_gene_name) %>% count()
filtered_tbl %>% group_by(alt_5SS_alt_3SS) %>%  filter(log10_prp18_effect_SE_with_pseudo_counts > 0) %>% count()
# mutate( site_used = ifelse(upf1_prp18_FAnS > 0.0001, TRUE, FALSE ) 
min(filtered_tbl$BP_3SS_dist)
filtered_tbl %>% filter(BP_3SS_dist == 7)
alternative_3SS_only %>% group_by(fiveSS_type) %>% count()
alternative_5SS_only %>% group_by(fiveSS_type) %>% count()
annotated_only %>% group_by(fiveSS_type) %>% count()
my_comparisons <- list( c(FALSE, TRUE) )
all_trint$threeSS_type_factor
## %>% filter(alt_BP_3SS_dist > 20)

dev.off()
all_trint$RPG
## ## BEGIN accessibility of motifs
my_comparisons <- list( c(TRUE, FALSE) )
ggboxplot(all_trint %>% filter(alt_BP_3SS_dist > 7, annotated_to_alternative_3SS_dist < 100, threeSS_type_factor != 'AC', threeSS_type_factor != 'other' ), x = "site_used", y = "accessibility", color = "threeSS_type_annotated_factor",
           width = .8, facet.by = c("threeSS_type_annotated_factor"), nrow = 1 ) +
  geom_jitter(width = 0.2, size = 1, alpha = 0.1) + 
  geom_abline(intercept = 0, slope = 0, linetype = 'dashed', alpha = .5)  + 
  # geom_text(label = filtered_tbl$common_gene_name, nudge_x  = NUDGE_X, check_overlap = TRUE, size = 5) +
  stat_compare_means(comparisons = my_comparisons) + # Add pairwise comparisons p-value
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(filename = paste(FIGURE_DIR, 'accessibility_by_site_used_with_50nt_DS.svg'), device = 'svg', width = 6, height = 6) 

##U-score of motifs
my_comparisons <- list( c(TRUE, FALSE) )
ggboxplot(all_trint %>% filter(alt_BP_3SS_dist > 7, annotated_to_alternative_3SS_dist < 100, threeSS_type_factor != 'AC', threeSS_type_factor != 'other' ), x = "site_used", y = "U_score", color = "threeSS_type_annotated_factor",
          width = .8, facet.by = c("threeSS_type_annotated_factor"), nrow = 1) +
  geom_jitter(width = 0.2, size = 1, alpha = 0.1) + 
  geom_abline(intercept = 0, slope = 0, linetype = 'dashed', alpha = .5)  + 
  # geom_text(label = filtered_tbl$common_gene_name, nudge_x  = NUDGE_X, check_overlap = TRUE, size = 5) +
  stat_compare_means(comparisons = my_comparisons) + # Add pairwise comparisons p-value
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(filename = paste(FIGURE_DIR, 'U_score_by_site_used_with_50nt_DS.svg'), device = 'svg', width = 6, height = 6) 

## BEGIN accessibility of motifs
my_comparisons <- list( c("AAG", "GAG"), c("GAG", "BG"),c("GAG", "HAT"), c("TAG", "GAG") , c("CAG", "GAG") )
ggboxplot(all_trint %>% filter(alt_BP_3SS_dist > 6, annotated_to_alternative_3SS_dist < 5), x = "threeSS_type_factor", y = "accessibility", color = "threeSS_type_factor",
          width = .4, nrow = 1) +
  geom_jitter(width = 0.1, size = 1, alpha = 0.1) + 
  geom_abline(intercept = 0, slope = 0, linetype = 'dashed', alpha = .5)  + 
  # geom_text(label = filtered_tbl$common_gene_name, nudge_x  = NUDGE_X, check_overlap = TRUE, size = 5) +
  stat_compare_means(comparisons = my_comparisons) + # Add pairwise comparisons p-value
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(filename = paste(DIR, 'accessibility_by_site_used.svg'), device = 'svg', width = 6, height = 6) 

GAG <- all_trint %>% filter(alt_BP_3SS_dist > 10, annotated_to_alternative_3SS_dist < 5, motif == 'GAG') %>% arrange(annotated_to_alternative_3SS_dist ) %>% print(n=100)
ggdensity(all_trint %>% filter(alt_BP_3SS_dist > 10, annotated_to_alternative_3SS_dist < 5, annotated_to_alternative_3SS_dist != 0), x = "annotated_to_alternative_3SS_dist", facet.by = c("threeSS_type_factor") ) + 
  xlim(-100, 5)

ggdensity(all_trint %>% filter(alt_BP_3SS_dist > 10, annotated_to_alternative_3SS_dist != 0), x = "annotated_to_alternative_3SS_dist", facet.by = c("threeSS_type_factor") ) + 
  xlim(-50, 50)

## load all trinucleotide file with 50 bp downstream for figures below
ggscatter(all_trint %>% filter(alt_BP_3SS_dist > 10), x = "annotated_to_alternative_3SS_dist", y = "U_score", facet.by = c("threeSS_type_factor"), alpha = .2, color = "threeSS_type_factor" ) + 
  xlim(-50, 50)

ggscatter(all_trint %>% filter(alt_BP_3SS_dist > 10), x = "annotated_to_alternative_3SS_dist", y = "accessibility", facet.by = c("threeSS_type_factor"), alpha = .2, color = "threeSS_type_factor" ) + 
  xlim(-50, 50)

my_comparisons <- list( c("AAG", "GAG"), c("AAG", "BG"),c("AAG", "HAT"), c("TAG", "AAG") , c("CAG", "AAG") )
ggscatter(all_trint %>% filter(alt_BP_3SS_dist > 10), x = "alt_BP_3SS_dist", y = "accessibility", color = "threeSS_type_factor", facet.by = c("threeSS_type_factor")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + scale_x_log10()
ggsave(filename = paste(DIR, 'accessibility_by_site_used.svg'), device = 'svg', width = 6, height = 6) 

my_comparisons <- list( c("AAG", "GAG"), c("AAG", "BG"),c("AAG", "HAT"), c("TAG", "AAG") , c("CAG", "AAG") )
ggboxplot(filtered_tbl, x = "threeSS_type_factor", y = "accessibility", color = "threeSS_type_factor",
          width = .4, nrow = 1) +
  geom_jitter(width = 0.1, size = 1, alpha = 0.1) + 
  geom_abline(intercept = 0, slope = 0, linetype = 'dashed', alpha = .5)  + 
  # geom_text(label = filtered_tbl$common_gene_name, nudge_x  = NUDGE_X, check_overlap = TRUE, size = 5) +
  stat_compare_means(comparisons = my_comparisons) + # Add pairwise comparisons p-value
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(filename = paste(DIR, 'accessibility_by_site_used.svg'), device = 'svg', width = 6, height = 6) 
## END accessibility of motifs

## BEGIN exonic nt analysis
my_comparisons <- list( c("G", "A"), c("C", "A"),c("T", "A") )
p2 <- ggboxplot(alternative_3SS_only, x = "DS_3SS_1", y = "log10_upf1_prp18_mean_SE",
                color = "threeSS_type_factor", facet.by = c("threeSS_type_factor"), width = BOX_WIDTH)+ 
  geom_hline(yintercept = 0, linetype = 'dotted') + 
  geom_jitter(width = 0.15, size = 1, alpha = .15) +
  geom_text_repel(aes(label=ifelse(common_gene_name %in% genes_to_highlight & total_gapped_reads > 5 & annotated_to_alternative_3SS_dist < 0 & upf1_prp18_mean_FAnS_with_pseudo_fraction >= 0.1, as.character(paste(common_gene_name, threeSS_seq, sep = ' ')),'')) , nudge_x = .5, size = 3, force = 1) +
  stat_compare_means(comparisons = my_comparisons) # Add pairwise comparisons p-value stat_compare_means(label.y = 3.5)     # Add global p-value
p2 

ggboxplot(filtered_tbl, x = "US_5SS_1", y = "log10_upf1_mean_SE",
          color = "alt_5SS_alt_3SS", width = .4, facet.by = c("alt_5SS_alt_3SS"), nrow = 1) +
  geom_jitter(width = 0.1, size = 1, alpha = 0.1) + 
  geom_abline(intercept = 0, slope = 0, linetype = 'dashed', alpha = .5)  # + scale_y_log10() +
  # geom_text(label = filtered_tbl$common_gene_name, nudge_x  = NUDGE_X, check_overlap = TRUE, size = 5) +
  stat_compare_means(comparisons = my_comparisons) # Add pairwise comparisons p-value
## END exonic nt analysis


## Fig S7A
ggscatter(filtered_tbl_min_6_reads, y = "upf1_prp18_mean_SE_with_pseudo_fraction", x = "upf1_mean_SE_with_pseudo_fraction", 
          size = 2, alpha = .3, 
          facet.by = c("fiveSS_type_factor", "threeSS_type_factor"), color = c("alt_5SS_alt_3SS") ) + 
  geom_abline(slope = 1, linetype = 'dotted') +
#  geom_text_repel(aes(label=ifelse(common_gene_name %in% genes_to_highlight & total_gapped_reads > 5 & upf1_prp18_mean_FAnS_with_pseudo_fraction >= 0.1, as.character(paste(common_gene_name, threeSS_seq, sep = ' ')),'')) , nudge_x = .5, size = 3, force = 1) +
  scale_x_log10(limits = c(-Inf, 1) ) + scale_y_log10(limits = c(-Inf, 1))
filename <- (paste(FIGURE_DIR,'all_introns_5SS_by_3SS_type_SE.svg', sep = '') )
ggsave(filename, width = 14, height = 6)


## BEGIN FIG 3D on RNA structure
marker = list(color = brewer.pal(7, "Dark2"))
marker[[1]][2]

filtered_tbl_min_6_reads_with_structure_annotated <- filtered_tbl_min_6_reads_with_structure %>%
  filter(alt_5SS_alt_3SS == 'FALSE FALSE')

max_value_index <-which.max(density(filtered_tbl_min_6_reads_with_structure_annotated$effective_BP_3SS_dist)$y)
density(filtered_tbl_min_6_reads_with_structure_annotated$effective_BP_3SS_dist)$x[max_value_index]

filtered_tbl_min_6_reads_with_structure_alt_3SS <- filtered_tbl_min_6_reads_with_structure %>%
  filter(alt_5SS_alt_3SS == 'FALSE TRUE')

max_value_index <-which.max(density(filtered_tbl_min_6_reads_with_structure_alt_3SS$effective_BP_3SS_dist)$y)
density(filtered_tbl_min_6_reads_with_structure_alt_3SS$effective_BP_3SS_dist)$x[max_value_index]

filtered_tbl_min_6_reads_without_structure_alt_3SS <- filtered_tbl_min_6_reads_with_no_structure %>%
  filter(alt_5SS_alt_3SS == 'FALSE TRUE')

max_value_index <-which.max(density(filtered_tbl_min_6_reads_without_structure_alt_3SS$effective_BP_3SS_dist)$y)
density(filtered_tbl_min_6_reads_without_structure_alt_3SS$effective_BP_3SS_dist)$x[max_value_index]

filtered_tbl_min_6_reads_without_structure_alt_3SS <- filtered_tbl_min_6_reads_with_no_structure %>%
  filter(alt_5SS_alt_3SS == 'TRUE TRUE')
max_value_index <-which.max(density(filtered_tbl_min_6_reads_without_structure_alt_3SS$effective_BP_3SS_dist)$y)
density(filtered_tbl_min_6_reads_without_structure_alt_3SS$effective_BP_3SS_dist)$x[max_value_index]

filtered_tbl_min_6_reads_with_structure_alt_3SS <- filtered_tbl_min_6_reads_with_structure %>%
  filter(alt_5SS_alt_3SS == 'TRUE TRUE')
max_value_index <-which.max(density(filtered_tbl_min_6_reads_with_structure_alt_3SS$effective_BP_3SS_dist)$y)
density(filtered_tbl_min_6_reads_with_structure_alt_3SS$effective_BP_3SS_dist)$x[max_value_index]

filtered_tbl_min_6_reads_with_structure$threeSS_type
ggplot() + geom_density( data = filtered_tbl_min_6_reads_with_structure %>% filter(threeSS_type != 'other',threeSS_type != 'AC', threeSS_type_annotated_factor != 'alternative_5SS_YAG'), aes(BP_3SS_dist), linetype = 'solid', color = marker[[1]][1]) + 
  geom_density(data = filtered_tbl_min_6_reads_with_structure %>% filter(threeSS_type != 'other',threeSS_type != 'AC', threeSS_type_annotated_factor != 'alternative_5SS_YAG'), aes(effective_BP_3SS_dist), linetype = 'dashed', color =marker[[1]][2])  +  
  geom_density(data = filtered_tbl_min_6_reads_with_no_structure %>% filter(threeSS_type != 'other',threeSS_type != 'AC', threeSS_type_annotated_factor != 'alternative_5SS_YAG'), aes(BP_3SS_dist), linetype = 'twodash', color = marker[[1]][3]) + 
  scale_x_continuous(limits = c(0,75), breaks = seq(  0, 1000, 10) ) +
  facet_grid(threeSS_type_annotated_factor ~ .)
filename <- (paste(FIGURE_DIR,'density_plots_effective_and_linear_BP_dist.eps', sep = '') )
ggsave(filename, width = 4, height = 8)
## END FIG 3D on RNA structure


df <- left_join(alt_junctions_per_intron_min_one_thousandth_FAnS, alternative_3SS_only, by = "annotated_junction")
# %>% filter(RPG == TRUE)
df1 <- alternative_3SS_only #  %>% filter(log10_upf1_mean_SE > -6  & log10_upf1_prp18_mean_SE > -6 )

fit <- lm(log10_upf1_mean_SE ~  annotated_U_score, data = df1 %>% filter(log10_upf1_mean_SE > -6) )
p1 <- ggplotRegression(fit, 'alternative_3SS_only')
p1

ggplot(alternative_3SS_only %>% filter(log10_upf1_mean_SE > -6), aes(annotated_U_score, upf1_mean_SE,  colour = RPG) ) +
  geom_point(alpha = .5) + scale_y_log10() + geom_smooth(method='lm', size = 1, alpha = 0.1)
filename <- (paste(FIGURE_DIR,'log10_upf1_mean_SE_by_annotated_U_score.svg', sep = '') )
ggsave(filename, width = 8, height = 6)

fit <- lm(log10_upf1_mean_SE ~  annotated_U_score, data = df1 %>% filter(RPG, log10_upf1_mean_SE > -6) )
p1 <- ggplotRegression(fit, 'alternative_3SS_only')
p1

fit <- lm(log10_upf1_mean_SE ~  annotated_U_score, data = df1 %>% filter(!RPG, log10_upf1_mean_SE > -6) )
p1 <- ggplotRegression(fit, 'alternative_3SS_only')
p1


fit <- lm(log10_upf1_prp18_mean_SE ~  annotated_U_score, data = df1 %>% filter(log10_upf1_prp18_mean_SE > -6))
p2 <- ggplotRegression(fit, 'alternative_3SS_only')
plot_grid( p1,p2,  ncol =1 )
filename <- (paste(FIGURE_DIR,'log10_upf1_mean_SE_by_annotated_U_score.eps', sep = '') )
ggsave(filename, width = 4, height = 8)

fit <- lm(log10_upf1_mean_SE ~  log2_BP_3SS_dist, data = df1 %>% filter(log10_upf1_mean_SE > -6) )
p1 <- ggplotRegression(fit, 'alternative_3SS_only')
fit <- lm(log10_upf1_prp18_mean_SE ~  log2_BP_3SS_dist, data = df1 %>% filter(log10_upf1_prp18_mean_SE > -6))
p2 <- ggplotRegression(fit, 'alternative_3SS_only')
plot_grid( p1,p2,  ncol =1 )
filename <- (paste(FIGURE_DIR,'log10_upf1_mean_SE_by_log2_BP_3SS_dist.eps', sep = '') )
ggsave(filename, width = 4, height = 8)

df1 <- alternative_3SS_only %>% filter(RPG == FALSE) %>% filter(log10_upf1_mean_SE > -6  & log10_upf1_prp18_mean_SE > -6 )
fit <- lm(log10_upf1_mean_SE ~  annotated_U_score, data = df1)
p1 <- ggplotRegression(fit, 'alternative_3SS_only')
fit <- lm(log10_upf1_prp18_mean_SE ~  annotated_U_score, data = df1)
p2 <- ggplotRegression(fit, 'alternative_3SS_only')
fit <- lm(log10_upf1_mean_SE ~ alternative_U_score, data = df1)
p3 <- ggplotRegression(fit, 'alternative_3SS_only')
fit <- lm(log10_upf1_prp18_mean_SE ~  alternative_U_score, data = df1)
p4 <- ggplotRegression(fit, 'alternative_3SS_only')
fit <- lm(log10_upf1_mean_SE ~ alt_minus_ann_U_score, data = df1)
p5 <- ggplotRegression(fit, 'alternative_3SS_only')
fit <- lm(log10_upf1_prp18_mean_SE ~  alt_minus_ann_U_score, data = df1)
p6 <- ggplotRegression(fit, 'alternative_3SS_only')
plot_grid( p1,p2,p3,p4,p5,p6,  ncol =2 )

fit <- lm(log10_upf1_mean_SE ~  annotated_poly_U_count, data = df1)
p1 <- ggplotRegression(fit, 'alternative_3SS_only')
fit <- lm(log10_upf1_prp18_mean_SE ~  annotated_poly_U_count, data = df1)
p2 <- ggplotRegression(fit, 'alternative_3SS_only')
plot_grid( p1,p2,  ncol =1 )

fit <- lm(log10_upf1_mean_SE ~  annotated_poly_Y_count, data = df1)
p1 <- ggplotRegression(fit, 'alternative_3SS_only')
fit <- lm(log10_upf1_prp18_mean_SE ~  annotated_poly_Y_count, data = df1)
p2 <- ggplotRegression(fit, 'alternative_3SS_only')
plot_grid( p1,p2,  ncol =1 )

fit <- lm(log10_upf1_mean_SE ~  log2_annotated_BP_3SS_dist, data = df1)
p1 <- ggplotRegression(fit, 'alternative_3SS_only')
fit <- lm(log10_upf1_prp18_mean_SE ~  log2_annotated_BP_3SS_dist, data = df1)
p2 <- ggplotRegression(fit, 'alternative_3SS_only')
plot_grid( p1,p2,  ncol =1 )


df1 <- alternative_3SS_only %>% filter(RPG == TRUE) %>% filter(log10_upf1_mean_SE > -6)

fit <- lm(log10_upf1_mean_SE ~  annotated_U_score, data = df1)
p1 <- ggplotRegression(fit, 'alternative_3SS_only')
fit <- lm(log10_upf1_prp18_mean_SE ~  annotated_U_score, data = df1)
p2 <- ggplotRegression(fit, 'alternative_3SS_only')
fit <- lm(log10_upf1_mean_SE ~ alternative_U_score, data = df1)
p3 <- ggplotRegression(fit, 'alternative_3SS_only')
fit <- lm(log10_upf1_prp18_mean_SE ~  alternative_U_score, data = df1)
p4 <- ggplotRegression(fit, 'alternative_3SS_only')
fit <- lm(log10_upf1_mean_SE ~ alt_minus_ann_U_score, data = df1)
p5 <- ggplotRegression(fit, 'alternative_3SS_only')
fit <- lm(log10_upf1_prp18_mean_SE ~  alt_minus_ann_U_score, data = df1)
p6 <- ggplotRegression(fit, 'alternative_3SS_only')
plot_grid( p1,p2,p3,p4,p5,p6,  ncol =2 )

df1 <- alternative_3SS_only %>% filter(log10_upf1_mean_SE > -6  & log10_upf1_prp18_mean_SE > -6 ) #%>% filter(RPG == FALSE)

fit <- lm(log10_upf1_mean_SE ~  log2_BP_3SS_dist, data = alternative_3SS_only %>% filter(log10_upf1_mean_SE > -6)  )
p1 <- ggplotRegression(fit, 'alt_3SS')
fit <- lm(log10_upf1_prp18_mean_SE ~  log2_BP_3SS_dist, data = alternative_3SS_only %>% filter(log10_upf1_prp18_mean_SE > -6) )
p2 <- ggplotRegression(fit, 'alt_3SS')
fit <- lm(log10_upf1_mean_SE ~  log2_effective_BP_3SS_dist, data = alternative_3SS_only %>% filter(log10_upf1_mean_SE > -6) )
p3 <- ggplotRegression(fit, 'alt_3SS')
fit <- lm(log10_upf1_prp18_mean_SE ~  log2_effective_BP_3SS_dist, data = alternative_3SS_only %>% filter(log10_upf1_prp18_mean_SE > -6) )
p4 <- ggplotRegression(fit, 'alt_3SS')
plot_grid( p1,p2,p3,p4,  ncol = 1 )
filename <- (paste(FIGURE_DIR,'log10_upf1_mean_SE_by_log2_BP_3SS_dist.eps', sep = '') )
ggsave(filename, width = 4, height = 14)


## Fig S
## linear gives better correlation than effective BP distance
fit <- lm(log10_prp18_effect_SE_with_pseudo_counts ~ log2_BP_3SS_dist, data = df1)
p4 <- ggplotRegression(fit, 'alternative_3SS_only')
fit <- lm(log10_prp18_effect_SE_with_pseudo_counts ~ log2_effective_BP_3SS_dist, data = df1)
p5 <- ggplotRegression(fit, 'alternative_3SS_only')

plot_grid(p4, p5, ncol =1 )
filename <- (paste(FIGURE_DIR,'BP_3SS_dist_and_U_score.svg', sep = '') )
ggsave(filename, width = 12, height = 16)

df1 <- alternative_3SS_only %>% filter(log10_upf1_mean_SE > -5  & log10_upf1_prp18_mean_SE > -5 )

## Fig S?
## SE and prp18 effect by BP 3SS distance all threeSS together
fit <- lm(log10_upf1_mean_SE ~ log2_BP_3SS_dist, data = df1)
p1 <- ggplotRegression(fit, 'log10_upf1_SE')
fit <- lm(log10_upf1_prp18_mean_SE ~ log2_BP_3SS_dist, data = df1)
p2 <- ggplotRegression(fit, 'log10_upf1_prp18_SE')
fit <- lm(log10_prp18_effect_SE_with_pseudo_counts ~ log2_BP_3SS_dist, data = df1)
p3 <- ggplotRegression(fit, 'log10_prp18_effect_SE')
fit <- lm(log10_prp18_effect_FAnS_with_pseudo_counts ~ log2_BP_3SS_dist, data = df1)
p4 <- ggplotRegression(fit, 'log10_prp18_effect_FAnS')
plot_grid(p1, p2, p3, p4)
filename <- (paste(FIGURE_DIR,'PRP18_effect_by_BP_dist_unused_sites_not_included.svg', sep = '') )
ggsave(filename, width = 12, height = 12)


## prp18 effect by BP 3SS distance by threeSS type

df <- alternative_3SS_only %>% filter(threeSS_type == 'GAG')
YAG_fit <- lm(log10_prp18_effect_SE_with_pseudo_counts ~ log2_BP_3SS_dist, data = df)
p1 <- ggplotRegression(YAG_fit, 'GAG')

df <- alternative_3SS_only %>% filter(threeSS_type == 'AAG')
RAG_fit <- lm(log10_prp18_effect_SE_with_pseudo_counts ~ log2_BP_3SS_dist, data = df)
p2 <- ggplotRegression(RAG_fit, 'AAG')

df <- alternative_3SS_only %>% filter(threeSS_type == 'CAG')
BG_fit <- lm(log10_prp18_effect_SE_with_pseudo_counts ~ log2_BP_3SS_dist, data = df)
p3 <- ggplotRegression(BG_fit, 'CAG')

df <- alternative_3SS_only %>% filter(threeSS_type == 'TAG')
nonG_fit <- lm(log10_prp18_effect_SE_with_pseudo_counts ~ log2_BP_3SS_dist, data = df)
p4 <- ggplotRegression(nonG_fit, 'TAG')

df <- alternative_3SS_only %>% filter(threeSS_type == 'BG')
nonG_fit <- lm(log10_prp18_effect_SE_with_pseudo_counts ~ log2_BP_3SS_dist, data = df)
p5 <- ggplotRegression(nonG_fit, 'BG')

df <- alternative_3SS_only %>% filter(threeSS_type == 'HAT')
nonG_fit <- lm(log10_prp18_effect_SE_with_pseudo_counts ~ log2_BP_3SS_dist, data = df)
p6 <- ggplotRegression(nonG_fit, 'HAT')

plot_grid(p1, p2, p3, p4, p5, p6)
ggsave(paste(FIGURE_DIR, 'ggplot_Regression_linear_BP_dist_by_3SS_type.svg', sep = ''), width = 8, height = 8)


## nt by nt prp18 effective by BP distance, drop-off at 12 and 47 nt
min(alternative_3SS_only$BP_3SS_dist  ) # alternative_3SS_only%>% filter(BP_3SS_dist < 55)
effective_fit <- lm(log10_prp18_effect_FAnS_with_pseudo_counts ~ log2_BP_3SS_dist, data = alternative_3SS_only )
p1 <- ggplotRegression(effective_fit, 'effective_fit')
p1
filename <- (paste(FIGURE_DIR,'PRP18_effect_by_BP_dist_scatter.svg', sep = '') )
ggsave(filename, width = 6, height = 4)


alternative_3SS_only$alt_minus_ann_U_score
fit <- lm(log10_upf1_prp18_mean_SE ~ alt_minus_ann_U_score, data = alternative_3SS_only %>% filter(log10_upf1_prp18_mean_SE > -6) )
p1 <- ggplotRegression(fit, 'alternative_3SS_only')
alternative_3SS_only$alt_minus_ann_U_score
fit <- lm(log10_upf1_mean_SE ~ alt_minus_ann_U_score, data = alternative_3SS_only %>% filter(log10_upf1_mean_SE > -6) )
p2 <- ggplotRegression(fit, 'alternative_3SS_only')
plot_grid(p1, p2)
filename <- (paste(FIGURE_DIR,'SE_mean_by_alt_to_ann_U_score.svg', sep = '') )
ggsave(filename, width = 6, height = 4)

alternative_3SS_only$alt_minus_ann_U_score
fit <- lm(log10_upf1_prp18_mean_SE ~ annotated_U_score, data = alternative_3SS_only %>% filter(log10_upf1_prp18_mean_SE > -6) )
p1 <- ggplotRegression(fit, 'alternative_3SS_only')
alternative_3SS_only$alt_minus_ann_U_score
fit <- lm(log10_upf1_mean_SE ~ annotated_U_score, data = alternative_3SS_only %>% filter(log10_upf1_mean_SE > -6) )
p2 <- ggplotRegression(fit, 'alternative_3SS_only')
plot_grid(p1, p2)
filename <- (paste(FIGURE_DIR,'SE_mean_by_ann_U_score.svg', sep = '') )
ggsave(filename, width = 6, height = 4)
## Fig 2
LIMITS = c(-50, 50)
#  
ggplot(filtered_tbl %>% filter(total_gapped_reads > 5, alt_3SS, threeSS_type_factor != 'AC', threeSS_type_factor != 'other'), aes(BP_3SS_dist , prp18_effect_FAnS_with_pseudo_counts ) ) + 
  geom_point(size = 2, alpha = 0.3, stroke = 2, aes(colour = threeSS_type)) +  scale_x_log10() + scale_y_log10() +
  facet_grid(threeSS_type_factor ~ .) +
  geom_abline(slope = 0, linetype = 'dotted') + geom_smooth(method='lm', size = .5, alpha = 0.1) +
  geom_vline(xintercept = 7, linetype = 'dotted') +
  geom_vline(xintercept = 10, linetype = 'dotted') +
  geom_text_repel(aes(label=ifelse(common_gene_name %in% genes_to_highlight & threeSS_type_factor %in% c('BG', 'HAT') & annotated_to_alternative_3SS_dist != 0 & upf1_prp18_mean_FAnS_with_pseudo_fraction > .1 , as.character(paste(common_gene_name, threeSS_seq, sep = ' ')),'')) , nudge_y = 4, size = 5, force = 1, fontface = 'italic') +
  geom_text_repel(aes(label=ifelse(common_gene_name %in% genes_to_highlight & substr(threeSS_type_factor, 2,3) == 'AG' & annotated_to_alternative_3SS_dist != 0 & upf1_prp18_mean_FAnS_with_pseudo_fraction > .1 , as.character(common_gene_name),'')) , nudge_y = 4, size = 5, force = 1, fontface = 'italic') 
  

filename <- (paste(FIGURE_DIR,'prp18_effect_FAnS_by_BP_dist_alt_only.svg', sep = '') )
ggsave(filename, width = 10, height = 12)


filtered_tbl_min_6_reads_max_50_alt_dist <- filtered_tbl_min_6_reads_min_one_thousandth_FAnS %>%
  mutate(annotated_to_alternative_3SS_dist_max_50 = ifelse(annotated_to_alternative_3SS_dist < -50, -50,
                                                           ifelse(annotated_to_alternative_3SS_dist > 50, 50,
                                                                  annotated_to_alternative_3SS_dist    )
  ) )


filtered_tbl_max_50_alt_dist <- filtered_tbl %>%
  mutate(annotated_to_alternative_3SS_dist_max_50 = ifelse(annotated_to_alternative_3SS_dist < -50, -50,
                                                           ifelse(annotated_to_alternative_3SS_dist > 50, 50,
                                                                  annotated_to_alternative_3SS_dist    )
  ) )

## Fig 1B
LIMITS = c(-50, 50)
genes_to_highlight <- c('PHO85', 'MUD1', 'SEC14', 'SPT14', 'MAF1', 'UBC12', 'BET4', 'YCL002C', 'NYV1', 'RAD14', 'RPL30', 'PRE3', 'RPL24B')

filtered_tbl_min_6_reads_max_50_alt_dist$annotated_to_alternative_3SS_dist_max_50
ggplot(filtered_tbl_max_50_alt_dist %>% filter(total_gapped_reads > 5, !alt_5SS, threeSS_type_factor != 'AC', threeSS_type_factor != 'other') , aes(annotated_to_alternative_3SS_dist_max_50 , prp18_effect_FAnS_with_pseudo_counts ) ) + 
  geom_point( size = 3, alpha = 0.3, stroke = 1, aes(colour = threeSS_type)) + scale_y_log10() +  
  scale_x_continuous( limits = LIMITS, breaks = seq(-50, 50, 10))+ 
  facet_grid(threeSS_type_factor ~ .) +
  geom_abline(slope = 0, linetype = 'dotted') + geom_smooth(method='lm', size = .5, alpha = 0.1) + 
  geom_vline(xintercept = 0, linetype = 'dotted') +
  geom_text_repel(aes(label=ifelse(common_gene_name %in% genes_to_highlight & threeSS_type_factor %in% c('BG', 'HAT') & annotated_to_alternative_3SS_dist != 0 & upf1_prp18_mean_FAnS_with_pseudo_fraction > .1 , as.character(paste(common_gene_name, threeSS_seq, sep = ' ')),'')) , nudge_y = 4, size = 5, force = 1, fontface = 'italic') +
  geom_text_repel(aes(label=ifelse(common_gene_name %in% genes_to_highlight & substr(threeSS_type_factor, 2,3) == 'AG' & annotated_to_alternative_3SS_dist != 0 & upf1_prp18_mean_FAnS_with_pseudo_fraction > .1 , as.character(common_gene_name),'')) , nudge_y = 4, size = 5, force = 1, fontface = 'italic') 


filename <- (paste(FIGURE_DIR,'prp18_effect_FAnS_by_alt_ann_dist_alt_only.svg', sep = '') )
ggsave(filename, width = 10, height = 12)

gghistogram(alternative_3SS_only %>% filter(total_gapped_reads > 5, substr(threeSS_seq, 2,3) == 'AG'),  "annotated_to_alternative_3SS_dist", binwidth = 20 ) + 
  facet_grid(threeSS_type_factor ~ .)

ggdensity(alternative_3SS_only %>% filter(total_gapped_reads > 5),  "annotated_to_alternative_3SS_dist", add = "median" ) + 
  facet_grid(threeSS_type_factor ~ .) +
  scale_x_continuous( limits = LIMITS, breaks = seq(-100, 300, 10)) 

filtered_tbl_min_6_reads_max_50_alt_dist$alt_5SS_alt_3SS
filtered_tbl_min_6_reads_max_50_alt_dist$annotated_to_alternative_3SS_dist_max_50

filtered_tbl_min_6_reads_max_50_alt_dist %>% filter( annotated_to_alternative_3SS_dist == 2 & threeSS_seq == 'GAG') %>% 
    select(common_gene_name, BP_3SS_dist, threeSS_seq, upf1_prp18_mean_FAnS, upf1_prp18_mean_SE, prp18_effect_FAnS ) %>% print(n=40) 

annotated_only %>% group_by(DS_3SS_1, DS_3SS_2) %>% count()

LIMITS = c(-10, 10)
ggplot(filtered_tbl_min_6_reads_max_50_alt_dist, aes( annotated_to_alternative_3SS_dist_max_50, prp18_effect_SE_with_pseudo_counts ) ) + 
  geom_jitter(width = .2, size = 2, alpha = 0.3, stroke = .5, aes(colour = threeSS_type)) + scale_y_log10() + 
  facet_grid(threeSS_type_factor ~ .) + scale_x_continuous( limits = LIMITS, breaks = seq(-10, 10, 2)) +
  geom_abline(slope = 0, linetype = 'dotted') +  
  geom_vline(xintercept = 0, linetype = 'dotted')
filename <- (paste(FIGURE_DIR,'prp18_SE_effect_by_ann_to_alt_dist_zoomed_in.svg', sep = '') )
ggsave(filename, width = 8, height = 6)

df <- filtered_tbl_min_6_reads %>% filter( annotated_to_alternative_3SS_dist > -11 & annotated_to_alternative_3SS_dist < 11)

## for annotated junctions
ggplot(annotated_only, aes(BP_3SS_dist , prp18_effect_SE_with_pseudo_counts ) ) + 
  geom_point(size = 2, alpha = 0.3, stroke = 2, aes(colour = threeSS_type, shape = alt_5SS_alt_3SS)) +  scale_x_log10() + scale_y_log10() + 
  facet_grid(threeSS_type ~ .) +
  geom_abline(slope = 0, linetype = 'dotted') + geom_smooth(method='lm', size = .5, alpha = 0.1) 
filename <- (paste(FIGURE_DIR,'prp18_effect_SE_by_BP_dist_ann_only.svg', sep = '') )
ggsave(filename, width = 8, height = 6)

ggplot(alternative_3SS_only, aes(accessibility, prp18_effect_SE_with_pseudo_counts ) ) + 
  geom_point(size = 2, alpha = 0.3, stroke = 2, aes(colour = threeSS_type, shape = alt_5SS_alt_3SS)) + scale_y_log10() + 
  facet_grid(threeSS_type ~ .) + 
  geom_abline(slope = 0, linetype = 'dotted') + geom_smooth(method='lm', size = .5, alpha = 0.1) 
filename <- (paste(FIGURE_DIR,'prp18_SE_effect_by_ann_to_alt_dist_alt_only.svg', sep = '') )
ggsave(filename, width = 8, height = 6)

## Fig S8
SIZE = 1
STROKE = .5

p1 <- ggplot(annotated_only, aes(BP_3SS_dist , upf1_mean_SE ) ) + 
  geom_point(size = SIZE, alpha = 0.3, stroke = STROKE, aes(colour = threeSS_type)) +  scale_x_log10() + scale_y_log10() + 
  facet_grid(threeSS_type ~ .) +
  geom_abline(slope = 0, linetype = 'dotted') + geom_smooth(method='lm', size = .5, alpha = 0.1) 

p2 <- ggplot(annotated_only, aes(BP_3SS_dist , upf1_prp18_mean_SE ) ) + 
  geom_point(size = SIZE, alpha = 0.3, stroke = STROKE, aes(colour = threeSS_type)) +  scale_x_log10() + scale_y_log10() + 
  facet_grid(threeSS_type ~ .) +
  geom_abline(slope = 0, linetype = 'dotted') + geom_smooth(method='lm', size = .5, alpha = 0.1) 

p3 <- ggplot(alternative_3SS_only, aes(BP_3SS_dist , upf1_mean_SE ) ) + 
  geom_point(size = SIZE, alpha = 0.3, stroke = STROKE, aes(colour = threeSS_type )) +  scale_x_log10() + scale_y_log10() + 
  facet_grid(threeSS_type ~ .) +
  geom_abline(slope = 0, linetype = 'dotted') + geom_smooth(method='lm', size = .5, alpha = 0.1) 

p4 <- ggplot(alternative_3SS_only, aes(BP_3SS_dist , upf1_prp18_mean_SE ) ) + 
  geom_point(size = SIZE, alpha = 0.3, stroke = STROKE, aes(colour = threeSS_type )) +  scale_x_log10() + scale_y_log10() + 
  facet_grid(threeSS_type ~ .) +
  geom_abline(slope = 0, linetype = 'dotted') + geom_smooth(method='lm', size = .5, alpha = 0.1) 

p5 <- ggplot(alternative_3SS_only, aes(annotated_to_alternative_3SS_dist, upf1_mean_SE ) ) + 
  geom_point(size = SIZE, alpha = 0.3, stroke = STROKE, aes(colour = threeSS_type)) + scale_y_log10() + 
  facet_grid(threeSS_type ~ .) + scale_x_continuous( limits = LIMITS) +
  geom_abline(slope = 0, linetype = 'dotted') + geom_smooth(method='lm', size = .5, alpha = 0.1) 

p6 <- ggplot(alternative_3SS_only, aes(annotated_to_alternative_3SS_dist, upf1_prp18_mean_SE ) ) + 
  geom_point(size = SIZE, alpha = 0.3, stroke = STROKE, aes(colour = threeSS_type)) + scale_y_log10() + 
  facet_grid(threeSS_type ~ .) + scale_x_continuous( limits = LIMITS) +
  geom_abline(slope = 0, linetype = 'dotted') + geom_smooth(method='lm', size = .5, alpha = 0.1) 

plot_grid(p1, p2, p3, p4, p5, p6, ncol = 2)
filename <- (paste(FIGURE_DIR,'SE_by_BP_dist_ann_and_alt_facet.svg', sep = '') )
ggsave(filename, width = 12, height = 10)


filtered_tbl_min_6_reads %>% filter(fiveSS_type == 'non-GT', annotated_to_alternative_5SS_dist != annotated_to_alternative_3SS_dist) %>% 
  select(prp18_effect_FAnS_with_pseudo_counts, prp18_effect_SE_with_pseudo_counts, annotated_to_alternative_5SS_dist, fiveSS_seq, annotated_to_alternative_3SS_dist, threeSS_seq, US_5SS_1, DS_3SS_1, US_5SS_2, DS_3SS_2, US_5SS_3, DS_3SS_3  ) %>%
  count()

filtered_tbl%>% filter(fiveSS_type == 'non-GT', annotated_to_alternative_5SS_dist != annotated_to_alternative_3SS_dist) %>% 
  select(prp18_effect_FAnS_with_pseudo_counts, prp18_effect_SE_with_pseudo_counts, annotated_to_alternative_5SS_dist, fiveSS_seq, annotated_to_alternative_3SS_dist, threeSS_seq, US_5SS_1, DS_3SS_1, US_5SS_2, DS_3SS_2, US_5SS_3, DS_3SS_3  ) %>%
  count()



filtered_tbl %>% filter(fiveSS_type == 'non-GT', annotated_to_alternative_5SS_dist != annotated_to_alternative_3SS_dist) %>% 
  select(common_gene_name, total_gapped_reads, total_reads_without_mismatch_near_junction, overhang_perfect_matches_upstream_junction, overhang_perfect_matches_downstream_junction, prp18_effect_FAnS_with_pseudo_counts, prp18_effect_SE_with_pseudo_counts, annotated_to_alternative_5SS_dist, fiveSS_seq, annotated_to_alternative_3SS_dist, threeSS_seq, US_5SS_1, DS_3SS_1, US_5SS_2, DS_3SS_2, US_5SS_3, DS_3SS_3  ) %>%
  write_tsv(paste(DIR, 'nonGT_fiveSS_events.xls'))

non_G_tbl <- filtered_tbl_min_6_reads %>% filter(threeSS_type == 'HAT' | threeSS_type == 'AC' |  threeSS_type == 'other', annotated_to_alternative_5SS_dist != annotated_to_alternative_3SS_dist)
non_G_tbl <- non_G_tbl %>% mutate( fiveSS1 = substr(fiveSS_seq, 1, 1), fiveSS2 = substr(fiveSS_seq, 2, 2),
                                   threeSS1 =  substr(threeSS_seq, 1, 1), threeSS2 =  substr(threeSS_seq, 2,2), threeSS3 =  substr(threeSS_seq, 3, 3)  )

non_G_tbl %>% filter( threeSS_seq == 'AGA') %>% group_by(annotated_to_alternative_3SS_dist) %>% tally() %>% arrange(desc(n))
counts <- non_G_tbl %>% group_by(threeSS_seq) %>% tally() %>% arrange(desc(n))  # %>% write_tsv(paste(DIR,'nonG_3SS_filter_1_2.txt') )
non_G_tbl %>% filter(threeSS_seq == 'TTG') %>% group_by(fiveSS_seq) %>% tally() %>% arrange(desc(n))
x <- expand.grid(rep(list(c('A', 'G', 'T', 'C')), 3))
x <- x %>% filter( Var3 != 'G')
combos <- as.tibble( do.call(paste0, x) ) %>% rename( "threeSS_seq" =  value )
df <- left_join(combos, counts, fill = 0) %>% arrange(desc(threeSS_seq))
ggplot( df, aes( as.factor(threeSS_seq), n)) + geom_col() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(paste(FIGURE_DIR, 'nonG_threeSS_seq_counts_filter_1_and_2.jpg'), width = 12, height = 8)


df <- filtered_tbl_min_6_reads %>% filter(threeSS_type == 'BG', annotated_to_alternative_5SS_dist != annotated_to_alternative_3SS_dist)
df <- df %>% mutate( fiveSS1 = substr(fiveSS_seq, 1, 1), fiveSS2 = substr(fiveSS_seq, 2, 2),
                     threeSS1 =  substr(threeSS_seq, 1, 1), threeSS2 =  substr(threeSS_seq, 2,2), threeSS3 =  substr(threeSS_seq, 3, 3)  )

## BEGIN FIG. 3B
## AG allowed
aggregate_FAnS <- filtered_tbl %>% filter(alt_3SS == TRUE & alt_5SS == FALSE & substr(threeSS_seq, 2,3) == "AG") %>% group_by(threeSS_seq) %>% tally(upf1_prp18_mean_FAnS) %>% arrange(desc(n)) ##  %>% write_tsv(paste(DIR,'aggreagate_FAnS_by_3SS_filter_1.txt') )
aggregate_SE <- filtered_tbl %>% filter(alt_3SS == TRUE & alt_5SS == FALSE & substr(threeSS_seq, 2,3) == "AG") %>% group_by(threeSS_seq) %>% tally(upf1_prp18_mean_SE) %>% arrange(desc(n)) ##  %>% write_tsv(paste(DIR,'aggreagate_SE_by_3SS_filter_1.txt') )
aggregate_upf1_FAnS <- filtered_tbl %>% filter(alt_3SS == TRUE & alt_5SS == FALSE & substr(threeSS_seq, 2,3) == "AG") %>% group_by(threeSS_seq) %>% tally(upf1_mean_FAnS) %>% arrange(desc(n)) ##  %>% write_tsv(paste(DIR,'aggreagate_FAnS_by_3SS_filter_1.txt') )
aggregate_upf1_SE <- filtered_tbl %>% filter(alt_3SS == TRUE & alt_5SS == FALSE & substr(threeSS_seq, 2,3) == "AG") %>% group_by(threeSS_seq) %>% tally(upf1_mean_SE) %>% arrange(desc(n)) ##  %>% write_tsv(paste(DIR,'aggreagate_SE_by_3SS_filter_1.txt') )

x <- expand.grid(rep(list(c('A', 'G', 'T', 'C')), 3))
combos <- as.tibble( do.call(paste0, x) ) %>% rename( "threeSS_seq" =  value ) %>% arrange(threeSS_seq) 

left_join(combos,aggregate_SE,  fill = 0) %>% arrange(n)  %>% print(n=40)## %>% filter( n == 0 ) 

df <- left_join(combos, aggregate_FAnS,  fill = 0) %>% arrange(desc(n)) %>% slice(1:20) # %>% filter(aggregate_FAnS > 1)
p1 <- ggplot( df, aes( reorder(threeSS_seq, -n, sum), n)) + geom_col() + ggtitle("aggregate prp18upf1 FAnS on alt 3SS") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) #+ scale_y_log2()

df <- left_join(combos,aggregate_SE,  fill = 0) %>% arrange(desc(n)) %>% slice(1:20)
p2<- ggplot( df, aes( reorder(threeSS_seq, -n, sum), n))  + geom_col() + ggtitle("aggregate prp18upf1 SE on alt 3SS") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) #+ scale_y_log2()

df <- left_join(combos,aggregate_upf1_FAnS,  fill = 0) %>% arrange(desc(n)) %>% slice(1:20)
p3 <-ggplot( df, aes( reorder(threeSS_seq, -n, sum), n))  + geom_col() + ggtitle("aggregate upf1 FAnS on alt 3SS")+ ylim(0, 25) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) #+ scale_y_log2()

df <- left_join(combos,aggregate_upf1_SE,  fill = 0) %>% arrange(desc(n)) %>% slice(1:20)
p4 <-ggplot( df, aes( reorder(threeSS_seq, -n, sum), n))  + geom_col() + ggtitle("aggregate upf1 SE on alt 3SS")+ 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) #+ scale_y_log2()

plot_grid(p1,p3,p2,p4, ncol = 1)
ggsave(paste(FIGURE_DIR, 'aggregate_FAnS_SE_prp18_effect_AG.svg'), width = 8, height = 8)

## non-AG allowed
aggregate_FAnS <- filtered_tbl %>% filter(alt_3SS == TRUE & alt_5SS == FALSE & threeSS_type %in% c("BG", "non-G") ) %>% group_by(threeSS_seq) %>% tally(upf1_prp18_mean_FAnS) %>% arrange(desc(n)) ##  %>% write_tsv(paste(DIR,'aggreagate_FAnS_by_3SS_filter_1.txt') )
aggregate_SE <- filtered_tbl %>% filter(alt_3SS == TRUE & alt_5SS == FALSE &  threeSS_type %in% c("BG", "non-G")) %>% group_by(threeSS_seq) %>% tally(upf1_prp18_mean_SE) %>% arrange(desc(n)) ##  %>% write_tsv(paste(DIR,'aggreagate_SE_by_3SS_filter_1.txt') )
aggregate_upf1_FAnS <- filtered_tbl %>% filter(alt_3SS == TRUE & alt_5SS == FALSE &  threeSS_type %in% c("BG", "non-G")) %>% group_by(threeSS_seq) %>% tally(upf1_mean_FAnS) %>% arrange(desc(n)) ##  %>% write_tsv(paste(DIR,'aggreagate_FAnS_by_3SS_filter_1.txt') )
aggregate_upf1_SE <- filtered_tbl %>% filter(alt_3SS == TRUE & alt_5SS == FALSE &  threeSS_type %in% c("BG", "non-G")) %>% group_by(threeSS_seq) %>% tally(upf1_mean_SE) %>% arrange(desc(n)) ##  %>% write_tsv(paste(DIR,'aggreagate_SE_by_3SS_filter_1.txt') )
x <- expand.grid(rep(list(c('A', 'G', 'T', 'C')), 3))
combos <- as.tibble( do.call(paste0, x) ) %>% rename( "threeSS_seq" =  value ) %>% arrange(threeSS_seq) 

left_join(combos,aggregate_SE,  fill = 0) %>% arrange(n)  %>% print(n=40)## %>% filter( n == 0 ) 

df <- left_join(combos, aggregate_FAnS,  fill = 0) %>% arrange(desc(n)) %>% slice(1:20) # %>% filter(aggregate_FAnS > 1)
p1 <- ggplot( df, aes( reorder(threeSS_seq, -n, sum), n)) + geom_col() + ggtitle("aggregate prp18upf1 FAnS on alt 3SS") + ylim(0, 1.5) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) #+ scale_y_log2()

df <- left_join(combos,aggregate_SE,  fill = 0) %>% arrange(desc(n)) %>% slice(1:20)
p2<- ggplot( df, aes( reorder(threeSS_seq, -n, sum), n))  + geom_col() + ggtitle("aggregate prp18upf1 SE on alt 3SS")+ ylim(0, 0.5) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) #+ scale_y_log2()

df <- left_join(combos,aggregate_upf1_FAnS,  fill = 0) %>% arrange(desc(n)) %>% slice(1:20)
p3 <-ggplot( df, aes( reorder(threeSS_seq, -n, sum), n))  + geom_col() + ggtitle("aggregate upf1 FAnS on alt 3SS")+ ylim(0, 1.5) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) #+ scale_y_log2()

df <- left_join(combos,aggregate_upf1_SE,  fill = 0) %>% arrange(desc(n)) %>% slice(1:20)
p4 <-ggplot( df, aes( reorder(threeSS_seq, -n, sum), n))  + geom_col() + ggtitle("aggregate upf1 SE on alt 3SS")+ ylim(0, 0.5) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) #+ scale_y_log2()

plot_grid(p1,p3,p2,p4, ncol = 1)
ggsave(paste(FIGURE_DIR, 'aggregate_FAnS_SE_prp18_effect_nonAG.svg'), width = 8, height = 8)

## END FIG. 3B

## BEGIN FIG. 4D
df <- alternative_3SS_only ##  %>% filter(alt_3SS == TRUE)
my_comparisons <- list( c("non-G", "YAG"), c("non-G", "RAG"),c("non-G", "BG"), c("RAG", "YAG"), c("BG", "YAG"), c("BG", "RAG") )
ggboxplot(df, x = "threeSS_type_factor", y = "accessibility",
          color = "threeSS_type_factor")+ 
  geom_jitter(width = 0.1, size = 1, alpha = 0.3) +
  stat_compare_means(comparisons = my_comparisons) # Add pairwise comparisons p-value
#  stat_compare_means(label.y = 1.5)     # Add global p-value
filename <- (paste(FIGURE_DIR,'alternative_3SS_only_accessibility_by_3SS_motif.svg', sep = '') )
ggsave(filename, width = 4, height = 4)
## END FIG. 4D

## Fig S6
## BEGIN normalized read counts for replicate #1 of upf1 and upf1/prp18
p1 <- ggplot(annotated_only, aes(upf1_mean_reads, upf1_prp18_mean_reads) ) + geom_point( aes(colour = factor(threeSS_type)), size = .5, alpha = .5 )  + 
  geom_abline(intercept = 0, slope = 1, linetype = 'dotted', alpha = .5) +
  scale_x_log10( limits = c(0.1, 100000)) + scale_y_log10(limits = c(0.1, 100000))
p2 <- ggplot(annotated_only, aes(upf1_mean_intron_reads, upf1_prp18_mean_intron_reads) ) + geom_point( aes(colour = factor(threeSS_type)), size = .5, alpha = .5 ) + 
  geom_abline(intercept = 0, slope = 1, linetype = 'dotted', alpha = .5) +
  scale_x_log10( limits = c(0.1, 10000)) + scale_y_log10(limits = c(0.1, 10000))
p3 <- ggplot(alternative_3SS_only, aes(upf1_mean_reads, upf1_prp18_mean_reads) ) + geom_point( aes(colour = factor(threeSS_type)), size = .5, alpha = .5)  + 
  geom_abline(intercept = 0, slope = 1, linetype = 'dotted', alpha = .5) +
  scale_x_log10( limits = c(0.1, 10000)) + scale_y_log10(limits = c(0.1, 10000))
p4 <- ggplot(alternative_5SS_only, aes(upf1_mean_reads, upf1_prp18_mean_reads) ) + geom_point( aes(colour = factor(fiveSS_type)), size = .5, alpha = .5  ) + 
  geom_abline(intercept = 0, slope = 1, linetype = 'dotted', alpha = .5) +
  scale_x_log10( limits = c(0.1, 10000)) + scale_y_log10(limits = c(0.1, 10000))
plot_grid(p1, p2, p3, p4, ncol=2)
ggsave(filename = paste(DIR, 'mean_reads_rpm_upf1_by_upf1_prp18.svg'), device = 'svg', width = 12, height = 6) 
## END normalized read counts for replicate #1 of upf1 and upf1/prp18

ggplot(annotated_only, aes(BP_5SS_dist, annotated_U_score, color = RPG )) + geom_point()
ggplot(annotated_only, aes(log2_BP_3SS_dist, annotated_U_score, color = RPG )) + geom_point()
annotated_only %>% filter( BP_3SS_dist > 100)

## Very similar to figure 4b, shows BP_3SS distance by prp18 effect, with strong correlation for RAG and YAG, RPG vs non-RPG show no major difference
ggplot(alternative_3SS_only %>% filter(total_gapped_reads > 5), aes(x = BP_3SS_dist, y = prp18_effect_SE_with_pseudo_counts  , colour=  RPG) ) + geom_point( ) + scale_x_continuous(trans='log2') + scale_y_continuous(trans='log2') + facet_wrap(~threeSS_type) # 
ggsave( paste(DIR,'alt_3SS_prp18_effect_vs_BP_3SS_dist.eps', sep = ''), height = HEIGHT, width = WIDTH )

## replicate counts for alt. 3SS
alt_3SS_junction_log2_cts <- tbl %>% 
  filter(alt_5SS == FALSE & alt_3SS == TRUE) %>%
  select( contains(':unique_junction_reads') ) %>%
  log10()
eps( paste(DIR,'alt_3SS_junction_log2_cts.eps', sep = ''), height = HEIGHT, width = WIDTH)
g <- ggpairs(alt_3SS_junction_log2_cts)
#warnings()
print(g)
dev.off()

## BEGIN DS 3SS nt for upf1
my_comparisons <- list( c("A", "T"), c("A", "C"),c("A", "G"),c("T", "C") ,c("T", "G") ,c("C", "G")  )
p1 <- ggboxplot(alternative_3SS_only, x = "DS_3SS_1", y = "upf1_mean_SE",
                width = BOX_WIDTH) + 
  geom_jitter(width = 0.1, size = 1, alpha = 0.3) + scale_y_log10() +
  #geom_text(label = alternative_3SS_only$common_gene_name, nudge_x  = NUDGE_X, check_overlap = TRUE, size = SIZE) +
  stat_compare_means(comparisons = my_comparisons) +
  ggtitle("DS_3SS_1" )

p2 <- ggboxplot(alternative_3SS_only, x = "DS_3SS_2", y = "upf1_mean_SE",
                width = BOX_WIDTH) + 
  geom_jitter(width = 0.1, size = 1, alpha = 0.3) + scale_y_log10() +
  #geom_text(label = alternative_3SS_only$common_gene_name, nudge_x  = NUDGE_X, check_overlap = TRUE, size = SIZE) +
  stat_compare_means(comparisons = my_comparisons) +
  ggtitle("DS_3SS_2" )

p3 <- ggboxplot(alternative_3SS_only, x = "DS_3SS_3", y = "upf1_mean_SE",
                width = BOX_WIDTH) + 
  geom_jitter(width = 0.1, size = 1, alpha = 0.3) + scale_y_log10() +
  #geom_text(label = alternative_3SS_only$common_gene_name, nudge_x  = NUDGE_X, check_overlap = TRUE, size = SIZE) +
  stat_compare_means(comparisons = my_comparisons) +
  ggtitle("DS_3SS_3" )
plot_grid(p1, p2, p3, nrow =1 )
filename <- (paste(FIGURE_DIR,'ds_3SS_nt_upf1.svg', sep = '') )
ggsave(filename, width = 9, height = 6)
## END DS 3SS nt for upf1
alternative_3SS_only$log10_prp18_effect_FAnS_with_pseudo_counts
## DS 3SS nt for upf1 prp18
p1 <- ggboxplot(alternative_3SS_only, x = "DS_3SS_1", y = "log10_prp18_effect_FAnS_with_pseudo_counts",
                width = BOX_WIDTH) + 
  geom_jitter(width = 0.1, size = 1, alpha = 0.3) + scale_y_log10() +
  #geom_text(label = alternative_3SS_only$common_gene_name, nudge_x  = NUDGE_X, check_overlap = TRUE, size = SIZE) +
  stat_compare_means(comparisons = my_comparisons) +
  ggtitle("DS_3SS_1" )

p2 <- ggboxplot(alternative_3SS_only, x = "DS_3SS_2", y = "log10_prp18_effect_FAnS_with_pseudo_counts",
                width = BOX_WIDTH) + 
  geom_jitter(width = 0.1, size = 1, alpha = 0.3) + scale_y_log10() +
  #geom_text(label = alternative_3SS_only$common_gene_name, nudge_x  = NUDGE_X, check_overlap = TRUE, size = SIZE) +
  stat_compare_means(comparisons = my_comparisons) +
  ggtitle("DS_3SS_2" )

p3 <- ggboxplot(alternative_3SS_only, x = "DS_3SS_3", y = "log10_prp18_effect_FAnS_with_pseudo_counts",
                width = BOX_WIDTH) + 
  geom_jitter(width = 0.1, size = 1, alpha = 0.3) + scale_y_log10() +
  #geom_text(label = alternative_3SS_only$common_gene_name, nudge_x  = NUDGE_X, check_overlap = TRUE, size = SIZE) +
  stat_compare_means(comparisons = my_comparisons) +
  ggtitle("DS_3SS_3" )
plot_grid(p1, p2, p3, nrow =1 )
filename <- (paste(FIGURE_DIR,'ds_3SS_nt_upf1_prp18.svg', sep = '') )
ggsave(filename, width = 9, height = 6)
## END DS 3SS nt for upf1 prp18

## MOTIF ANALYSIS FIG 3C and 3D
filtered_tbl  %>% group_by(alt_5SS, alt_3SS, DS_3SS_1) %>% tally()

filtered_tbl %>% filter(alt_3SS == TRUE) %>% group_by(threeSS_type, DS_3SS_1) %>% tally()
annotated_only %>% group_by(DS_3SS_1) %>% tally()

alternative_3SS_only %>% filter(total_gapped_reads > 5) %>% group_by(US_5SS_1) %>% tally()
alternative_5SS_only %>% filter(total_gapped_reads > 5) %>% group_by(US_5SS_1) %>% tally()
annotated_only %>% group_by(US_5SS_1) %>% tally()

alternative_3SS_only %>% filter(total_gapped_reads > 5) %>% group_by(US_5SS_2) %>% tally()
alternative_5SS_only %>% filter(total_gapped_reads > 5) %>% group_by(US_5SS_2) %>% tally()
annotated_only %>% group_by(US_5SS_2) %>% tally()

alternative_3SS_only %>% filter(total_gapped_reads > 5) %>% group_by(US_5SS_3) %>% tally()
alternative_5SS_only %>% filter(total_gapped_reads > 5) %>% group_by(US_5SS_3) %>% tally()
annotated_only %>% group_by(US_5SS_3) %>% tally()

alternative_3SS_only %>% group_by(DS_3SS_2) %>% tally()
annotated_only %>% group_by(DS_3SS_2) %>% tally()

alternative_3SS_only %>% group_by(DS_3SS_3) %>% tally()
annotated_only %>% group_by(DS_3SS_3) %>% tally()

# '{'A':32.79, 'C':19.15, 'G':20.43, 'T':27.63}' 
for (fiveSS_motif in c('GTATGT', 'GT',  'non-GT') ) {
  seq_filename <-paste(DIR, fiveSS_motif, '_alt5SS_intron_filtered.txt', sep = '')  
  ## run the command printed on command line to generate logo files
  print( paste("weblogo -S .3 -X FALSE -n 100 -i -10 -t ", fiveSS_motif , " < ", seq_filename, ' >', DIR, fiveSS_motif, 'filtered.eps', sep = ''  ) )
  filtered_tbl %>% 
    filter(alt_5SS == TRUE, alt_3SS == FALSE, fiveSS_type == fiveSS_motif) %>%
    mutate(combined_seq = paste(upstream_5SS, downstream_3SS, sep = '') ) %>%
    select( combined_seq ) %>%
    write_delim( seq_filename, col_names = FALSE )
}


for (threeSS_motif in c('TAG', 'CAG',  'AAG', 'GAG', 'HAT', 'BG') ) {
  seq_filename <-paste(DIR, threeSS_motif, '_alt3SS_intron_filtered.txt', sep = '')  
  ## run the command printed on command line to generate logo files
  print( paste("weblogo -X FALSE -n 100 -i -10 -t ", threeSS_motif , " < ", seq_filename, ' >', DIR, threeSS_motif, 'filtered.eps', sep = ''  ) )
  filtered_tbl %>% 
    filter(alt_3SS == FALSE, alt_5SS == FALSE, fiveSS_type == 'GTATGT', threeSS_type == threeSS_motif) %>%
    mutate(combined_seq = paste(fiveSS_seq, downstream_5SS, upstream_3SS,  threeSS_seq, sep = '') ) %>%
    select( combined_seq ) %>%
    write_delim( seq_filename, col_names = FALSE )
}

for (threeSS_motif in c('TAG', 'CAG',  'AAG', 'GAG', 'HAT', 'BG') ) {
  seq_filename <-paste(DIR, threeSS_motif, '_ann3SS_intron_filtered.txt', sep = '')  
  ## run the command printed on command line to generate logo files
  print( paste("weblogo -X FALSE -n 100 -i -10 -t ", threeSS_motif , " < ", seq_filename, ' >', DIR, threeSS_motif, 'filtered.eps', sep = ''  ) )
  filtered_tbl %>% 
    filter(alt_3SS == FALSE, alt_5SS == FALSE, threeSS_type == threeSS_motif) %>%
    mutate(combined_seq = paste(upstream_5SS, downstream_3SS, sep = '') ) %>%
    select( combined_seq ) %>%
    write_delim( seq_filename, col_names = FALSE )
}

for (threeSS_type_ann_factor in c('annotated_YAG', 'alternative_YAG', 'AAG', 'GAG', 'BG', 'HAT', 'other', 'alternative_5SS_YAG') ) {
  seq_filename <-paste(DIR, threeSS_type_ann_factor, '_intron_filtered.txt', sep = '')  
  ## run the command printed on command line to generate logo files
  print( paste("weblogo -S .3 -X FALSE -n 100 -i -10 -t ", threeSS_type_ann_factor , " < ", seq_filename, ' >', DIR, threeSS_type_ann_factor, 'filtered.eps', sep = ''  ) )
  filtered_tbl %>% 
    filter(fiveSS_type == 'GTATGT', threeSS_type_annotated_factor == threeSS_type_ann_factor) %>%
    mutate(combined_seq = paste(upstream_5SS, downstream_3SS, sep = '') ) %>%
    select( combined_seq ) %>%
    write_delim( seq_filename, col_names = FALSE )
}
## fiveSS_seq, downstream_5SS, upstream_3SS,  threeSS_seq

seq_name <- 'ann_3SS_intron'
seq_filename <-paste(DIR, seq_name, '.txt', sep = '') 
filtered_tbl %>% 
  filter(annotated) %>%
  mutate(combined_seq = paste(fiveSS_seq, downstream_5SS, upstream_3SS,  threeSS_seq, sep = '') ) %>%
  select( combined_seq ) %>%
  write_delim( seq_filename, col_names = FALSE )
paste("weblogo -X FALSE  -n 100 -i -10 -t ", seq_name , " < ", seq_filename, ' >', DIR, seq_name, '.eps', sep = ''  )

seq_name <- 'alt_3SS_intron'
seq_filename <-paste(DIR, seq_name, '.txt', sep = '') 
filtered_tbl %>% 
  filter(alt_3SS == TRUE) %>%
  mutate(combined_seq = paste(fiveSS_seq, downstream_5SS, upstream_3SS,  threeSS_seq, sep = '') ) %>%
  select( combined_seq ) %>%
  write_delim( seq_filename, col_names = FALSE )
paste("weblogo -X FALSE -n 100 -i -10 -t ", seq_name , " < ", seq_filename, ' >', DIR, seq_name, '.eps', sep = ''  )

## BEGIN plot SLOW, FAST, WT splicing data (2018 Genome Research paper)
slow_fast_RNAP_tbl <- read_tsv(paste(DIR,'slow_fast_RNAP_data.txt', sep = '') )
#SEF <- slow_fast_RNAP_tbl%>% select( SEF_FAST_A, SEF_FAST_B, SEF_SLOW_A, SEF_SLOW_B, SEF_WT_A, SEF_WT_B)
slow_fast_RNAP_tbl <- slow_fast_RNAP_tbl %>% mutate( SEF_FAST = (SEF_FAST_A + SEF_FAST_B )/2, SEF_SLOW = (SEF_SLOW_A + SEF_SLOW_B )/2,SEF_WT = (SEF_WT_A + SEF_WT_B )/2, chromosome = paste('chr',chromosome, sep = '') ) %>%
  rename( adjusted_start = start, adjusted_end = end, adjusted_strand = strand)
slow_fast_RNAP_tbl <- slow_fast_RNAP_tbl %>% mutate( fast_effect = SEF_FAST / SEF_WT, slow_effect = SEF_SLOW / SEF_WT )
slow_fast_RNAP_tbl <- slow_fast_RNAP_tbl %>% mutate( log10_fast_effect = log10(fast_effect),log10_slow_effect = log10(slow_effect))
slow_fast_RNAP_tbl <- slow_fast_RNAP_tbl %>% mutate( log10_SEF_WT = log10(SEF_WT), log10_SEF_FAST = log10(SEF_FAST), log10_SEF_SLOW = log10(SEF_SLOW) )
slow_fast_RNAP_tbl <- slow_fast_RNAP_tbl %>% left_join(filtered_tbl) #, by = c("chromosome", adjusted_start, adjusted_en )  )
slow_fast_RNAP_tbl$
  
  slow_fast_RNAP_tbl  %>% filter(RPG) %>% select(log10_fast_effect, log10_slow_effect, log10_prp18_effect_SE_with_pseudo_counts, log10_prp18_effect_FAnS_with_pseudo_counts )%>% ggpairs() 
# slow_fast_RNAP_tbl %>% mutate(fast_avg = )
slow_fast_RNAP_tbl  %>% filter(RPG) %>% select(SEF_WT, SEF_SLOW, SEF_FAST, upf1_mean_FAnS, upf1_prp18_mean_FAnS, upf1_mean_SE, upf1_prp18_mean_SE )%>% ggpairs() 

ggplot(slow_fast_RNAP_tbl, aes(upf1_mean_FAnS, SEF_WT)) + geom_point() + scale_x_log10() + scale_y_log10() + 
  geom_text_repel(aes(label=ifelse( upf1_mean_FAnS > .1, as.character(common_gene_name),'')) , nudge_x = .5, size = 5, force = 1, fontface = 'italic') +
  geom_abline(intercept = 0, slope = 1, linetype = 'dashed', alpha = .5) + geom_smooth(method='lm', size = .5, alpha = 0.1) 
filename <- (paste(FIGURE_DIR,'alt_3SS_FAnS_compared_to_FAST_SLOW_WT_SPLICING_PAPER.eps', sep = '') )
ggsave(filename, width = 5, height = 5)
p1 <- ggplotRegression( lm(log10_upf1_mean_FAnS_with_pseudo_fraction  ~ log10_SEF_WT, data = slow_fast_RNAP_tbl), 'alt_3SS_FAnS')
p2 <- ggplotRegression( lm(log10_upf1_mean_FAnS_with_pseudo_fraction ~ log10_SEF_FAST, data = slow_fast_RNAP_tbl), 'alt_3SS_FAnS')
p3 <- ggplotRegression( lm(log10_upf1_mean_FAnS_with_pseudo_fraction ~ log10_SEF_SLOW, data = slow_fast_RNAP_tbl), 'alt_3SS_FAnS')
p4 <- ggplotRegression( lm(log10_upf1_prp18_mean_FAnS_with_pseudo_fraction  ~ log10_SEF_WT, data = slow_fast_RNAP_tbl), 'alt_3SS_FAnS')
p5 <- ggplotRegression( lm(log10_upf1_prp18_mean_FAnS_with_pseudo_fraction ~ log10_SEF_FAST, data = slow_fast_RNAP_tbl), 'alt_3SS_FAnS')
p6 <- ggplotRegression( lm(log10_upf1_prp18_mean_FAnS_with_pseudo_fraction ~ log10_SEF_SLOW, data = slow_fast_RNAP_tbl), 'alt_3SS_FAnS')
plot_grid(p1, p2, p3, p4, p5, p6)
## END plot SLOW, FAST, WT splicing data (2018 Genome Research paper)

