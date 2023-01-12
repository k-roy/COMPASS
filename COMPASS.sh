#!/bin/bash

##################### VARIABLE AND PATH DEFINITIONS BEGIN #####################
NUM_THREADS=32 
# Change to allowable number of threads for your system.

READS_TO_PROCESS=-1
# Subsample this many reads to test the pipeline.
# Set to -1 for all reads after tests pass.

# change $HOME variable to your desired directory
COMPASS_DIR=$HOME"/COMPASS/"
REFERENCE_DIR=$COMPASS_DIR"genome_references/"
FASTA=$REFERENCE_DIR"S288C_reference_sequence_R64-2-1_20150113_reformatted_chromosome_names.fasta"
GFF=$REFERENCE_DIR"S288C_reference_sequence_R64-2-1_20150113_reformatted_chromosome_names.gff"
GTF=$REFERENCE_DIR"S288C_reference_sequence_R64-2-1_20150113_reformatted_chromosome_names.gtf"
INTRONS_FILE=$REFERENCE_DIR"saccharomyces_cerevisiae_R64-2-1_20150113_introns.tsv"
ALIGNERS_FILE=$COMPASS_DIR"aligners.txt"

## ASSIGN AND CREATE SUBDIRECTORIES
NUMBERED_READS_DIR=$COMPASS_DIR"numbered_reads_fastq/"
mkdir $NUMBERED_READS_DIR
TRIMMED_DIR=$COMPASS_DIR"trimmed_fastq/"
mkdir $TRIMMED_DIR
SEPARATE_ALIGNMENTS_DIR=$COMPASS_DIR"separate_alignments/"
mkdir $SEPARATE_ALIGNMENTS_DIR
OPTIMAL_ALIGNMENTS_DIR=$COMPASS_DIR"COMPASS_alignments/"
mkdir $OPTIMAL_ALIGNMENTS_DIR
JUNCTION_ALIGNMENTS_DIR=$COMPASS_DIR"junction_read_alignments/"
mkdir $JUNCTION_ALIGNMENTS_DIR
COMPASS_JUNCTIONS_DIR=$COMPASS_DIR"COMPASS_junctions/"
mkdir $COMPASS_JUNCTIONS_DIR
LOG_DIR=$COMPASS_DIR"log/"
mkdir $LOG_DIR

###################### VARIABLE AND PATH DEFINITIONS END ######################

conda activate compass
gffread -T --force-exons --gene2exon $GFF -o $GTF $GFF

####################### BATCH PROCESSING COMMANDS BEGIN #######################
SAMPLE_SUFFIX="_name_sorted.bam"

# batch processing prp18/Roy et al.
ACCESSION="PRJNA544962"
READ_LENGTH=150
for i in {9130287..9130292};

do SAMPLE="SRR"$i

echo "READS_TO_PROCESS: "$READS_TO_PROCESS
echo "starting process_reads_and_align.sh for "$SAMPLE;
sh process_reads_and_align.sh -C "$COMPASS_DIR" -R "$REFERENCE_DIR" \
-S "$SAMPLE" -F "$FASTA" -G "$GTF" -P "$NUM_THREADS" \
-L $READ_LENGTH -Z $READS_TO_PROCESS \
> $LOG_DIR$SAMPLE"_COMPASS_process_reads_and_align.log" 2>&1

# echo "starting compare_splice_junctions_from_multiple_aligners.py for "$SAMPLE;
# python -u compare_splice_junctions_from_multiple_aligners.py \
# "$SEPARATE_ALIGNMENTS_DIR" "$OPTIMAL_ALIGNMENTS_DIR" "$FASTA" "$INTRONS_FILE" \
# "$ALIGNERS_FILE" "$SAMPLE" "$SAMPLE_SUFFIX" "$READS_TO_PROCESS" \
# > $LOG_DIR$SAMPLE"_COMPASS.log" 2>&1

done
####################### BATCH PROCESSING COMMANDS END #######################
