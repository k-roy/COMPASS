#!/bin/bash

######################### HARD DRIVE SPACE WARNING ############################
# The COMPASS directory needs to have substantial free space for the alignment files,
# on the order of 10X the size of the fastq.gz files when run on many aligners. 
# These can be deleted at the end of the pipeline for each sample to avoid accumulating excess storage.

##################### VARIABLE AND PATH DEFINITIONS #####################
NUM_THREADS=32 
# Change to allowable number of threads for your system.

READS_TO_PROCESS=-1
# Subsample this many reads to test the pipeline (typically 100,000 reads).
# Set to -1 for all reads after tests pass.

# change $HOME variable to your desired directory
COMPASS_DIR="$HOME/COMPASS/"

# download genome references and put into this folder
REFERENCE_DIR=$COMPASS_DIR"genome_references/"

# path to samfixcigar java program
SAMFIXCIGAR=$COMPASS_DIR"jvarkit/dist/samfixcigar.jar"

# S. cerevisiae genome reference
GENOME_VERSION="S288C_reference_sequence_R64-2-1_20150113_reformatted_chromosome_names"

FASTA=$REFERENCE_DIR$GENOME_VERSION".fasta"
GFF=$REFERENCE_DIR$GENOME_VERSION".gff"
GTF=$REFERENCE_DIR$GENOME_VERSION".gtf"
NUM_THREADS=16
MIN_INTRON_LENGTH=20
MAX_INTRON_LENGTH=2000 # typically set to 200,000 for human introns

ACCESSION= # set to "NA" for datasets not on SRA
READ_LENGTH=100
READS_TO_PROCESS=-1 # 1000000 # 
# Subsample this many reads if SUBSAMPLE_READS is set to true.
# set to -1 for all reads

ALIGNERS_FILE=$COMPASS_DIR"sample_aligner_info.tsv"
HISAT2_GENOME_DIR=$REFERENCE_DIR"HISAT2_annotated_index"

# the introns file is needed for HISAT2
INTRONS_FILE=$REFERENCE_DIR"saccharomyces_cerevisiae_R64-2-1_20150113_introns.tsv"

# The sample_aligner_info.txt file tells COMPASS which alignment programs to use from among
# BBMap, STAR (both default and noncanonical splicing modes), 
# HISAT2 (both default and noncanonical splicing modes), Magic-BLAST, and GSNAP.
ALIGNERS_FILE=$COMPASS_DIR"sample_aligner_info.txt"

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

################# ACTIVATE COMPASS ENVIRONMENT IN CONDA ###################
conda activate compass

# gffread is needed to convert the S. cerevisiae GFF file from SGD to GTF format
gffread -T --force-exons --gene2exon $GFF -o $GTF $GFF

####################### SAMPLE PROCESSING COMMANDS #######################
echo "File name is "$0 # holds the current script
echo "Sample name is "$1
SAMPLE=$1
SAMPLE_NAME=$SAMPLE

if [[ $READS_TO_PROCESS -lt 1 ]];
then
echo "processing all reads for "$SAMPLE
else
echo "subsampling "$READS_TO_PROCESS" reads for "$SAMPLE
SAMPLE=$SAMPLE"_subsampled"
fi

echo "READS_TO_PROCESS: "$READS_TO_PROCESS
echo "starting process_reads_and_align.sh for "$SAMPLE_NAME;
sh process_reads_and_align.sh -C "$COMPASS_DIR" -R "$REFERENCE_DIR" \
-S "$SAMPLE_NAME" -F "$FASTA" -G "$GTF" -P $NUM_THREADS \
-M $MIN_INTRON_LENGTH -N $MAX_INTRON_LENGTH -A $SAMFIXCIGAR \
-L $READ_LENGTH -Z $READS_TO_PROCESS \
> $LOG_DIR$SAMPLE"_process_reads_and_align.log" 2>&1

echo "starting compare_splice_junctions_from_multiple_aligners.py for "$SAMPLE;
python compare_splice_junctions_from_multiple_aligners.py \
"$COMPASS_DIR" "$REFERENCE_DIR" "$SAMPLE" "$FASTA" "$INTRONS_FILE" "$NUM_THREADS" \
"$MIN_INTRON_LENGTH" "$MAX_INTRON_LENGTH" "$ALIGNERS_FILE" "$READS_TO_PROCESS" \
> $LOG_DIR$SAMPLE"_compare_splice_junctions_from_multiple_aligners_sh_output.log" 2>&1

echo "starting analyze_exonic_and_intronic_sequence.py for "$SAMPLE;
python analyze_exonic_and_intronic_sequence.py \
"$COMPASS_DIR" "$REFERENCE_DIR" "$SAMPLE" "$FASTA" "$INTRONS_FILE" "$NUM_THREADS" \
"$MIN_INTRON_LENGTH" "$MAX_INTRON_LENGTH" "$ALIGNERS_FILE" "$READS_TO_PROCESS" \
> $LOG_DIR$SAMPLE"_analyze_exonic_and_intronic_sequence_sh_output.log" 2>&1

echo "starting create_splice_site_bed.py for "$SAMPLE;
python create_splice_site_bed.py $COMPASS_JUNCTIONS_DIR $SAMPLE

echo "starting add_unspliced_read_counts_to_junctions.py for "$SAMPLE;
for strand in "plus" "minus";
do echo $strand;
for splice_site in "fiveSS" "threeSS";
do echo $splice_site;
BED="$COMPASS_JUNCTIONS_DIR""$SAMPLE""_"$strand"_strand_"$splice_site"_splice_site.bed"
BAM="$OPTIMAL_ALIGNMENTS_DIR""$SAMPLE""_COMPASS_UGA_"$strand"_strand_sorted.bam"
OUT="$COMPASS_JUNCTIONS_DIR""$SAMPLE""_"$strand"_strand_"$splice_site"_splice_site_coverage.bed"
samtools depth -b $BED $BAM > $OUT 
done
done

python add_unspliced_read_counts_to_junctions.py $COMPASS_JUNCTIONS_DIR $SAMPLE
