
## batch processing Aslanzadeh et al.
COMPASS_DIR="/mnt/mindrinos/kevinroy/projects/COMPASS/"
LOG_DIR=$COMPASS_DIR"log/"
ACCESSION="PRJNA387451"
READ_LENGTH=150
for i in {5582776..5582781};
do SAMPLE="SRR"$i
echo $SAMPLE;
LOG_FILE=$LOG_DIR$SAMPLE"_COMPASS_specific_USJ_extraction.log"
echo '' > $LOG_FILE
nohup sh extract_COMPASS_USJ_alignments_from_multiple_aligner_BAM.sh -A "$ACCESSION" -S "$SAMPLE" -L "$READ_LENGTH" > $LOG_FILE &
done

# check log files
cd /mnt/mindrinos/kevinroy/projects/COMPASS/log/
for i in {5582776..5582781};
do SAMPLE="SRR"$i
LOG_FILE=$LOG_DIR$SAMPLE"_COMPASS_specific_USJ_extraction.log"
echo $LOG_FILE;
tail -15 $LOG_FILE;
done

#download splice junction files
IN_DIR="/Users/kevinroy/Desktop/qian/mnt/mindrinos/kevinroy/projects/COMPASS/processed_data/alignments/COMPASS_integration_full_output/"
OUT_DIR="/Volumes/SPxDrive/COMPASS/"
rsync -azcP $IN_DIR"XXX" $OUT_DIR

## revise script
cd /mnt/mindrinos/kevinroy/projects/COMPASS/
echo '' > extract_COMPASS_USJ_alignments_for_specific_junction_from_multiple_aligner_BAM.sh
nano extract_COMPASS_USJ_alignments_for_specific_junction_from_multiple_aligner_BAM.sh


########################################################################
while getopts A:S:L: flag
do
    case "${flag}" in
		A) ACCESSION=${OPTARG};;
        S) SAMPLE_NAMES=${OPTARG};;
		L) READ_LENGTH=${OPTARG};;
    esac
done
echo "ACCESSION: $ACCESSION";
echo "SAMPLE: $SAMPLE";
echo "READ_LENGTH: $READ_LENGTH";

COMPASS_DIR="/mnt/mindrinos/kevinroy/projects/COMPASS/"
OUT_DIR=$COMPASS_DIR"processed_data/"
ALIGNMENTS_DIR=$OUT_DIR"alignments/"
COMPASS_OUT_DIR=$ALIGNMENTS_DIR"COMPASS_integration_full_output/"
mkdir $COMPASS_OUT_DIR
GENOME_DIR=$COMPASS_DIR"S288C_reference_genome_R64-2-1_20150113/"
FASTA=$GENOME_DIR"S288C_reference_sequence_R64-2-1_20150113_reformatted_chromosome_names.fasta"
GTF=$GENOME_DIR"saccharomyces_cerevisiae_R64-2-1_20150113_exon_features.gtf"
INTRONS_FILE=$GENOME_DIR"saccharomyces_cerevisiae_R64-2-1_20150113_introns.tsv"
ALIGNERS_FILE=$OUT_DIR"aligners.txt"
SAMPLE_NAMES="SRR5582776 SRR5582777 SRR5582778 SRR5582779 SRR5582780 SRR5582781"
COMBINED_SAMPLES_NAME="Aslanzadeh_combined"
SAMPLE_SUFFIX="_COMPASS_USJ_coord_sorted.bam"
# chrVIII    75728   75832   105 +      
QUERY_CHROM="chrVIII"
QUERY_START="75728"
QUERY_END="75832"
QUERY_STRAND="+"
READS_TO_PROCESS=100000000000
cd $COMPASS_DIR
python -u extract_COMPASS_USJ_alignments_for_specific_junction_from_multiple_aligner_BAM.py \
"$ALIGNMENTS_DIR" "$COMPASS_OUT_DIR" "$FASTA" "$INTRONS_FILE" "$ALIGNERS_FILE" "$SAMPLE_NAMES" \
 "$SAMPLE_SUFFIX" "$COMBINED_SAMPLES_NAME" "$QUERY_CHROM" "$QUERY_START" "$QUERY_END" "$QUERY_STRAND"

# query_chrom = 'chrXVI'
# query_start = 218647 - 1# 218647 
# query_end = 218724 -1  # 218724
# query_strand = '+'
# script, ALIGNMENTS_DIR, OUT_DIR, FASTA, INTRONS_FILE, ALIGNERS_FILE, SAMPLE_NAMES, SAMPLE_SUFFIX, COMBINED_SAMPLES_NAME, QUERY_CHROM, QUERY_START, QUERY_END, QUERY_STRAND = argv

########################################################################

COMBINED_SAMPLES_NAME="Aslanzadeh_combined"
QUERY_CHROM="chrVIII"
QUERY_START="75728"
QUERY_END="75832"
QUERY_STRAND="+"

ALIGNMENTS_DIR=~/Desktop/qian/mnt/mindrinos/kevinroy/projects/COMPASS/processed_data/alignments/
for RUN_MODE in bbmap STAR_default STAR_noncanonical HISAT2_default HISAT2_noncanonical;
do echo $RUN_MODE;
IN_DIR=$ALIGNMENTS_DIR$RUN_MODE"/"
OUT_DIR="/Volumes/SPxDrive/COMPASS/processed_data/alignments/"$RUN_MODE"/"
rsync -azvP $IN_DIR$RUN_MODE"_Aslanzadeh_combined_"$QUERY_CHROM"_"$QUERY_START"_"QUERY_END"_sorted.bam"* $OUT_DIR
done

# bbmap_Aslanzadeh_combined_chrXVI_218646_218723_sorted.bam