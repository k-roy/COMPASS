
## batch processing Aslanzadeh et al.
#COMPASS_DIR="/mnt/mindrinos/kevinroy/projects/COMPASS/"
#cd $COMPASS_DIR
#LOG_DIR=$COMPASS_DIR"log/"
#ACCESSION="PRJNA387451"
#READ_LENGTH=150
#for i in {5582776..5582781};
#for i in {5582780..5582780};
#do SAMPLE="SRR"$i
#echo $SAMPLE;
#LOG_FILE=$LOG_DIR$SAMPLE"_COMPASS.log"
#echo '' > $LOG_FILE
#nohup sh COMPASS_compare_splice_junctions_from_multiple_aligner_BAM.sh -A "$ACCESSION" -S "$SAMPLE" -L "$READ_LENGTH" > $LOG_FILE &
#done

# batch processing prp18/Roy et al.
#COMPASS_DIR="/mnt/mindrinos/kevinroy/projects/COMPASS/"
#cd $COMPASS_DIR
#LOG_DIR=$COMPASS_DIR"log/"
#ACCESSION="PRJNA544962"
#READ_LENGTH=150
#for i in {9130287..9130292};
#do SAMPLE="SRR"$i
#echo $SAMPLE;
#LOG_FILE=$LOG_DIR$SAMPLE"_COMPASS.log"
#echo '' > $LOG_FILE
#nohup sh COMPASS_compare_splice_junctions_from_multiple_aligner_BAM.sh -A "$ACCESSION" -S "$SAMPLE" -L "$READ_LENGTH" > $LOG_FILE &
#done

## batch processing Aslanzadeh et al.  USJ only (will be much faster)
#COMPASS_DIR="/mnt/mindrinos/kevinroy/projects/COMPASS/"
#cd $COMPASS_DIR
#LOG_DIR=$COMPASS_DIR"log/"
#ACCESSION="PRJNA387451"
#READ_LENGTH=150
#for i in {5582776..5582781};
#do SAMPLE="SRR"$i
#echo $SAMPLE;
#LOG_FILE=$LOG_DIR$SAMPLE"_COMPASS_USJ_only.log"
#echo '' > $LOG_FILE
#nohup sh COMPASS_compare_splice_junctions_from_multiple_aligner_BAM.sh -A "$ACCESSION" -S "$SAMPLE" -L "$READ_LENGTH" > $LOG_FILE &
#done

## batch processing for MPE seq dataset
#COMPASS_DIR="/mnt/mindrinos/kevinroy/projects/COMPASS/"
#cd $COMPASS_DIR
#LOG_DIR=$COMPASS_DIR"log/"
#ACCESSION="PRJNA472800"
#READ_LENGTH=60
#for SAMPLE in SRR7208762 SRR7208766 SRR7208767 SRR7208768 SRR7208769 SRR7208771 SRR7208772 SRR7208773 SRR7208774 SRR7208775 SRR7208776 SRR7208777;
#do echo $SAMPLE;
#LOG_FILE=$LOG_DIR$SAMPLE"_COMPASS.log"
#echo '' > $LOG_FILE
#nohup sh COMPASS_compare_splice_junctions_from_multiple_aligner_BAM.sh -A "$ACCESSION" -S "$SAMPLE" -L "$READ_LENGTH" > $LOG_FILE &
#done

## batch processing Talkish et al.
#COMPASS_DIR="/mnt/mindrinos/kevinroy/projects/COMPASS/"
#cd $COMPASS_DIR
#LOG_DIR=$COMPASS_DIR"log/"
#ACCESSION="PRJNA354419"
#READ_LENGTH=150
#for i in {5041706..5041709};
#do SAMPLE="SRR"$i
#echo $SAMPLE;
#LOG_FILE=$LOG_DIR$SAMPLE"_COMPASS.log"
#echo '' > $LOG_FILE
#nohup sh COMPASS_compare_splice_junctions_from_multiple_aligner_BAM.sh -A "$ACCESSION" -S "$SAMPLE" -L "$READ_LENGTH" > $LOG_FILE &
#done

## single sample processing 
# SRR5582776_subsampled_name_sorted.bam
#COMPASS_DIR="/mnt/mindrinos/kevinroy/projects/COMPASS/"
#cd $COMPASS_DIR
#LOG_DIR=$COMPASS_DIR"log/"
#ACCESSION="PRJNA387451"
#READ_LENGTH=150
#SAMPLE="SRR5582779_subsampled"
#echo '' > $LOG_DIR$SAMPLE"_COMPASS.log"
#nohup sh COMPASS_compare_splice_junctions_from_multiple_aligner_BAM.sh -A "$ACCESSION" -S "$SAMPLE" -L "$READ_LENGTH" > $LOG_DIR$SAMPLE"_COMPASS.log" &

# check log files
#cd /mnt/mindrinos/kevinroy/projects/COMPASS/
#for i in {5582776..5582781};
#do SAMPLE="SRR"$i
#echo "$SAMPLE.log";
#tail -15 "$SAMPLE.log";
#done

#download splice junction files
#IN_DIR="/Users/kevinroy/Desktop/qian/mnt/mindrinos/kevinroy/projects/COMPASS/processed_data/alignments/COMPASS_integration/"
#OUT_DIR="/Volumes/SPxDrive/COMPASS/"
#rsync -azcP $IN_DIR"SRR5582776_COMPASS_splice_junctions.tsv" $OUT_DIR


## revise script
#cd /mnt/mindrinos/kevinroy/projects/COMPASS/
#echo '' > COMPASS_compare_splice_junctions_from_multiple_aligner_BAM.sh
#nano COMPASS_compare_splice_junctions_from_multiple_aligner_BAM.sh


########################################################################
#while getopts A:S:L: flag
#do
    #case "${flag}" in
		#A) ACCESSION=${OPTARG};;
        #S) SAMPLE=${OPTARG};;
		#L) READ_LENGTH=${OPTARG};;
    #esac
#done
#echo "ACCESSION: $ACCESSION";
#echo "SAMPLE: $SAMPLE";
#echo "READ_LENGTH: $READ_LENGTH";

COMPASS_DIR="/mnt/mindrinos/kevinroy/projects/COMPASS/"
OUT_DIR=$COMPASS_DIR"processed_data/"
ALIGNMENTS_DIR=$OUT_DIR"alignments/"
COMPASS_OUT_DIR=$ALIGNMENTS_DIR"COMPASS_integration/"
mkdir $COMPASS_OUT_DIR
GENOME_DIR='/u/project/guillom/kevinh97/COMPASS/S288C_reference_genome_R64-2-1_20150113/'
FASTA='/u/project/guillom/kevinh97/COMPASS/S288C_reference_genome_R64-2-1_20150113/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa'
GTF='/u/project/guillom/kevinh97/saccharomyces_cerevisiae_R64-1-1_20110208_no_chr.gtf'
INTRONS_FILE='/u/project/guillom/kevinh97/COMPASS/S288C_reference_genome_R64-2-1_20150113/saccharomyces_cerevisiae_R64-2-1_20150113_introns_no_chr.tsv'
# for RUN_MODE in bbmap/ STAR_default/ STAR_noncanonical/ HISAT2_default/ HISAT2_noncanonical/;
ALIGNERS_FILE=$OUT_DIR"aligners.txt"
# SAMPLE_SUFFIX="_COMPASS_USJ.bam"
SAMPLE_SUFFIX="_name_sorted.bam"
READS_TO_PROCESS=100000000000
cd $COMPASS_DIR
python -u COMPASS_compare_splice_junctions_from_multiple_aligner_BAM.py "$ALIGNMENTS_DIR" "$COMPASS_OUT_DIR" "$FASTA" "$INTRONS_FILE" "$ALIGNERS_FILE" "$SAMPLE" "$SAMPLE_SUFFIX" "$READS_TO_PROCESS"

