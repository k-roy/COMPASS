#!/bin/bash

# genome fasta downloaded from SGD (chromosome names reformatted to standard chrI, chrII etc manually)

# conda install -c bioconda gffread
# gffread-0.12.1
# also change CDS and noncoding_exon in feature field to exon

# conda install -c bioconda cutadapt
# cutadapt-1.18

# conda install -c bioconda bbmap
# BBMap version 38.18

# conda install -c bioconda star
# star-2.7.6a

# conda install -c bioconda hisat2
# hisat2-2.2.1

# git clone "https://github.com/lindenb/jvarkit.git"
# cd jvarkit
# ./gradlew samfixcigar

## revise script
cd /mnt/mindrinos/kevinroy/projects/COMPASS/
echo '' > COMPASS_process_reads_and_align_server_version.sh
nano COMPASS_process_reads_and_align_server_version.sh

# batch processing Aslanzadeh et al.
ACCESSION="PRJNA387451"
READ_LENGTH=150
for i in {5582776..5582781};
do SAMPLE="SRR"$i
echo $SAMPLE;
echo '' > "$SAMPLE.log"
nohup sh COMPASS_process_reads_and_align_server_version.sh -A "$ACCESSION" -S "$SAMPLE" -L "$READ_LENGTH" > "$SAMPLE.log" &
done

## batch processing Talkish et al.
ACCESSION="PRJNA354419"
READ_LENGTH=150
for i in {5041706..5041709};
do SAMPLE="SRR"$i
echo $SAMPLE;
echo '' > "$SAMPLE.log"
nohup sh COMPASS_process_reads_and_align_server_version.sh -A "$ACCESSION" -S "$SAMPLE" -L "$READ_LENGTH" > "$SAMPLE.log" &
done

# batch processing prp18/Roy et al.
DIR="/mnt/mindrinos/kevinroy/projects/COMPASS/"
ACCESSION="PRJNA544962"
LOG_DIR=$DIR$ACCESSION"_log_files/"
mkdir $LOG_DIR
READ_LENGTH=150
for i in {9130287..9130292};
do SAMPLE="SRR"$i
log_out=$LOG_DIR$SAMPLE"_COMPASS_process_reads_and_align.log"
echo $SAMPLE;
echo '' > $log_out
nohup sh COMPASS_process_reads_and_align_server_version.sh -A "$ACCESSION" -S "$SAMPLE" -L "$READ_LENGTH" > $log_out &
done

## single sample processing
# ACCESSION="PRJNA387451"
# READ_LENGTH=150
# SAMPLE="SRR5582779_subsampled"
# echo '' > "$SAMPLE.log"
# nohup sh COMPASS_process_reads_and_align_server_version.sh -A "$ACCESSION" -S "$SAMPLE" -L "$READ_LENGTH" > "$SAMPLE.log" &

##################################################################################################################################

while getopts A:S:L: flag
do
    case "${flag}" in
		A) ACCESSION=${OPTARG};;
        S) SAMPLE=${OPTARG};;
		L) READ_LENGTH=${OPTARG};;
    esac
done
echo "ACCESSION: $ACCESSION";
echo "SAMPLE: $SAMPLE";
echo "READ_LENGTH: $READ_LENGTH";

NUM_THREADS=24

COMPASS_DIR="/mnt/mindrinos/kevinroy/projects/COMPASS/"
OUT_DIR=$COMPASS_DIR"processed_data/"
mkdir $OUT_DIR
## ASSIGN AND CREATE SUBDIRECTORIES
NUMBERED_READS_DIR=$OUT_DIR"numbered_reads_fastq/"
mkdir $NUMBERED_READS_DIR
TRIMMED_DIR=$OUT_DIR"trimmed_fastq/"
mkdir $TRIMMED_DIR
ALIGNMENTS_DIR=$OUT_DIR"alignments/"
mkdir $ALIGNMENTS_DIR

GENOME_DIR=$COMPASS_DIR"S288C_reference_genome_R64-2-1_20150113/"
FASTA=$GENOME_DIR"S288C_reference_sequence_R64-2-1_20150113_reformatted_chromosome_names.fasta"
GTF=$GENOME_DIR"saccharomyces_cerevisiae_R64-2-1_20150113_exon_features.gtf"

STAR_GENOME_DIR=$GENOME_DIR"STAR_annotated_"$READ_LENGTH"_bp_SJDB_index/"
STAR_OVERHANG=$(expr $READ_LENGTH - 1)

HISAT2_GENOME_DIR=$GENOME_DIR"HISAT2_annotated_index/"
GENOME_NAME="Scer_R64_2_1"

SPLICE_SITES=$HISAT2_GENOME_DIR"splicesites.txt"
EXONS=$HISAT2_GENOME_DIR"exons.txt"

SAMFIXCIGAR="/mnt/mindrinos/kevinroy/projects/COMPASS/jvarkit/dist/samfixcigar.jar"

RAW_FASTQ_DIR="/mnt/mindrinos/kevinroy/projects/COMPASS/raw_data/"$ACCESSION"/"$SAMPLE"/"

raw_R1=$RAW_FASTQ_DIR$SAMPLE"_1.fastq.gz"
raw_R2=$RAW_FASTQ_DIR$SAMPLE"_2.fastq.gz"
trimmed_R1=$TRIMMED_DIR$SAMPLE"_trimmed_R1.fastq"
trimmed_R2=$TRIMMED_DIR$SAMPLE"_trimmed_R2.fastq"
numbered_R1=$NUMBERED_READS_DIR$SAMPLE"_numbered_R1.fastq"
numbered_R2=$NUMBERED_READS_DIR$SAMPLE"_numbered_R2.fastq"
## TRIM TRUSEQ ADAPTER AND POLY-A TAILS IN ONE STEP WITH CUTADAPT
## cutadapt version 1.18

cutadapt --overlap 2 -j 1 -q 20,20 -g "T{100}" -a AGATCGGAAGAGC -A AGATCGGAAGAGC -A "A{100}" \
-n 2 --trim-n --minimum-length 50 --max-n 4 -o $trimmed_R1 -p $trimmed_R2 $raw_R1 $raw_R2

# ## --overlap 2 specifies that a partial match of 2 bases of the adapter SAMPLE will allow trimming
## -j 0 specifies auto-detection for optimal number of cores to use
# -q 20,20 removes bases with lower than q score of 20 at 5p and 3p ends of reads
## --trim-n trims N's from the end, --max-n 4 discards reads with more than 4 N's,
## g means trim sequence at beginning of read1
## a means trim sequence at end of read1
## -A means trim sequence at end of read2
## -n 2 means trim first adapter specified by -A, then second adapter specified by -A in that order

cat < $trimmed_R1 | awk '{print (NR%4 == 1) ? "@" ++i "_R1": $0}' > $numbered_R1
cat < $trimmed_R2 | awk '{print (NR%4 == 1) ? "@" ++i "_R2": $0}' > $numbered_R2

cd $ALIGNMENTS_DIR
## MAP READS WITH BBMAP
##BBMap bbmap.sh version 38.32
BBMAP_DIR=$ALIGNMENTS_DIR"bbmap/"
mkdir $BBMAP_DIR
out=$BBMAP_DIR$SAMPLE
bbmap.sh in1=$numbered_R1 in2=$numbered_R2 ref=$FASTA out=$out".bam" \
nhtag=t mdtag=t nhtag=t xmtag=t  amtag=t nmtag=t tipsearch=2000 pairlen=10000

## MAP READS WITH STAR USING DEFAULT SETTINGS
## STAR version 2.7.0d
STAR_DIR=$ALIGNMENTS_DIR"STAR_default/"
mkdir $STAR_DIR
out=$STAR_DIR$SAMPLE
STAR --runThreadN $NUM_THREADS --genomeDir $STAR_GENOME_DIR --sjdbOverhang $STAR_OVERHANG \
--readFilesIn $numbered_R1 $numbered_R2 --outFileNamePrefix $out \
--alignEndsType EndToEnd --outSAMattributes NH HI NM MD AS nM jM jI XS
samtools view -bS -o $out".bam" $out"Aligned.out.sam" 
rm $out"Aligned.out.sam" 

## MAP READS WITH STAR USING MODIFIED SETTINGS TO ALLOW FOR NON-CANONICAL JUNCTIONS
## STAR version 2.7.0d

STAR_DIR=$ALIGNMENTS_DIR"STAR_noncanonical/"
mkdir $STAR_DIR
out=$STAR_DIR$SAMPLE
STAR --runThreadN $NUM_THREADS --genomeDir $STAR_GENOME_DIR --sjdbOverhang $STAR_OVERHANG \
--readFilesIn $numbered_R1 $numbered_R2 --outFileNamePrefix $out \
--alignEndsType EndToEnd --outSAMattributes NH HI NM MD AS nM jM jI XS \
--scoreGapNoncan 0 --scoreGapGCAG 0 --scoreGapATAC 0
samtools view -bS -o $out".bam" $out"Aligned.out.sam" 
rm $out"Aligned.out.sam" 

HISAT2_DIR=$ALIGNMENTS_DIR"HISAT2_default/"
mkdir $HISAT2_DIR
## MAP READS WITH HISAT2
out=$HISAT2_DIR$SAMPLE
hisat2 \
	--known-splicesite-infile $SPLICE_SITES \
	--no-softclip \
	--threads $NUM_THREADS \
	--time \
	--reorder \
	-x $HISAT2_GENOME_DIR$GENOME_NAME \
	-1 $numbered_R1 \
	-2 $numbered_R2 \
	-S $out".sam" \
	--min-intronlen 20 \
	--max-intronlen 10000 \
	--rna-strandness RF \
	--novel-splicesite-outfile $HISAT2_DIR$SAMPLE"_HISAT2_splice_junctions.txt" \
	--summary-file $HISAT2_DIR$SAMPLE"_HISAT2_summary.txt" \
	--new-summary
samtools view -bS -o $out".bam" $out".sam" 
rm $out".sam"

HISAT2_DIR=$ALIGNMENTS_DIR"HISAT2_noncanonical/"
mkdir $HISAT2_DIR
## MAP READS WITH HISAT2
out=$HISAT2_DIR$SAMPLE
hisat2 \
	--known-splicesite-infile $SPLICE_SITES \
	--no-softclip \
	--threads $NUM_THREADS \
	--time \
	--reorder \
	-x $HISAT2_GENOME_DIR$GENOME_NAME \
	-1 $numbered_R1 \
	-2 $numbered_R2 \
	-S $out".sam" \
	--min-intronlen 20 \
	--max-intronlen 10000 \
	--pen-noncansplice 0 \
	--rna-strandness RF \
	--novel-splicesite-outfile $HISAT2_DIR$SAMPLE"_HISAT2_splice_junctions.txt" \
	--summary-file $HISAT2_DIR$SAMPLE"_HISAT2_summary.txt" \
	--new-summary
samtools view -bS -o $out".bam" $out".sam" 
rm $out".sam"

for RUN_MODE in bbmap/ STAR_default/ STAR_noncanonical/ HISAT2_default/ HISAT2_noncanonical/;
do echo $RUN_MODE;
RUN_MODE_DIR=$ALIGNMENTS_DIR$RUN_MODE;
out=$RUN_MODE_DIR$SAMPLE;
# samtools sort -o $out"_coord_sorted.bam" $out".bam"
# rm $out".bam"
java -jar $SAMFIXCIGAR -r $FASTA -o $out"_reformatted_cigar.bam" --samoutputformat BAM $out"_coord_sorted.bam"  
samtools sort -n -o $out"_name_sorted.bam" -@ $NUM_THREADS $out"_reformatted_cigar.bam"
done

