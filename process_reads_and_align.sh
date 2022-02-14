#!/bin/bash

## STAR (default mode and non-canonical mode), HISAT2 (default mode and non-canonical mode), and GSNAP are guided by splice site annotations.
## BBMap and MAGIC-BLAST are splice-site agnostic.
## Minimap and graphmap are recently developed, popular aligners, 
## but are geared towards long, error prone nanopore reads and not splice junction discovery.

while getopts C:R:S:F:G:A:M:N:P:L:Z: flag;
do
    case "${flag}" in
		C) COMPASS_DIR=${OPTARG};;
		R) REFERENCE_DIR=${OPTARG};;
        S) SAMPLE=${OPTARG};;
		F) FASTA=${OPTARG};;
		G) GTF=${OPTARG};;
        A) SAMFIXCIGAR=${OPTARG};;
		M) MIN_INTRON_LENGTH=${OPTARG};;
		N) MAX_INTRON_LENGTH=${OPTARG};;
        P) NUM_THREADS=${OPTARG};;
		L) READ_LENGTH=${OPTARG};;
        Z) READS_TO_PROCESS=${OPTARG};;
    esac
done

echo "COMPASS_DIR: $COMPASS_DIR";
echo "REFERENCE_DIR: $REFERENCE_DIR";
echo "SAMPLE: $SAMPLE";
echo "FASTA: $FASTA";
echo "GTF: $GTF";
echo "MIN_INTRON_LENGTH: $MIN_INTRON_LENGTH";
echo "MAX_INTRON_LENGTH: $MAX_INTRON_LENGTH";
echo "NUM_THREADS: $NUM_THREADS";
echo "READ_LENGTH: $READ_LENGTH";
echo "READS_TO_PROCESS: $READS_TO_PROCESS";

##################### VARIABLE AND PATH DEFINITIONS BEGIN #####################
BBMAP_GENOME_DIR=$REFERENCE_DIR"bbmap/"

HISAT2_GENOME_DIR=$REFERENCE_DIR"HISAT2_annotated_index"
mkdir $HISAT2_GENOME_DIR
GENOME_VERSION=$(basename "$FASTA" .fasta)
HISAT2_INDEX=$HISAT2_GENOME_DIR"/"$GENOME_VERSION
SPLICE_SITES=$HISAT2_GENOME_DIR"/"$GENOME_VERSION"_splice_sites.txt"
EXONS=$HISAT2_GENOME_DIR"/"$GENOME_VERSION"_exons.txt"

STAR_GENOME_DIR=$REFERENCE_DIR"STAR_annotated_"$READ_LENGTH"_bp_SJDB_index"

BLAST_GENOME_DIR=$REFERENCE_DIR"BLAST/"
mkdir $BLAST_GENOME_DIR
BLAST_INDEX=$BLAST_GENOME_DIR$GENOME_VERSION

GSNAP_GENOME_DIR=$REFERENCE_DIR"GSNAP/"
mkdir $GSNAP_GENOME_DIR
###################### VARIABLE AND PATH DEFINITIONS END ######################

RAW_FASTQ_DIR=$COMPASS_DIR"fastq/"
NUMBERED_READS_DIR=$COMPASS_DIR"numbered_reads_fastq/"
TRIMMED_DIR=$COMPASS_DIR"trimmed_fastq/"
SEPARATE_ALIGNMENTS_DIR=$COMPASS_DIR"separate_alignments/"

raw_R1=$RAW_FASTQ_DIR$SAMPLE"_R1.fastq.gz"
raw_R2=$RAW_FASTQ_DIR$SAMPLE"_R2.fastq.gz"

if [[ $READS_TO_PROCESS -lt 1 ]];
then
echo "processing all reads for "$SAMPLE
else
echo "subsampling "$READS_TO_PROCESS" reads for "$SAMPLE
SAMPLE=$SAMPLE"_subsampled"
# For initial testing of the COMPASS workflow, 
# it is recommended to process only a single sample 
# and to subsample the reads in the range from 100,000 to 1 million 
subsampled_R1=$RAW_FASTQ_DIR$SAMPLE"_1.fastq.gz"
subsampled_R2=$RAW_FASTQ_DIR$SAMPLE"_2.fastq.gz"
# reformat.sh in1=$raw_R1 in2=$raw_R2 out1=$subsampled_R1 out2=$subsampled_R2 samplereadstarget=$READS_TO_PROCESS
raw_R1=$subsampled_R1
raw_R2=$subsampled_R2
fi

trimmed_R1=$TRIMMED_DIR$SAMPLE"_trimmed_R1.fastq"
trimmed_R2=$TRIMMED_DIR$SAMPLE"_trimmed_R2.fastq"
numbered_R1=$NUMBERED_READS_DIR$SAMPLE"_numbered_R1.fastq"
numbered_R2=$NUMBERED_READS_DIR$SAMPLE"_numbered_R2.fastq"

############################### CUTADAPT BEGIN ################################
## TRIM TRUSEQ ADAPTER AND POLY-A TAILS FROM 3' ENDS 
## AND LOW QUALITY BASES FROM BOTH ENDS IN ONE STEP WITH CUTADAPT
## cutadapt version 1.18

cutadapt --overlap 2 --cores $NUM_THREADS -q 20,20 -g "T{100}" -a AGATCGGAAGAGC -A AGATCGGAAGAGC -A "A{100}" \
-n 2 --trim-n --minimum-length 50 --max-n 4 -o $trimmed_R1 -p $trimmed_R2 $raw_R1 $raw_R2

# ## --overlap 2 specifies that a partial match of 2 bases of the adapter SAMPLE will allow trimming
## -j 0 specifies auto-detection for optimal number of cores to use
# -q 20,20 removes bases with lower than q score of 20 at 5p and 3p ends of reads
## --trim-n trims N's from the end, --max-n 4 discards reads with more than 4 N's,
## g means trim sequence at beginning of read1
## a means trim sequence at end of read1
## -A means trim sequence at end of read2
## -n 2 means trim first adapter specified by -A, then second adapter specified by -A in that order
################################ CUTADAPT END #################################

########################### RENUMBER READS BEGIN ##############################
cat < $trimmed_R1 | awk '{print (NR%4 == 1) ? "@" ++i "_R1": $0}' > $numbered_R1
cat < $trimmed_R2 | awk '{print (NR%4 == 1) ? "@" ++i "_R2": $0}' > $numbered_R2
rm $trimmed_R1
rm $trimmed_R2
############################ RENUMBER READS END ###############################

cd $SEPARATE_ALIGNMENTS_DIR

############################# BBMAP BEGIN #####################################
# ## MAP READS WITH BBMAP
# ##BBMap bbmap.sh version 38.32
BBMAP_DIR=$SEPARATE_ALIGNMENTS_DIR"bbmap/"
mkdir $BBMAP_DIR
out=$BBMAP_DIR$SAMPLE
# # BBMAP will automatically generate a reference index
bbmap.sh in1=$numbered_R1 in2=$numbered_R2 ref=$FASTA out=$out".bam" \
path=$BBMAP_GENOME_DIR unpigz=t pigz=t \
nhtag=t mdtag=t nhtag=t xmtag=t amtag=t nmtag=t \
maxindel=$MAX_INTRON_LENGTH pairlen=$MAX_INTRON_LENGTH intronlen=$MIN_INTRON_LENGTH

# from the BBMap manual:
# To map vertebrate RNA-seq reads to a genome:
# bbmap.sh in=reads.fq out=mapped.sam maxindel=200k ambig=random intronlen=20 xstag=us

######################################### BBMAP END ################################################

######################################### HISAT2 BEGIN ##############################################
HISAT2_DIR=$SEPARATE_ALIGNMENTS_DIR"HISAT2_default/"
mkdir $HISAT2_DIR
## CHECK IF GENOME INDEX ALREADY EXISTS FOR HISAT2
## IF NOT, GENOME INDEX FOR HISAT2
# hisat2-2.2.1

DIR_TO_CHECK=$HISAT2_GENOME_DIR
if [ -n "$(find "$DIR_TO_CHECK" -maxdepth 0 -type d -empty 2>/dev/null)" ]; then
    echo "HISAT2_GENOME_DIR is empty"
    samtools faidx $FASTA
    picard CreateSequenceDictionary -R $FASTA
    hisat2_extract_splice_sites.py $GTF > $SPLICE_SITES
    hisat2_extract_exons.py $GTF > $EXONS
    hisat2-build -p $NUM_THREADS --ss $SPLICE_SITES --exon $EXONS $FASTA $HISAT2_INDEX
else
    echo "The HISAT2 genome database is already built."
fi

# MAP READS WITH HISAT2
out=$HISAT2_DIR$SAMPLE
hisat2 \
	--known-splicesite-infile $SPLICE_SITES \
	--no-softclip \
	--threads $NUM_THREADS \
	--time \
	--reorder \
	-x $HISAT2_INDEX \
	-1 $numbered_R1 \
	-2 $numbered_R2 \
	-S $out".sam" \
	--min-intronlen $MIN_INTRON_LENGTH \
	--max-intronlen $MAX_INTRON_LENGTH \
	--rna-strandness RF \
	--novel-splicesite-outfile $HISAT2_DIR$SAMPLE"_HISAT2_splice_junctions.txt" \
	--summary-file $HISAT2_DIR$SAMPLE"_HISAT2_summary.txt" \
	--new-summary
samtools view -bS -@ $NUM_THREADS -o $out".bam" $out".sam" 
rm $out".sam"

HISAT2_DIR=$SEPARATE_ALIGNMENTS_DIR"HISAT2_noncanonical/"
mkdir $HISAT2_DIR
## MAP READS WITH HISAT2
out=$HISAT2_DIR$SAMPLE
hisat2 \
	--known-splicesite-infile $SPLICE_SITES \
	--no-softclip \
	--threads $NUM_THREADS \
	--time \
	--reorder \
	-x $HISAT2_INDEX \
	-1 $numbered_R1 \
	-2 $numbered_R2 \
	-S $out".sam" \
	--min-intronlen $MIN_INTRON_LENGTH \
	--max-intronlen $MAX_INTRON_LENGTH \
	--pen-noncansplice 0 \
	--rna-strandness RF \
	--novel-splicesite-outfile $HISAT2_DIR$SAMPLE"_HISAT2_splice_junctions.txt" \
	--summary-file $HISAT2_DIR$SAMPLE"_HISAT2_summary.txt" \
	--new-summary
samtools view -bS -@ $NUM_THREADS -o $out".bam" $out".sam" 
rm $out".sam"
######################################### HISAT2 END ##############################################

######################################### STAR BEGIN ##############################################

## CHECK IF GENOME INDEX ALREADY EXISTS FOR STAR
## IF NOT, GENERATE GENOME INDEX FOR STAR
# star-2.7.6a
## generate genome index for 100-bp paired end reads - uncomment below code if needed - only needs to be done once
STAR_OVERHANG=$(expr $READ_LENGTH - 1)
mkdir $STAR_GENOME_DIR

DIR_TO_CHECK=$STAR_GENOME_DIR
if [ -n "$(find "$DIR_TO_CHECK" -maxdepth 0 -type d -empty 2>/dev/null)" ]; then
    echo "STAR_GENOME_DIR is empty"
    # STAR needs GTF: https://github.com/alexdobin/STAR/issues/387
    STAR --runThreadN $NUM_THREADS --runMode genomeGenerate --genomeDir $STAR_GENOME_DIR \
    --genomeFastaFiles $FASTA --sjdbGTFfile $GTF --sjdbOverhang $STAR_OVERHANG --sjdbGTFfeatureExon exon
else
    echo "The STAR genome database is already built."
fi

## MAP READS WITH STAR USING DEFAULT SETTINGS
## STAR version 2.7.9a
STAR_DIR=$SEPARATE_ALIGNMENTS_DIR"STAR_default/"
mkdir $STAR_DIR
out=$STAR_DIR$SAMPLE
STAR --runThreadN $NUM_THREADS --genomeDir $STAR_GENOME_DIR --sjdbOverhang $STAR_OVERHANG \
--readFilesIn $numbered_R1 $numbered_R2 --outFileNamePrefix $out \
--alignEndsType EndToEnd --outSAMattributes NH HI NM MD AS nM jM jI XS
samtools view -bS -@ $NUM_THREADS -o $out".bam" $out"Aligned.out.sam" 
rm $out"Aligned.out.sam" 

## MAP READS WITH STAR USING MODIFIED SETTINGS TO ALLOW FOR NON-CANONICAL JUNCTIONS
## STAR version 2.7.9a
STAR_DIR=$SEPARATE_ALIGNMENTS_DIR"STAR_noncanonical/"
mkdir $STAR_DIR
out=$STAR_DIR$SAMPLE
STAR --runThreadN $NUM_THREADS --genomeDir $STAR_GENOME_DIR --sjdbOverhang $STAR_OVERHANG \
--readFilesIn $numbered_R1 $numbered_R2 --outFileNamePrefix $out \
--alignEndsType EndToEnd --outSAMattributes NH HI NM MD AS nM jM jI XS \
--scoreGapNoncan 0
samtools view -bS -@ $NUM_THREADS -o $out".bam" $out"Aligned.out.sam" 
rm $out"Aligned.out.sam" 
######################################### STAR END ##############################################

######################################### MAGIC-BLAST START ##############################################
MAGIC_BLAST_DIR=$SEPARATE_ALIGNMENTS_DIR"MAGIC_BLAST/"
mkdir $MAGIC_BLAST_DIR
## CHECK IF BLAST DATABASE ALREADY EXISTS FOR MAGIC-BLAST
## IF NOT, GENERATE BLAST DATABASE FOR MAGIC-BLAST

DIR_TO_CHECK=$BLAST_GENOME_DIR
if [ -n "$(find "$DIR_TO_CHECK" -maxdepth 0 -type d -empty 2>/dev/null)" ]; then
    echo "BLAST_GENOME_DIR is empty"
    makeblastdb -in $FASTA -dbtype nucl -parse_seqids -out $BLAST_INDEX -title $GENOME_VERSION
else
    echo "The BLAST genome database is already built."
fi

out_prefix=$MAGIC_BLAST_DIR$SAMPLE
# threads reduced for MAGIC BLAST to avoid crashing due to excessive memory usage
magicblast -query $numbered_R1 -query_mate $numbered_R2 -db $BLAST_INDEX \
-md_tag -fr  -no_query_id_trim -infmt fastq -num_threads 12 \
-max_db_word_count 10 -out $out_prefix".sam"
samtools view -bS -@ $NUM_THREADS -o $out_prefix".bam" $out_prefix".sam" 
rm $out_prefix".sam" 
######################################### MAGIC-BLAST END ##############################################

################################## GSNAP BEGIN ################################
GSNAP_DIR=$SEPARATE_ALIGNMENTS_DIR"GSNAP/"
mkdir $GSNAP_DIR
out_prefix=$GSNAP_DIR$SAMPLE

DIR_TO_CHECK=$GSNAP_GENOME_DIR
if [ -n "$(find "$DIR_TO_CHECK" -maxdepth 0 -type d -empty 2>/dev/null)" ]; then
    echo "GSNAP_GENOME_DIR is empty"
    gmap_build -D $GSNAP_GENOME_DIR -d $GENOME_VERSION $FASTA
    cat $GTF | gtf_splicesites | awk '$4>9' > $GSNAP_GENOME_DIR$GENOME_VERSION".splicesites"
    # There are issues with how gtf_splicesites parses the S.cerevisiae GTF.
    # It assumes even single exon genes have splice sites on either side of the exon, resulting in negative length introns.
	# It also processes annotations like +1 translational frameshifts as exons, resuling in 1 nt long introns.
    # Filtering the intron size for at least 10 nt seems to solve the issue .
    cat $GSNAP_GENOME_DIR$GENOME_VERSION".splicesites" | iit_store -o $GSNAP_GENOME_DIR$GENOME_VERSION
    cp $GSNAP_GENOME_DIR$GENOME_VERSION".iit" $GSNAP_GENOME_DIR$GENOME_VERSION"/"$GENOME_VERSION".maps/"
else
    echo "The GMAP/GSNAP genome database is already built."
fi

gsnap -D $GSNAP_GENOME_DIR -d $GENOME_VERSION --use-splicing=$GENOME_VERSION \
$numbered_R1 $numbered_R2 --output-file $out_prefix".sam" \
--nthreads=$NUM_THREADS \
--ambig-splice-noclip --novelsplicing=1 \
--add-paired-nomappers --sam-extended-cigar --format=sam
samtools view -bS -@ $NUM_THREADS -o $out_prefix".bam" $out_prefix".sam"
rm $out_prefix".sam"
################################## GSNAP END ################################
# rm $numbered_R1
# rm $numbered_R2
############################### REFORMAT ALIGNMENTS BEGIN ##############################
## HISAT2, STAR, MAGIC-BLAST, and GSNAP map to sam files.
## Sam files are converted to bam and then deleted to save disk space.

for RUN_MODE in bbmap/ HISAT2_default/ HISAT2_noncanonical/ STAR_default/ STAR_noncanonical/ MAGIC_BLAST/ GSNAP/;  
do echo $RUN_MODE;
RUN_MODE_DIR=$SEPARATE_ALIGNMENTS_DIR$RUN_MODE;
out=$RUN_MODE_DIR$SAMPLE;
## sort bam by coordinate, as this is required for samfixcigar to run properly
samtools sort -@ $NUM_THREADS -o $out"_coord_sorted.bam" $out".bam"
rm $out".bam"
## samfixcigar produces consistent representation of mismatches as X and matches as M in SAM format 1.4
java -jar $SAMFIXCIGAR --reference $FASTA --out $out"_reformatted_cigar.bam" --samoutputformat BAM $out"_coord_sorted.bam" 
# sort by read name, which is also read number in this case 
samtools sort -n -@ $NUM_THREADS -o $out"_name_sorted.bam" $out"_reformatted_cigar.bam"
rm $out"_coord_sorted.bam"
# The reformatted cigar bam can be useful for accessing alignments 
# for individual aligners for specific regions later on.
done
############################### REFORMAT ALIGNMENTS END ##############################
