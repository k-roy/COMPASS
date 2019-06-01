#!/bin/bash

DIR="/Volumes/LaCie/upf1_prp18_TruSeq_triplicates/"
FASTQ_DIR=$DIR"raw_fastq/"
SUBSAMPLED_DIR=$DIR"subsampled_fastq/"
mkdir $SUBSAMPLED_DIR
NUMBERED_READS_DIR=$DIR"numbered_reads_fastq/"
mkdir $NUMBERED_READS_DIR
TRIMMED_DIR=$DIR"trimmed_fastq/"
mkdir $TRIMMED_DIR
UNALIGNED_BAM_DIR=$DIR"unaligned_bam/"
mkdir $UNALIGNED_BAM_DIR
ALIGNMENTS_DIR=$DIR"alignments/"
mkdir $ALIGNMENTS_DIR

## STAR version 2.7.0d
YEAST_GENOME_DIR="/Volumes/MyPassportforMac_5TB/Dropbox/datasets/yeast_genome/"
GENOME_FASTA_FILE=$YEAST_GENOME_DIR"S288C_reference_genome_R64-2-1_20150113/S288C_reference_sequence_R64-2-1_20150113_reformatted_chromosome_names.fasta"
GENOME_DIRECTORY=$DIR"STAR_indexed_R64-2-1_20150113_annotated_genome_100_bp_SJDB"
mkdir $GENOME_DIRECTORY
GTF_FILE=$YEAST_GENOME_DIR"S288C_reference_genome_R64-2-1_20150113/saccharomyces_cerevisiae_R64-2-1_20150113_no_seq.gtf"
## generate genome index for 100-bp paired end reads - uncomment below code if needed - only needs to be done once
STAR --runThreadN 8 --runMode genomeGenerate --genomeDir $GENOME_DIRECTORY --genomeFastaFiles $GENOME_FASTA_FILE --sjdbGTFfile $GTF_FILE --sjdbOverhang 99 --sjdbGTFfeatureExon CDS


#SAMPLES="upf1 upf1_prp18";
#for sample in $SAMPLES;
#do for i in {1..3}; 

SAMPLES="upf1";
for sample in $SAMPLES;
do for i in {1..1}; 
do prefix=$sample"_rep"$i
R1=$prefix"_R1.fastq.gz"; 
R2=$prefix"_R2.fastq.gz";
echo $R1 $R2


prefix="subsampled_"$sample"_rep"$i


## BBMap reformat.sh version 38.32
reformat.sh reads=10000 in1=$FASTQ_DIR$R1 in2=$FASTQ_DIR$R2 \
out1=$SUBSAMPLED_DIR$prefix"_R1.fastq.gz" out2=$SUBSAMPLED_DIR$prefix"_R2.fastq.gz" \
overwrite=true

zcat < $SUBSAMPLED_DIR$prefix"_R1.fastq.gz" | awk '{print (NR%4 == 1) ? "@" ++i "_R1": $0}' | gzip -c > $NUMBERED_READS_DIR$prefix"_R1.fastq.gz"
zcat < $SUBSAMPLED_DIR$prefix"_R2.fastq.gz" | awk '{print (NR%4 == 1) ? "@" ++i "_R2": $0}' | gzip -c > $NUMBERED_READS_DIR$prefix"_R2.fastq.gz"
done; done


## cutadapt version 1.18

in_R1=$NUMBERED_READS_DIR$prefix"_R1.fastq.gz"
in_R2=$NUMBERED_READS_DIR$prefix"_R2.fastq.gz"
out_R1=$TRIMMED_DIR$prefix"_trimmed_R1.fastq"
out_R2=$TRIMMED_DIR$prefix"_trimmed_R2.fastq"
cutadapt --overlap 2 -j 0 -q 20,20 -g "T{100}" -a AGATCGGAAGAGC -A AGATCGGAAGAGC -A "A{100}" \
 -n 2 --trim-n  --minimum-length 50 --max-n 4 -o $out_R1  -p $out_R2 $in_R1 $in_R2

## --overlap 2 specifies that a partial match of 2 bases of the adapter prefix will allow trimming
## -j 0 specifies auto-detection for optimal number of cores to use
# -q 20,20 removes bases with lower than q score of 20 at 5p and 3p ends of reads
## --trim-n trims N's from the end, --max-n 4 discards reads with more than 4 N's,
## g means trim sequence at beginning of read1
## a means trim sequence at end of read1
## -A means trim sequence at end of read2
## -n 2 means trim first adapter specified by -A, then second adapter specified by -A in that order

## STAR version 2.7.0d
STAR_DIR=$ALIGNMENTS_DIR"STAR_default_annotated/"
mkdir $STAR_DIR
in_R1=$out_R1
in_R2=$out_R2
out=$STAR_DIR$prefix
STAR --runThreadN 8 --genomeDir $GENOME_DIRECTORY --sjdbOverhang 99 \
--readFilesIn $in_R1 $in_R2 --outFileNamePrefix $out \
--alignEndsType EndToEnd --outSAMattributes NH HI NM MD AS nM jM jI XS
java -jar /Applications/jvarkit/dist/samfixcigar.jar -r $GENOME_FASTA_FILE $out"Aligned.out.sam" -o $out"_reformatted_cigar.sam"
## samtools 1.9
samtools sort -n -o $out"_reformatted_cigar_name_sorted.bam" $out"_reformatted_cigar.sam" 

# 
# ## STAR version 2.7.0d
STAR_DIR=$ALIGNMENTS_DIR"STAR_noncanonical_annotated/"
mkdir $STAR_DIR
out=$STAR_DIR$prefix
STAR --runThreadN 8 --genomeDir $GENOME_DIRECTORY --sjdbOverhang 99 \
--readFilesIn $in_R1 $in_R2 --outFileNamePrefix $out \
--alignEndsType EndToEnd --scoreGapNoncan 0 --scoreGapGCAG 0 --scoreGapATAC 0 \
--outSAMattributes NH HI NM MD AS nM jM jI XS
java -jar /Applications/jvarkit/dist/samfixcigar.jar -r $GENOME_FASTA_FILE $out"Aligned.out.sam" -o $out"_reformatted_cigar.sam"
samtools sort -n -o $out"_reformatted_cigar_name_sorted.bam" $out"_reformatted_cigar.sam" 

##BBMap bbmap.sh version 38.32
BBMAP_DIR=$ALIGNMENTS_DIR"bbmap/"
mkdir $BBMAP_DIR
out=$BBMAP_DIR$prefix"_bbmap.sam"
bbmap.sh in1=$in_R1 in2=$in_R2 ref=$GENOME_FASTA_FILE out=$out \
intronlen=20 nhtag=t mdtag=t nhtag=t xmtag=t  amtag=t  nmtag=t tipsearch=1000 pairlen=10000
samtools sort -n -o $BBMAP_DIR$prefix"_bbmap_name_sorted.bam" $out 

done; done


