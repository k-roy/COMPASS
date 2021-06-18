COMPASS_DIR="/u/project/guillom/kevinh97/COMPASS/"
GENOME_DIR=$COMPASS_DIR"S288C_reference_genome_R64-2-1_20150113/"
READ_LENGTH=150
STAR_OVERHANG=$(expr $READ_LENGTH - 1)
FASTA=$GENOME_DIR"Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa"
GTF=$GENOME_DIR"Saccharomyces_cerevisiae.R64-1-1.104.gtf"
NUM_THREAD=8

samtools faidx $FASTA
java -jar /u/project/guillom/kevinh97/picard/build/libs/picard.jar CreateSequenceDictionary -R $FASTA

## GENERATE GENOME INDEX FOR STAR
# star-2.7.6a
## generate genome index for 100-bp paired end reads - uncomment below code if needed - only needs to be done once
STAR_GENOME_DIR=$GENOME_DIR"STAR_annotated_"$READ_LENGTH"_bp_SJDB_index/"
# STAR needs GTF: https://github.com/alexdobin/STAR/issues/387
STAR_OVERHANG=$(expr $READ_LENGTH - 1)
STAR --runThreadN $NUM_THREAD --runMode genomeGenerate --genomeSAindexNbases 10 --genomeDir $STAR_GENOME_DIR \
--genomeFastaFiles $FASTA --sjdbGTFfile $GTF --sjdbOverhang $STAR_OVERHANG --sjdbGTFfeatureExon exon

# GENERATE GENOME INDEX FOR HISAT2
# hisat2-2.2.1
HISAT2_GENOME_DIR=$GENOME_DIR"HISAT2_annotated_index/"
GENOME_NAME="Scer_R64_2_1"
SPLICE_SITES=$HISAT2_GENOME_DIR"splicesites.txt"
EXONS=$HISAT2_GENOME_DIR"exons.txt"
hisat2_extract_splice_sites.py $GTF > $SPLICE_SITES
hisat2_extract_exons.py $GTF > $EXONS
hisat2-build -p $NUM_THREAD --ss $SPLICE_SITES --exon $EXONS $FASTA $HISAT2_GENOME_DIR$GENOME_NAME