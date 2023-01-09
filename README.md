# COMPASS

This repository includes scripts for the Comparison of Multiple alignment Programs for Alternative Splice Site discovery (COMPASS).

Here are the steps needed to run COMPASS.

1) COMPASS_install_required_programs.sh
Creates the compass environment in conda.
Installs samfixcigar from jvarkit.

2) COMPASS.sh (run for each sample):
process_reads_and_align.sh
Read 3′ trimming for polyA tails and base calling quality scores (cutadapt)
Assigns consecutive numbers to each read in the fastq file (awk).
Alignment with multiple aligners: BBMap and STAR (both default and noncanonical splicing modes).
Genome indices are built for each aligner only if necessary (i.e. when running the first time).
Sam files are converted to .bam format (samtools view), sorted by mapped reference coordinates (samtools sort), cigars reformatted to SAM format 1.4 (samfixcigar), and then sorted by read numbers (samtools sort).
compare_splice_junctions_from_multiple_aligners.py
Comparison of Multiple Programs for Accurate Splice Site discovery (COMPASS) core program
Adjustment of ambiguous junctions to most likely splice sites based on species-specific splice signals
add_exonic_intronic_sequence.py
Adjustment of ambiguous junctions to most likely splice sites based on species-specific splice signals
create_splice_site_bed.py
add_unspliced_read_counts_to_junctions.py

3) COMPASS_combine_junction_tables_from_multiple_samples.py
Alternative splice junction calling
Junction quality filtering
