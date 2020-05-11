#!/bin/bash
#
#$ -cwd # running the job at the current working directory
#$ -j y # Stdout and stderr from this script will be merged (-j y) and written to job.sh.o<job_id>.
#$ -S /bin/bash
#$ -pe smp 4
#          ^ this is the number of CPU cores your program will use,
# or the amount of memory it needs / 4GB, whichever is larger, but
# not more than 16. E.g. 4 is for up to 4 CPUs and 16GB RAM.
# These limits aren't enforced, so it's up to you to avoid oversubscription.
#
# To get an e-mail when the job is done:
#$ -m e
#$ -M your@email.com
#

##################################
# Author: Zhiyi Sun
# Date: last updated on Nov, 19, 2018
# Description: 
#	This pipeline script is for analyzing the Illumina sequencing data of circle PAM assay (developed by Ryan Fuchs).
#	It takes pair-end fastq files as input, after merging the pair-end reads into full-length fragments and mapping them to the reference sequence (need to be prepared in advance and indexed using the bowtie2-build utility) 
#	it determines the cleavage pattern and PAM sequences and reports the results to multiple output files.
# In details, this scripts performs the following tasks: 			 
# 1, trim adapters and merge forward and reverse read* (seqprep)
# 	*NOTE: Before merging the reads, run fastqc to check data quality read length (seqprep works best if both R1 and R2 reads contains the adapter sequence (partial is fine) at the 3' end and there is a good overlap between R1 and R2).
# 2, remove short reads (<100bp)
# 3, map merged reads to the reference (bowtie2)
# 4, characterize cutting pattern and extract PAM sequence from good alignments 
#        Criteria for determine cutting pattern
#        1, Read must map to the correct reference (XL or XXL)
#        2, pirmary alignment only (FLAG !=256)
#        3, MAPQ > 20
#        4, alignment left most position <=25 (upstream of the PAM region)
#        5, if the sequence before PAM region match the reference and
#              --if read length = reference original lenghth (XL=140; XXL=120):  blunt end
#              --if read length < reference original length: 3'overhang, overhang length= reference length - read length       
#              --if read length > reference original length and the extra bases at the 3' end wrt the reference is same as the beginning of the read: 5' overhang, overhang length= read length - reference length
#       calls a Perl script: [characterize_cut_pattern_from_circle_PAM.pl]
#
# This pipeline will generate the following intermediate and final output files:
# 1, adapter trimmed R1 fastq file (output of step 1): *_trimmed.fq.gz
# 2, adapter trimmed R2 fastq file (output of step 1): *_trimmed.fq.gz
# 3, fastq file for merged reads (output of step 1): *_merged.fq.gz
# 4, fastq file for merged reads with desired read length, default is >=100bp (output of step 2): *_merged_100bp+.fq
# 5, summary of read length trimming result (output of step 2): *_merged_length_trimming_report.txt
# 6, bowtie2 output sam file (output of step 3): *_merged_100bp+_bt2.sam
# 7, result of cleavage pattern (blunt, 5overhang, 3overhang); distance to PAM and overhang size and PAM sequences of each good alignment (output of step 4) : *_merged_100bp+_cut_pattern_result.txt
# 8, summary of the final results (output of step 4): *_merged_100bp+_cut_pattern_result_summary.txt
# 9, count and fraction matrices of blunt cuts (count for each distance to PAM: 1 to 25): *_merged_100bp+_cut_pattern_result_blunt_matrix.txt
# 10, distance (C) x overhang (R) matrices of 5'overhang (one for count, one for fraction): *_merged_100bp+_cut_pattern_result_5overhang_matrix.txt
# 11, distance (C) x overhang (R) matrices of 3'overhang (one for count, one for fraction): *_merged_100bp+_cut_pattern_result_3overhang_matrix.txt
# 12, top strand distrance (C) X bottom strand distance (R) matrics of all cutting patterns combined: 25 columns X 36 rows (-5 to 30): *_merged_100bp+_cut_pattern_result_combined_matrix.txt
#
# Usage (on a computing cluster):
#	qsub circle_PAM_data_analysis_pipeline.sh R1_fastq_file R2_fastq_file output_prefix bowtie2_index reference_name reference_sequence
#
# To run this pipeline, you will need to supply the following 6 arguments (separate the arguments by space):
# 1, R1_fastq_file: fastq of R1. File name should end with "gz", "fastq" or "fq"
# 2, R2_fastq_file: fastq of R2.  File name should end with "gz", "fastq" or "fq"  
# 3 output_prefix: prefix of all the output file names
# 4, bowtie2_index: The path plus basename of the index for the reference genome. The basename is the name of any of the index files up to but not including the final .1.bt2 / .rev.1.bt2 / etc.   
# 5, reference_name: name of the reference sequence (it should be the same as the header (the content after ">") of the reference fasta file) 
# 6, reference_sequence: sequence of the reference sequence
# For example:
# qsub circle_PAM_data_analysis_pipeline.sh 44.1.fastq 44.2.fastq cas9_XL ~/work/projects/CirclePAM/reference/XL_XLL XL TGAGTGGGAAAGGTAATCGAACTGTNNNNNNNNNNCACGTTCAGAATGGCTTGGACTCGACATGGACTCCAATTGTATGGGATGCTAACACCAATCAAATGACACAGACTCCAGTCATGACCCTCAAGAGGAGAAAGCAGTGAGTGGGAAAGGTAATCGAACTGT
#############################################################################

# Check command line
[ $# -lt 6 ] && { echo "Usage: $0 R1_fastq_file R2_fastq_file output_prefix bowtie2_index reference_name reference_sequence"; exit 1; }

# parse command line, read input file names and output directory
args=("$@")
fq1=${args[0]}  # fastq of R1. File name should end with "gz", "fastq" or "fq"
fq2=${args[1]}	# fastq of R2.  File name should end with "gz", "fastq" or "fq"  
base=${args[2]} # base name for all the output files
index=${args[3]} # bowtie2 index of the reference sequence
ref_name=${args[4]} # name of the reference sequence
ref_seq=${args[5]} # sequence of the reference sequence 
 
echo -e "Step 1, trim adapter sequences and merge R1 and R2 reads...\n"
# Check if input file exists and readable
if test -r "$fq1" -a -f "$fq1"; then echo "Reading from first input file $fq1"; else echo "Unable to read from input file $fq1! Exit"; exit 1; fi
if test -r "$fq2" -a -f "$fq2"; then echo "Reading from second input file $fq2"; else echo "Unable to read from input file $fq2! Exit"; exit 1; fi

prefix1=${fq1%.gz}
prefix1=${prefix1%.fastq}
prefix1=${prefix1%.fq}
f1=$(basename $prefix1)

prefix2=${fq2%.gz}
prefix2=${prefix2%.fastq}
prefix2=${prefix2%.fq}
f2=$(basename $prefix2)

# prepare output files
trim1=$f1"_trimmed.fq.gz"  # adapter trimmed R1 fastq file (output of step 1)
trim2=$f2"_trimmed.fq.gz"  # adapter trimmed R2 fastq file (output of step 1)
merge=$base"_merged.fq.gz"  # output file for merged reads (output of step 1)
merge_good=$base"_merged_100bp+.fq"  # output file for merged reads with desired read length (output of step 2)
report1=$base"_merged_length_trimming_report.txt"   # read length trimming result summary (output of step 2)
sam=$base"_merged_100bp+_bt2.sam" # bowtie2 output sam file (output of step 3)
out=$base"_merged_100bp+_cut_pattern_result"  # prefix of the output files from step 4 
report2=$out"_summary.txt"; # result summary of step 4

# Activate conda environment where all the required programs are installed. In this case, the environment is called "bioinfo", you need to change it to the environment under your account
# If you don't use Conda, make sure all the required tools are either added to your system path ($PATH) or specify the path to individual programs when calling them. 
# Then comment out the following Conda lines:
# 	source activate ***
# 	source/conda deactivate
source activate bioinfo

## SeqPrep merging arguments
# -s <perform merging and output the merged reads to this file> 
# -E <write pretty alignments to this file for visual Examination>
# -o <minimum overall base pair overlap to merge two reads; default = 15>
# -m <maximum fraction of good quality mismatching bases to overlap reads; default = 0.020000>
# -n <minimum fraction of matching bases to overlap reads; default = 0.900000>
## Arguments for Adapter/Primer Trimming (Optional):
# -A: <forward read primer/adapter sequence to trim as it would appear at the end of a read (recommend about 20bp of this); Default is AGATCGGAAGAGCGGTTCAG (outdated illumina adapter sequence)
# -B: <reverse read primer/adapter sequence to trim as it would appear at the end of a read (recommend about 20bp of this); (should validate by grepping a file); default (genomic non-multiplexed adapter2) = AGATCGGAAGAGCGTCGTGT>
# -O <minimum overall base pair overlap with adapter sequence to trim; default = 10> 

overlap=50 # minimum overall base pair overlap to merge two reads. If the template length is between 100~150bp and the sequencing length is 150bp,  use 50 for minium overlap 
SeqPrep -f $fq1 -r $fq2 -1 $trim1 -2 $trim2 -A AGATCGGAAGAGCACACG -B AGATCGGAAGAGCGTCGT -O 5 -o $overlap -s $merge

echo -e "\nStep 2: Remove short merged reads <100bp..."
minlen=100; # minimum length of merged reads to keep. 
if [ ${merge: -3} == ".gz" ]; then
        gunzip -c $merge | awk -v len=$minlen -v file=$report1 'BEGIN {OFS = "\n"} {header = $0 ; getline seq ; getline qheader ; getline qseq ; total++; if (length(seq) >= len) {print header, seq, qheader, qseq; N++}};END{print "processed "total" records; "N" past the length filter of "len > file}' > $merge_good
else
        awk -v len=$minlen -v file=$report1 'BEGIN {OFS = "\n"} {header = $0 ; getline seq ; getline qheader ; getline qseq ; total++; if (length(seq) >= len) {print header, seq, qheader, qseq; N++}};END{print "processed "total" records; "N" past the length filter of "len > file}' $merge > $merge_good
fi
echo -e "check $report1 for read length trimming result.\n";

echo -e "Step 3: Map good-length merged reads to the reference sequence $index...\n"
# Suppress SAM records for reads that failed to align
bowtie2 --no-unal -x $index -U $merge_good -S $sam

echo -e "\nStep 4: Characterize cleavage pattern from the alignments..."
# You need to provide the path to the perl script (characterize_cut_pattern_from_circle_PAM.pl) in your system. Replace "/home/sunz/work/perl/" with your local path to the script
perl /mnt/home/sunz/perl/characterize_cut_pattern_from_circle_PAM.pl -i $sam -o $out -r $ref_name -s $ref_seq -f sam
echo -e "Output cut pattern and PAM sequences of individual reads to $out.txt.\nSee $report2 for a summary of the cleavage pattern results\n";

conda deactivate



