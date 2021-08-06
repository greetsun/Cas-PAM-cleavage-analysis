#!/bin/bash
#
#$ -cwd # running the job at the current working directory
#$ -j y # Stdout and stderr from this script will be merged (-j y) and written to job.sh.o<job_id>.
#$ -S /bin/bash
#$ -pe smp 1
#          ^ this is the number of CPU cores your program will use,
# or the amount of memory it needs / 4GB, whichever is larger, but
# not more than 16. E.g. 4 is for up to 4 CPUs and 16GB RAM.
# These limits aren't enforced, so it's up to you to avoid oversubscription.
#
# To get an e-mail when the job is done:
#$ -m e
#$ -M your@email.address
#


##################################
# Author: Zhiyi Sun (sunz@neb.com)
# Description: 
#	This pipeline script is for analyzing the Illumina sequencing data of the circle PAM control assay using BstXI (developed by Ryan Fuchs).
#	It takes pair-end fastq files as input, after merging the pair-end reads into full-length fragments and mapping them to the reference sequence (need to be prepared in advance and indexed using the bowtie2-build utility) 
#	it extracts the randomized region from the reads and calculates the background base frequencies, which will be used to normalize PAM sequence enrichment scores.
# In details, this scripts performs the following tasks: 			 
# 1, trim adapters and merge forward and reverse read* (seqprep)
# 	*NOTE: Before merging the reads, run fastqc to check data quality and read length (seqprep works best if both R1 and R2 reads contains the adapter sequence (partial is fine) at the 3' end and there is a good overlap between R1 and R2).
# 2, remove short reads (<100bp)
# 3, map merged reads to the reference (bowtie2)
# 4, Extract random sequence (10nt) from the good alignments*
#	 *criteria of a good alignment: primary alignment; MAPQ>=20; the entire read map to the reference (from start to end)
#	 This step will call another perl script: [extract_seq_from_sam.pl]
# 5, calculate the base frequency of the extracted random sequences
#
# Usage (on a computing cluster):
#	qsub circle_PAM_contral_BstXI_data_analysis_pipeline.sh R1_fastq_file R2_fastq_file output_prefix bowtie2_index reference_name reference_sequence_length start_position_random_sequence random_sequence_size
# For example:
# qsub circle_PAM_contral_BstXI_data_analysis_pipeline.sh 46.1.fastq 46.2.fastq BstXI_XL_control ~/work/projects/CirclePAM/reference/XL_XLL_BstXI XL_BstXI 136 90 10 
#
# To run this pipeline, you will need to supply the following 6 arguments (separate the arguments by space):
# 1, R1_fastq_file: fastq of R1. File name should end with "gz", "fastq" or "fq"
# 2, R2_fastq_file: fastq of R2.  File name should end with "gz", "fastq" or "fq"  
# 3, output_prefix: prefix of all the output file names
# 4, bowtie2_index: The path plus basename of the index for the reference genome. The basename is the name of any of the index files up to but not including the final .1.bt2 / .rev.1.bt2 / etc.   
# 5, reference_name: name of the reference sequence (it should be the same as the header (the content after ">") of the reference fasta file) 
# 6, reference_sequence_length: length of the reference sequence
# 7, start_position_random_sequence: 1-based left most position of the random sequence in the reference to be extracted
# 8, random_sequence_size: (optional) size of the random sequence to be extracted. Default=10 bases
#############################################################################


# Check command line
[ $# -lt 7 ] && { echo "Usage: $0 R1_fastq_file R2_fastq_file output_prefix bowtie2_index reference_name reference_length start_position_random_sequence random_sequence_size(optional)"; exit 1; }

# parse command line, read input file names and output directory
args=("$@")
fq1=${args[0]}  # fastq of R1. File name should end with "gz", "fastq" or "fq"
fq2=${args[1]}	# fastq of R2.  File name should end with "gz", "fastq" or "fq"  
base=${args[2]} # base name for all the output files
index=${args[3]} # bowtie2 index of the reference sequence
ref_name=${args[4]} # name of the reference sequence
ref_len=${args[5]} # length of the reference sequence
pos=${args[6]} # 1-based left most position of the random sequence in the reference to be extracted
len=${args[7]} # optional:length of the random sequence to be extracted. Default=10

if [ -z "$len" ]; then
        len=10;
fi
 
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
seqout=$base"_merged_100bp+_random_sequence.txt"  # output file of the extracted random sequences from the BstXI control libraries  (output of step 4)
logodata=$base"_merged_100bp+_random_sequence_logodata.txt" # logo data file of the weblogo program, it contains the base frequencies of each position of the input sequences (output of step 5)
basefreq=$base"_merged_100bp+_random_sequence_base_freq.txt" # output file of the nucleotide composition of the random sequences. (output of step 5)

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

echo -e "\nStep 4: Extract the random sequences from the good alignments..."
# You need to provide the path to the perl script (exract_seq_from_sam.pl) in your system. Replace "/home/sunz/work/perl/" with your local path to the script
perl /home/sunz/work/perl/extract_seq_from_sam.pl -i $sam -o $seqout -r $ref_name -l $ref_len -p $pos -n $len

echo -e "\nStep 5: Calculate nucleotide composition of the random sequences..."
weblogo -o $logodata -F logodata < $seqout

awk 'BEGIN{print "A\tC\tG\tT"};{if (NR>=9 && NR<=18) {nA=nA+$2; nC=nC+$3; nG=nG+$4; nT=nT+$5}}; END{total=nA+nC+nG+nT; printf("%.2f\t%.2f\t%.2f\t%.2f",nA/total*100, nC/total*100, nG/total*100, nT/total*100)}' $logodata > $basefreq
echo -e "Output base composition of the random sequence to $basefreq. DONE!"

# Deactivate the current Conda environment
conda deactivate

