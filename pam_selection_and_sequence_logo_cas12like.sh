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
#$ -M your@email.com
#

# Check command line
#[ $# -eq 0 ] && { echo "Usage: $0 inputfile"; exit 1; }
# or [ $# -ne N ] && { echo "Usage: $0 argument"; exit 1; }

###########################
# Name: pam_selection_and_sequence_logo_cas12like.sh 
# Author and contact: Zhiyi Sun (sunz@neb.com)
#
# NOTE: This script is modified from pam_selection_and_sequence_logo.sh; It is for Cas12 like proteins that cut downstream of the PAM sequences. Therefore, the PAM sequence is reverse complement of the reference sequence. 
#	(1) in step 1, output the reverse complemented sequence	
#	(2) use the reverse-complemented background base frequencies for correction in Step 2 and Step 3. 
#	(3) in step 2, modify the positions of base-frequency matrix from (1..10) to (-10..-1)
#
#----------------
# Description: 
#----------------
#	This program takes the output file (*cut_pattern_result_summary.txt) of the circle_PAM_data_analysis_pipeline as input file, 
#	selects the PAM sequences for a specific cut pattern (specified by user), then 
#	calculates the per-position base frequency as well as the Position Weight Matrix (PWM); and finally 
#	generates a regular sequence logo and a PWM sequence logo.  
#	The PWM is created by comparing the base frequencies of the PAM sequence to the background base frequencies.
#	You can provide the background base frequencies for A,C,G,T (as percentage) in the command line;
#	or the program will assume equal base composition as A=C=G=T=25(%)
#
# An example of the content of the input file (*cut_pattern_result_summary.txt) 
# RID     Reference       strand  type    distance        overhang        PAM
# NS500355:NS500355:HVWMYAFXX:1:11101:10028:1465  XL      +       blunt   4       0       GGTTCCAGGG
# NS500355:NS500355:HVWMYAFXX:1:11101:10066:14846 XL      +       5overhang       5       1       GGTTAAAGGT
# NS500355:NS500355:HVWMYAFXX:1:11101:10124:13112 XL      -       blunt   3       0       TGGAAGGATT
#
# Selection of the PAM sequence are based on the combination of "type", "distance" and "overhang"
#
#--------------
# Output files
#--------------
# This program will produce the following 6 output files:
# 1, a sequence file for the selected PAM: *_pam_seq.txt
# 2, sequence logo of the selected PAM: *_pam_seq_logo.eps
# 3, data used for generating the sequence logo: *_pam_seq_logodata.txt
# 4, a table of per-position base frequency of the selected PAM sequences (this table is used for making the position weight matrix in R): *_pam_seq_base_freq.txt
# 5, a regular sequence logo of the selected PAM sequences: *_pam_seq_logo.eps
# 6, a Position Weight Matrix sequence logo of the selected PAM sequences: *_position_weight_matrix_logo.pdf

#--------------
# Usage
#--------------
# qsub pam_selection_and_sequence_logo.sh inputfile output_prefix end_type distance_to_PAM overhang_size background_A_perc(optional) background_C_perc(optional) background_G_perc(optional) background_T_perc(optional)
# Where:
# inputfile: input file name, which is the output file (*cut_pattern_result_summary.txt) of the circle_PAM_data_analysis_pipeline 
# output_prefix: prefix of all the output file names
# end_type: type of the cleavage ends, choose from "blunt", "5overhang" and "3overhang" (corresponds to the 4th column in the input file)
# distance_to_PAM: the distance from the cut site w.r.t the top strand of the reference to the 1st base of the randomized region (PAM) in the substrate (corresponds to the 5th column in the input file)
# overhang_size: overhang size (corresponds to the 6th column in the input file)
# background_A_perc(optional):optional: background A percentage. value ranges between 0 to 100*
# background_C_perc(optional):optional: background C percentage. value ranges between 0 to 100*
# background_G_perc(optional):optional: background G percentage. value ranges between 0 to 100*
# background_T_perc(optional):optional: background T percentage. value ranges between 0 to 100*
# * NOTE: percentages of the 4 background bases have to add up to 100
############################

# Check command line
[ $# -lt 5 ] && { echo "Usage: $0 inputfile output_prefix end_type distance_to_PAM overhang_size background_A_perc(optional) background_C_perc(optional) background_G_perc(optional) background_T_perc(optional)"; exit 1; }

# parse command line, read input file names and output directory
args=("$@")
input=${args[0]} # input file name, which is the output file (*cut_pattern_result_summary.txt) of the circle_PAM_data_analysis_pipeline
prefix=${args[1]} # prefix of the output file names
end=${args[2]} # choose from "blunt", "5overhang" and "3overhang" (corresponds to the 4th column in the input file)
distance=${args[3]} # the distance from the cut site w.r.t the top strand of the reference to the 1st base of the randomized region (PAM) in the substrate (corresponds to the 5th column in the input file)
overhang=${args[4]} # overhang size (corresponds to the 6th column in the input file)
pA=${args[5]}  # optional: background A percentage. value ranges between 0 to 100
pC=${args[6]}  # optional: background C percentage. value ranges between 0 to 100
pG=${args[7]}  # optional: background G percentage. value ranges between 0 to 100
pT=${args[8]}  # optional: background T percentage. value ranges between 0 to 100
#bgf=${args[5]} # this argument is optional. If missing, sequence logo will be generated using equiprobable background base composition. If using a custom base composition, follow this format: '"{'A':10, 'C':40, 'G':40, 'T':10}"'

# Check if input file exists and readable
if test -r "$input" -a -f "$input"; then echo "Reading from input file $input"; else echo "Unable to read from input file $input! Exit"; exit 1; fi

# Prepare output file names
seq=$prefix"_pam_seq.txt";
logo=$prefix"_pam_seq_logo.eps";
logodata=$prefix"_pam_seq_logodata.txt"
matrix=$prefix"_pam_seq_base_freq.txt"
pwmlogo=$matrix"_position_weight_matrix_logo.pdf" 

echo "Step 1: select the PAM sequences for the following cut pattern:"
echo -e "End type=$end\nDistance to PAM=$distance\nOverhang size=$overhang"
echo "NOTE: Output the reverse complemented sequence for Cas12-like proteins"

# NOTE: for Cas12 like proteins, the PAM sequence is the reverse complement of the reference sequence, therefore output the converted sequence

awk -v E=$end -v D=$distance -v O=$overhang '{if ($4==E && $5==D && $6==O) print toupper($NF)}' $input | rev | tr ACGT TGCA > $seq


# Activate conda environment where all the required programs are installed. In this case, the environment is called "bioinfo", you need to change it to the environment under your account
# If you don't use Conda, make sure all the required tools are either added to your system path ($PATH) or specify the path to individual programs when calling them. 
# Then comment out the following Conda lines:
# 	source activate ***
# 	source/conda deactivate

source activate bioinfo

echo -e "\nStep 2:  calculate base frequency at each position and generate sequence logo for the selected PAM sequence..."

if [ -z "$pA" ] 
then
	echo "Make sequence logo using equiprobable background base composition" 
	weblogo -o $logodata -F logodata < $seq
	weblogo < $seq > $logo
else
	sum=$(echo $pA + $pC + $pG + $pT | bc); 
	sumint=$( printf "%.0f" $sum )
	if [[ "$sumint" -eq 100 ]]
	then
		## NOTE: for Cas12-like proteins, the PAM sequence is reverse complement to the reference sequence, therefore, need to use the reverse complement base frequency for correction
		bgf="{'A':$pT, 'C':$pG, 'G':$pC, 'T':$pA}"
	#	bgf="{'A':$pA, 'C':$pC, 'G':$pG, 'T':$pT}"  
	#	echo $bgf

    	echo -e "Make sequence logo using custom background base composition as: $bgf"
		weblogo --composition "$bgf" --first-index -10 < $seq > $logo  # number the first position in the sequence as "-10"
		weblogo --composition "$bgf" -o $logodata -F logodata < $seq
	else 
		echo "Percentages of background A, C, G, T do not add up to 100. Please rerun the program with valid percentage values. Exit"; 
		exit 1;
	fi	
fi

# Convert the logodata to a base-frequency matrix. For Cas12 like proteins, modify the positions from (1..10 ) to (-10..-1)
awk 'BEGIN{print "Pos\tA\tC\tG\tT"};{if (NR>9) print last; last=$1-11"\t"$2"\t"$3"\t"$4"\t"$5}' $logodata > $matrix
echo -e "Output per-position base composition of the PAM sequence to $matrix"
echo "See $logo for the PAM sequence logo";

echo -e "\nStep 3: generate Position Weight Matrix sequence logo for the selected PAM sequence..."
if [ -z "$pA" ] 
then 
	pA=25; pC=25; pG=25; pT=25;
fi 

## NOTE: for Cas12-like proteins, the PAM sequence is reverse complement to the reference sequence, therefore, need to use the reverse complement base frequency for correction
R CMD BATCH --no-save --quiet --slave "--args $matrix $pT $pG $pC $pA" /home/sunz/work/projects/CirclePAM/script/generate_PWM_seqlogo.r
echo -e "Done! See $pwmlogo for the Position Weight Matrix sequence logo\n" 

conda deactivate