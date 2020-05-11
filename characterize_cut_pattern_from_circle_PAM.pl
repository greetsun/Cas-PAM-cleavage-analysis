#!/usr/bin/perl
#################################
## Name: characterize_cut_pattern_from_PAM_specific.pl
## AUTHOR: Zhiyi Sun
## Date: 11/8/2018
## Description:
## This script is part of the [circle_PAM_data_analysis_pipeline] which is designed for the circle PAM NGS assay (developed by Ryan F.) 
# to characterize PAM sequence and cleavage patterns of Class II CRISPR-Cas nuclease. 
# It reads a sam* file, determines the cleavage pattern and extract PAM sequences from good alignments; and finally reports the results to multiple output files.
#        Criteria for determine cutting pattern
# 		1, Read must map to the right reference (specified by user: $opt_r)
# 		2, primary alignment only (FLAG !=256; not $flag & 256)  
# 		3, MAPQ > 20
# 		4, alignment left most position <=25 (upstream of the PAM region)
# 		5, if the sequence before PAM region match the reference and
#       --if read length = reference original lenghth (XL=140; XXL=120):  blunt end
#       --if read length < reference original length: 3'overhang, overhang length= reference length - read length
#       --if read length > reference original length and the extra bases at the 3' end wrt the reference is same as the beginning of the read: 5' overhang, overhang length= read length - reference length
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
# *Sam file content
# 1, read name
# 2, Flag (0=map to the top strand; 16= map to the reverse strand; 256= not primary alignment): all the reported alignments have either 0 or 16 flag value
# 3, Reference sequence name (XL or XXL in this case)
# 4, POS: 1-based leftmost mapping position
# 5, MAPQ
# 6, CIGAR (the length of the sequence must equal the sum of lengths of M/I/S/=/X operations in CIGAR)
# 7, name of the mate (NA)
# 8, position of the mate (NA)
# 9, template length (* for single reads)
# 10, sequence of the read: All mapped reads are represented on the forward genomic strand. The bases are reverse complemented from the unmapped read sequence and the quality scores and cigar strings are recorded consistently with the bases. 
#
#
## Usage: perl characterize_cut_pattern_from_circle_PAM.pl -i input -o outputfile prefix; it will create 6 output files: the first is the cutting pattern results and PAM sequences for each read (prefix.txt); and the 2nd is a result summary (prefix_summary.txt); the 3rd to 6th files are cutting pattern matrices for blunt, 5overhang, 3overhang and combined. -r reference sequence name -s reference sequence -f file format: sam/bam, default :sam file
#
## Output matrix files:  
# For blunt end and 5'overhang, 3'overhang: Print top-strand cut distance to PAM (columns) X overhang (rows) matrices
# For combined matrix: print top-strand cut distance to PAM (columns) X bottom strand cut distance to PAM (rows). Note: for Cas9; PAM is downstream of the cut site, and D1 is the distance on the PAM strand, and D2 is the distance on the opposite strand; but for Cas12a, PAM is upstream of the cut site, therefore D2 is the distance to PAM on the PAM strand, and D1 is the distance on the bottom strand
# Each matrix file, there are 2 matrices: one is for the count; the 2nd one is for the fraction.
#####################################

use strict;
use Getopt::Std;
# parse command line
our ($opt_i, $opt_o, $opt_r, $opt_f, $opt_s);

getopt('iorsf');
unless($opt_i and $opt_o and $opt_r and $opt_s){
    print "Usage:\n";
    print "$0 -i input sam/bam file -o outputfile prefix; it will create 6 output files: the first is the cutting pattern results and PAM sequences for each read (prefix.txt); and the 2nd is a result summary (prefix_summary.txt); the 3rd to 6th files are cutting pattern matrices for blunt, 5overhang, 3overhang and combined. -r reference sequence name -s reference sequence -f file format: sam/bam, default :sam file\n\n";
    exit(1);
}

my $refseq=uc($opt_s);
my $reflen=length($refseq)-25; # Assume the reference is built by adding the beginning 25bp to the end of the sequence 
my $beforePAM=substr($refseq, 0, 25); # the sequence before PAM region

my $sam=$opt_i;
### if input is bam file, convert it to sam 
if ($opt_f=~/bam/i) {
	$sam =~s/$/\.sam/;
  	system("samtools view $opt_i > $sam");
}

open (IN, $sam)||die "Unable to read from temporary sam file: $!\n";
my $out1=$opt_o.".txt";
my $out2=$opt_o."_summary.txt";
my $m0=$opt_o."_blunt_matrix.txt";
my $m1=$opt_o."_5overhang_matrix.txt";
my $m2=$opt_o."_3overhang_matrix.txt";
my $m3=$opt_o."_combined_matrix.txt";

open (OUT,">$out1")||die "cant create file $out1\n";
open (SUMMARY,">$out2")||die "cant create file $out2\n";
print OUT "RID\tReference\tstrand\ttype\tdistance\toverhang\tPAM\n";
print SUMMARY "Reference sequence name: $opt_r\nReference sequence: $opt_s\nReference original length: $reflen\n";

open (M0,">$m0")||die "cant create file $m0\n";
open (M1,">$m1")||die "cant create file $m1\n";
open (M2,">$m2")||die "cant create file $m2\n";
open (M3,">$m3")||die "cant create file $m3\n";

my $total=0; # total sam alignment lines
my $good=0; # primary alignments whose MAPQ score is at least 20 and map to the correct reference sequence
my $upPAM=0; # good alignments whose start sites is upstream of the PAM region
my $rightseq=0; # among the $upPAM alignments, those upstream PAM sequences match the reference sequence
my $top=0; # of reads that map to the original top strand
my $bottom=0; # of reads that map to the original bottom strand
my $blunt=0; 
my $overhang5=0;
my $overhang5_putative=0; # putative 5' overhang counts (read length > reference length, but the last N bases do not match the first N bases)
my $overhang3=0;
my $overhang3_putative=0; # putative 3' overhang counts (read length < reference length, but the sequences do not match
my %blunt_hash=();  # key=cut distance, value=count
my %overhang5_hash=();
my %overhang3_hash=();		
my %combined_hash=();	# combine blunt and overhang patterns, characterize each cutting pattern by the distance to PAM on the top strand (D1) of the designed oligo; and the distance to PAM on the bottom strand (D2); if blunt: D1=D2; if 5'overhang; D1>D2; if 3'overhang, D1<D2. NOTE: for cas9; PAM is downstream of the cut site, and D1 is the distance on the PAM strand, and D2 is the distance on the opposite strand; but for Cas12a, PAM is upstream of the cut site, therefore D2 is the distance to PAM on the PAM strand, and D1 is the distance on the bottom strand
my $valid=0; # sum of $blunt, $overhang5 and $overhang3

while (<IN>){
	if (/^\@\w{2}\t/){
		#warn "skipping SAM header line:\t$_";
	  	next;
	}
	
	$total++;
	my ($rid, $flag,$rname,$pos, $mapq, $cigar, $seq) = (split (/\t/))[0, 1,2,3,4,5, 9];
	# Note: $pos is 1-based left-most mapping position
	$seq=uc($seq); 
	my $strand="";
	# check which strand the read mapped to
	
	if ($flag & 256 or $mapq <20 or $rname ne $opt_r) { next; }# skip non-primary alignment and low-quality alignments
	
	$good++;
		
	if ($flag & 16) { #'negative strand';
		$bottom++;	
		$strand="-";
	} else { # 'positive strand';
		$top++;
		$strand="+";
	}	
	
	my $read_length=length($seq);
	my $pam="";
	my $cutdistance=0; # distance of the cutting site from the 1st base of the PAM region
	my $overhang_length=0;

	if ($pos<=25 ) { # cutting site is upstream of the PAM region
		$upPAM++;
		# Check if the sequence match the reference sequence
		$cutdistance=25-$pos+1;
#		my $readbeforePAM= substr($seq, 0, $cutdistance);
#		my $refbeforePAM=substr($beforePAM, -$cutdistance);
#		print "This read is $rid\t$flag\t$pos\t$seq\n";
#		print "read length: $read_length\n";
#		print "read before pam: $readbeforePAM\n";
#		print "reference before pam: $refbeforePAM\n";
		if (substr($seq, 0, $cutdistance) eq substr($beforePAM, -$cutdistance)) {
			$rightseq++;
			if ($read_length eq $reflen) {
#				print "end type: blunt\n";
				# Case 1: blund end digestion
				$blunt++;
				$pam=substr($seq, $cutdistance, 10);
				if (defined $blunt_hash{$cutdistance} ) {
					$blunt_hash{$cutdistance}++;
				} else {
					$blunt_hash{$cutdistance}=1;
				}
				
				if (defined $combined_hash{$cutdistance}{$cutdistance}) {
					$combined_hash{$cutdistance}{$cutdistance}++;
				} else {
					$combined_hash{$cutdistance}{$cutdistance}=1;
				}
				print OUT "$rid\t$opt_r\t$strand\tblunt\t$cutdistance\t0\t$pam\n";
				
			} elsif ($read_length > $reflen) {
				# case 2: 5'overhang 
				$overhang_length=$read_length-$reflen;
#				my $overhang_seq5=substr($seq, 0, $overhang_length);
#				my $overhang_seq3=substr($seq, -$overhang_length);
#				print "end type: 5 overhang; overhang length=$overhang_length\n";
#				print "overhang sequence at 5 end: $overhang_seq5\n";
#				print "overhang sequence at 3 end: $overhang_seq3\n";
				if (substr($seq, 0, $overhang_length) eq substr($seq, -$overhang_length)) {
					$overhang5++;
					$pam=substr($seq, $cutdistance, 10);
					if (defined $overhang5_hash{$cutdistance}{$overhang_length} ) {
						$overhang5_hash{$cutdistance}{$overhang_length}++;
					} else {
						$overhang5_hash{$cutdistance}{$overhang_length}=1;
					}
					
					my $d2=$cutdistance-$overhang_length; # distance to PAM on the bottom strand
					if (defined $combined_hash{$cutdistance}{$d2} ) {
						$combined_hash{$cutdistance}{$d2}++;
					} else {
						$combined_hash{$cutdistance}{$d2}=1;
					}

					print OUT "$rid\t$opt_r\t$strand\t5overhang\t$cutdistance\t$overhang_length\t$pam\n";

				} else {
					$overhang5_putative++;
				}
			} else {
				# Case 3: 3' overhang
				$overhang_length=$reflen-$read_length;
				if (substr($refseq, 0, $pos-1-$overhang_length) eq substr($seq, -($pos-1-$overhang_length))) {
					$overhang3++; 				
				 	$pam=substr($seq, $cutdistance, 10);
				 	if (defined $overhang3_hash{$cutdistance}{$overhang_length} ) {
						$overhang3_hash{$cutdistance}{$overhang_length}++;
					} else {
						$overhang3_hash{$cutdistance}{$overhang_length}=1;
					}
					
					my $d2=$cutdistance+$overhang_length; # distance to PAM on the bottom strand
				 	if (defined $combined_hash{$cutdistance}{$d2} ) {
						$combined_hash{$cutdistance}{$d2}++;
					} else {
						$combined_hash{$cutdistance}{$d2}=1;
					}

					print OUT "$rid\t$opt_r\t$strand\t3overhang\t$cutdistance\t$overhang_length\t$pam\n";

				} else {
					$overhang3_putative++;
				}
			} # end of if ($read_length eq $reflen) 
		} # end of if (uc(substr($seq, 0, $cutdistance)) eq substr($ref, -$cutdistance)	 	
	
	} # end of if ($pos<=25 )
} # end of while (<IN>)


$valid=$blunt+$overhang5+$overhang3;

print SUMMARY "Processed $total alignments in the input file $opt_i; among them,\n";
print SUMMARY "Primary alignments that map to the correct reference and have a minimum MAPQ of 20: $good (top strand:$top; bottom strand: $bottom); among them,\n";
print SUMMARY "Good alignments whose start locations are upstream of the PAM region: $upPAM; among them,\n";
print SUMMARY "Good alignments whose upstream sequences match the reference sequences: $rightseq; among them,\n";
print SUMMARY "blunt cut: $blunt; among them,\n";
 
foreach my $key (sort { $a <=> $b } keys(%blunt_hash) )
{
    print SUMMARY "\tdistance to PAM: $key ($blunt_hash{$key})\n"
}

print SUMMARY "5' overhang cut:\n\tEnd sequence do NOT match: $overhang5_putative\n\tEnd sequence match: $overhang5; among them,\n";
print SUMMARY "\t\tDistance\tOverhang\tcount\n";
foreach (sort { $a <=> $b } keys(%overhang5_hash) )
{
    foreach my $thisoverhang (sort {$a<=>$b} keys %{$overhang5_hash{$_}})
    {
     	print SUMMARY "\t\t$_\t$thisoverhang\t$overhang5_hash{$_}{$thisoverhang}\n"
	}
}

print SUMMARY "3' overhang cut:\n\tEnd sequence do NOT match: $overhang3_putative\n\tEnd sequence match: $overhang3; among them,\n";
print SUMMARY "\t\tDistance\tOverhang\tcount\n";
foreach (sort { $a <=> $b } keys(%overhang3_hash) )
{
    foreach my $thisoverhang (sort {$a<=>$b} keys %{$overhang3_hash{$_}})
    {
     	print SUMMARY "\t\t$_\t$thisoverhang\t$overhang3_hash{$_}{$thisoverhang}\n"
	}
}

## Print distance (columns) X overhang (rows) matrices for blunt (overhang=0); 5'overhang (M1), 3'overhang (M2), and combined (M3) 
## For each category, print 2 matrices: one for count; one for fraction
#----------------
# Count matrix
#----------------
print M0 "distance";
print M1 "R:overhang/C:distance";
print M2 "R:overhang/C:distance";
for (my $i=1; $i<=25; $i++) {
	print M0 "\t$i";
	print M1 "\t$i";
	print M2 "\t$i";
	
}
print M0 "\ncount\t";
print M1 "\n";
print M2 "\n";

for (my $r=1; $r<=25; $r++) {
	print M1 "$r";
	print M2 "$r";
	
	if ($r == 25) {
		if (defined $blunt_hash{$r}) {
			print M0 "$blunt_hash{$r}\n\n";
		} else {
			print M0 "0\n\n";
		}
	} else {
		if (defined $blunt_hash{$r}) {
			print M0 "$blunt_hash{$r}\t";
		} else {
			print M0 "0\t";
		}
	}	
	
	for (my $c=1; $c<=25; $c++) {
		if (defined ($overhang5_hash{$c}{$r})) {
			print M1 "\t$overhang5_hash{$c}{$r}";
		} else {
			print M1 "\t0";
		}
		
		if (defined ($overhang3_hash{$c}{$r})) {
			print M2 "\t$overhang3_hash{$c}{$r}";
		} else {
			print M2 "\t0";
		}
				
		if ($c==25) {
			print M1 "\n";
			print M2 "\n";
		}
	}
}

# For combined matrix: column=distance on the top strand (1-25);row=distance on the bottom strand (-5 to 30);
print M3 "R:D(bottom)/C:D(top)";
# Print column label
for (my $i=1; $i<=25; $i++) {
	print M3 "\t$i";
}
print M3 "\n";

for (my $r=-5; $r<=30; $r++) {
	print M3 "$r";
	
	for (my $c=1; $c<=25; $c++) {
		if (defined ($combined_hash{$c}{$r})) {
			print M3 "\t$combined_hash{$c}{$r}";
		} else {
			print M3 "\t0";
		}
		
		if ($c==25) {
			print M3 "\n";
		}
	}
}


#----------------
# Fraction matrix
#----------------
print M0 "distance";
print M1 "\nR:overhang/C:distance";
print M2 "\nR:overhang/C:distance";
for (my $i=1; $i<=25; $i++) {
	print M0 "\t$i";
	print M1 "\t$i";
	print M2 "\t$i";
	
}
print M0 "\nFraction\t";
print M1 "\n";
print M2 "\n";

for (my $r=1; $r<=25; $r++) {
	print M1 "$r";
	print M2 "$r";
	
	if ($r == 25) {
		if (defined $blunt_hash{$r}) {
			my $f=sprintf("%.4f",$blunt_hash{$r}/$blunt);
			print M0 "$f\n";
		} else {
			print M0 "0\n";
		}
	} else {
		if (defined $blunt_hash{$r}) {
			my $f=sprintf("%.4f",$blunt_hash{$r}/$blunt);
			print M0 "$f\t";
		} else {
			print M0 "0\t";
		}
	}	
	
	for (my $c=1; $c<=25; $c++) {
		if (defined ($overhang5_hash{$c}{$r})) {
			my $f=sprintf("%.4f",$overhang5_hash{$c}{$r}/$overhang5);
			print M1 "\t$f";
		} else {
			print M1 "\t0";
		}
		
		if (defined ($overhang3_hash{$c}{$r})) {
			my $f=sprintf("%.4f",$overhang3_hash{$c}{$r}/$overhang3);
			print M2 "\t$f";
		} else {
			print M2 "\t0";
		}
				
		if ($c==25) {
			print M1 "\n";
			print M2 "\n";
		}
	}
}

# For combined matrix: column=distance on the top strand (1-25);row=distance on the bottom strand (-5 to 30);
print M3 "\nR:D(bottom)/C:D(top)";
for (my $i=1; $i<=25; $i++) {
	print M3 "\t$i";
	
}
print M3 "\n";

for (my $r=-5; $r<=30; $r++) {
	print M3 "$r";
	
	for (my $c=1; $c<=25; $c++) {
		if (defined ($combined_hash{$c}{$r})) {
			my $f=sprintf("%.4f",$combined_hash{$c}{$r}/$valid);
			print M3 "\t$f";
		} else {
			print M3 "\t0";
		}
		
		if ($c==25) {
			print M3 "\n";
		}
	}
}


close IN;
close OUT;
close SUMMARY;
close M0;
close M1;
close M2;
close M3;
exit;
