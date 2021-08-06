#!/usr/bin/perl
#################################
## Name: extract_seq_from_sam.pl
## Author and contact: Zhiyi Sun (sunz@neb.com)
## Description: this program will scan user input sam/bam file; use the good alignments to extract sequence from the read sequence
## Sam file content
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
## Criteria for a good alignment 
# 1, Read must map to the right reference (specified by user: $opt_r)
# 2, pirmary alignment only (FLAG !=256; not $flag & 256)  
# 3, MAPQ > 20
# 4, full alignment (leftmost mapping position=1, read length=reference length)

# Usage: perl extract_seq_from_sam.pl -i input -o outputfile -r reference sequence name -l reference sequence length -f file format: sam/bam, default :sam file
#
#####################################

use strict;
use Getopt::Std;
# parse command line
our ($opt_i, $opt_o, $opt_r, $opt_f, $opt_l, $opt_p, $opt_n);

getopt('iorlfpn');
unless($opt_i and $opt_o and $opt_r and $opt_l and $opt_p and $opt_n){
    print "Usage:\n";
    print "$0 -i input sam/bam file -o outputfile -r reference sequence name -l reference sequence length -p 1-based-left-most position of the sequence to be extracted -n length of the sequence to be extracted -f file format: sam/bam, default :sam file\n\n";
    exit(1);
}

my $sam=$opt_i;
### if input is bam file, convert it to sam 
if ($opt_f=~/bam/i) {
	$sam =~s/$/\.sam/;
  	system("samtools view $opt_i > $sam");
}

open (IN, $sam)||die "Unable to read from temporary sam file: $!\n";

open (OUT,">$opt_o")||die "cant create file $opt_o\n";

my $total=0; # total sam alignment lines
my $good=0; # primary alignments whose MAPQ score is at least 20 and map fully to the correct reference sequence

while (<IN>){
	if (/^\@\w{2}\t/){
		#warn "skipping SAM header line:\t$_";
	  	next;
	}
	
	$total++;
	my ($flag,$rname,$pos, $mapq, $read) = (split (/\t/))[1,2,3,4,9];
	# Note: $pos is 1-based left-most mapping position
	$read=uc($read);
	my $Rlength=length($read); 
	
	if ($flag & 256 or $mapq <20 or $rname ne $opt_r or $pos ne 1 or $Rlength ne $opt_l) { next; }# skip non-primary alignment and low-quality alignments
	
	$good++;
	
	my $seq=substr($read, $opt_p-1, $opt_n);
	print OUT "$seq\n";

} # end of while (<IN>)

print "Processed $total alignments in the input file $opt_i; found $good good alignments and extracted $opt_n bases from position $opt_p; output to $opt_o. \n\n"; 

close IN;
close OUT;
exit;