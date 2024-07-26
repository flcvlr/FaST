use warnings;
use strict;

#### reading sequence data as : readname\tread2_seq\tread_1_seq\tread2_score

my $bc;
my $name;
my $seq;
my $umi;
my $score;
while ( my $line= <STDIN>) {
	chomp $line;
	($name,$seq,$bc,$score)=split("\t",$line);
	$name =(split(" ",$name))[0];
	substr($name,0,1)="";
	$bc =substr($bc,2,25);
	$umi=substr($seq,0,9);
	if ($umi =~ /N/) {next;}
	$seq=substr($seq,9);
	$score=substr($score,9,length($seq));
		
	### output: 
	
	print "$name\t4\t*\t0\t0\t*\t*\t0\t0\t$seq\t$score\tCB:Z:$bc\tMI:Z:$umi\n";
	}


		
