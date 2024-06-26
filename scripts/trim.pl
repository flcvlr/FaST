use warnings;
use strict;

#### reading sequence data as : barcode, readname, read2, quality score


my @d;
while (<STDIN>) {
	chomp;
	@d=split("\t",$_);
	$d[4]=substr($d[2],0,9);
	if ($d[4] =~ /N/) {next;}
	$d[2] =~ s/A{5,}[CTGN]*A{5,}$//;
	if (length($d[2]) < 30) {next;}
	$d[2]=substr($d[2],9);
	$d[3]=substr($d[3],9,length($d[2]));
	### output: barcode, readname, umi, read2, quality score
	print "$d[0]\t$d[1]\t$d[4]\t$d[2]\t$d[3]\n";
	}


		
