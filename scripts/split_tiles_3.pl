use warnings;
use strict;

#### this script reads input from STAR in sam format, collects BC statistics and splits the output (after conversion to bed format) based on tile to multiple parallel processes.

my @A_of_fh;
my @A_of_pids;
my %H_of_pids;
my %H_of_fh;
my $cores=20;
open(my $in, "-|", "samtools view -@ 3 $ARGV[0]/Aligned.bam");
my @tiles = `cut -f 5 $ARGV[0]/barcodes_per_tile.txt `;
foreach my $n (@tiles) {
	chomp $n;
	system ("mkfifo $ARGV[0]/$n");
	$H_of_pids{$n} = open($H_of_fh{$n}, "|-", "perl $ARGV[1]/DGE.pl $ARGV[0] $n $ARGV[2] $ARGV[3]");
}

foreach my $i (0..($cores-1)) {
	$A_of_pids[$i] = open($A_of_fh[$i], "|-", "perl $ARGV[1]/intersect.pl $ARGV[0] $ARGV[1] $ARGV[2] $ARGV[3] $i" );
}
my $tile;
while (<$in>) {
	$tile=(substr($_,-3,2) % $cores);
	print {$A_of_fh[$tile]} $_;
}

foreach my $i (0..$#A_of_fh) {
	if (exists $A_of_fh[$i]) {
		close $A_of_fh[$i];
		waitpid($A_of_pids[$i],0);
	}
}

foreach my $n (@tiles) {
	close $H_of_fh{$n};
	waitpid ($H_of_pids{$n},0);
}
my $t = (split(" ",localtime))[3];
print "$t\tFinished preparing DGE data\n";
close $in;	
exit;	
	
	
