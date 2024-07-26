use warnings;
use strict;

### ARGV[0] = fasta

my $chr ="";
my $len = 0;


open(my $genome, "<","$ARGV[0]");

while (<$genome>) {
	chomp;
	if ($_ =~ /^>/) {
		if ($chr ne "") {
			print "$chr\t$len\n";
		} 
		$_ = substr($_,1); 
		$chr = (split(" ",$_))[0];
		$len=0;
	} 
	else {
		$len += length($_);
	}
}
if ($chr ne "") {
	print "$chr\t$len\n";
} 

exit;
