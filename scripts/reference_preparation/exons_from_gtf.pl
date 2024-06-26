use warnings;
use strict;

my @d;
my $name;
my $start;
while (my $linein = <STDIN>) {
	if ($linein !~/^chr/) {next;}	
	@d =split("\t",$linein);
	if ($d[2]eq "exon") {
	$linein =~ /gene_name "([^"]+)"/;
	$name=$1;
	$start = $d[3] - 1;
	print "$d[0]\t$start\t$d[4]\t$name\t.\t$d[6]\n";
	}
}
