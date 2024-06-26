use warnings;
use strict;

my @d;
my $name;
my $start;
open (my $out, ">","$ARGV[0]/$ARGV[1].introns.bed");
while (my $linein = <STDIN>) {
	if ($linein !~/^chr/) {next;}	
	@d =split("\t",$linein);
	$linein =~ /gene_name "([^"]+)"/;
	$name=$1."#_INTRON";
	$start = $d[3] - 1;
	print $out "$d[0]\t$start\t$d[4]\t$name\t.\t$d[6]\n";
	
}
close $out;
