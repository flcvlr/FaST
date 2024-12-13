use warnings;
use strict;

my @d;
my $name;
my $start;
open (my $gtf, "<", "$ARGV[2]");
open (my $out, ">","$ARGV[0]/$ARGV[1].exons.bed");
while (my $linein = <$gtf>) {
	if ($linein =~/^#/) {next;}
	@d =split("\t",$linein);
	if ($d[2]eq "exon") {
	$linein =~ /gene_name "([^"]+)"/;
	$name=$1;
	$start = $d[3] - 1;
	print $out "$d[0]\t$start\t$d[4]\t$name\t.\t$d[6]\n";
	}
}
close $gtf;
close $out;

