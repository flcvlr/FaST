use warnings; 
use strict;
my $pad="";
my $x;
my $y;
my $in;
open (my $tiles, "-|", " cat $ARGV[0]/FQ*/barcodes_per_tile.txt");
open (my $out, "|-", " sort -u > $ARGV[0]/tiles_info");
while (my $linein = <$tiles>) {
	chomp $linein;
	$in=(split("\t",$linein))[4];
	my $swath= substr($in,3,1);
	my $tile = substr($in,4,2);
	my $side = substr($in,2,1);
	$tile =~ s/^0//;
	$x= int($swath*3748.77777);
	if (($swath % 2) == 1) {$y = int(4031.33333 * $tile); };
	if (($swath % 2) == 0) {
		if ($side == 1) {$y = int(4031.33333 * $tile)-689;}
		if ($side == 2) {$y = int(4031.33333 * $tile)+689;}
		}
	print $out "$in,$x,$y\n";

}
close $tiles;
close $out;
exit;


	
