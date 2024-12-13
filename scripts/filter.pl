use warnings;
use strict;

open(my $barcodes, "<", "$ARGV[1]/$ARGV[0]_barcodes");
my $info_len;
if ($ARGV[2] eq "Illumina") {$info_len=27;}
if ($ARGV[2] eq "Stereo-seq") {$info_len=25;}
my @map;
my %kmers;
$map[ord("A")]="00";
$map[ord("C")]="01";
$map[ord("G")]="10";
$map[ord("T")]="11";

my @bc_array;
my $nts=9;
my $len=24-$nts;

while (my $line = <$barcodes>) {
	my $frag=substr($line,0,$nts,"");
	my $hashed = oct(join("","0b", @map[unpack("C$nts", $frag)]));
	$bc_array[$hashed] .= $line;
}

open (my $rep_out, ">>", "$ARGV[1]/input.sam");
$rep_out->autoflush(1);
my $aligned =0;
my $out ="";
my $pos;
open(my $in, "<","$ARGV[1]/$ARGV[0]_reads");
while (my $line=<$in>) {
#	if ($line eq "END OF_RECORDS\n") {last;}
	my $bc=substr($line,-25,24);
	my $frag=substr($bc,0,$nts,"");
	my $hashed = oct(join("","0b", @map[unpack("C$nts", $frag)]));
	if (exists $bc_array[$hashed] && ( ($pos = index($bc_array[$hashed], $bc)) != -1) ) {
		substr ($line,-1,1)=substr($bc_array[$hashed],($pos+$len),$info_len);
		$out .= $line;
		if (length($out) > 3000) {print $rep_out $out; $out ="";}
		$aligned ++;
	}
}

print $rep_out $out;
close $barcodes;
close $rep_out;
