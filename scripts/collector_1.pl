use warnings;
use strict;
my $scale_factor;
my $tile_w;
if ($ARGV[2] eq "Illumina") {
	$scale_factor=15.588;
	$tile_w=3000;
}

if ($ARGV[2] eq "Stereo-seq") {
	$scale_factor=1;
	$tile_w=8000;
}
	

my @Aobc;
my $line;
my $frag;
my @num;
my $h;
my @map;
my %kmers;

$map[ord("A")]="00";
$map[ord("C")]="01";
$map[ord("G")]="10";
$map[ord("T")]="11";


my %pam=("00"=>"A","01"=>"C","10"=>"G","11"=>"T");

### get sequences and build DB
open(my $in, "<", "$ARGV[1]/$ARGV[0]_barcodes");
my $nts = 9;
my $bits = $nts*2;
while ($line=<$in>)  {
	if ($line eq "END OF RECORDS\n") {last;}
	$frag=substr($line,0,$nts,"");
	my $hashed = oct(join("","0b", @map[unpack("C$nts", $frag)]));
	if (!exists $Aobc[$hashed] || index($Aobc[$hashed] ,$line) == -1) {
		$Aobc[$hashed] .= $line;
	}
}

close $in;
### collect tile barcodes and pass back how many were in DB

open (my $rep_out, ">>", "$ARGV[1]/reports");
$rep_out->autoflush(1);

my $out_report;
my $found =0;

while (($line=<STDIN>) ne "END OF RECORDS\n" ) {
	if ($line eq "END OF FILE\n") {
		print $rep_out "$found\n"; 
		$found=0; 
		next;
	}
	$frag=substr($line,0,$nts,"");
	my $hashed = oct(join("","0b", @map[unpack("C$nts", $frag)]));	
	if (exists $Aobc[$hashed] && index($Aobc[$hashed] ,$line) != -1) {$found ++;}
}


my $retained=0;
my @bc_array;
my %amb;
my $ambiguous =0;
my $tile;
my $pos;
while (($line=<STDIN>) ne "END OF RECORDS\n" ) {
	my ($bc,$x,$y) = split("\t",$line);
	if (exists $amb{$bc}) {$ambiguous ++; next;}
	$frag=substr($bc,0,$nts,"");
	if ($frag eq "NEW_TILE_") {$tile = $x;next;}
	my $hashed = oct(join("","0b", @map[unpack("C$nts", $frag)]));
	if (!exists $Aobc[$hashed] || index($Aobc[$hashed] ,$bc) == -1) {next;}
	my $coord = sprintf("%08d",int(($x/$scale_factor)+.5) + (int(($y/$scale_factor)+.5))*$tile_w);
	if (!exists($bc_array[$hashed])) {$bc_array[$hashed] = "$bc\tTL:Z:$tile\tCO:Z:$coord\n"; $retained ++; next;}
	if ( ($pos = index($bc_array[$hashed] ,$bc)) == -1) {$bc_array[$hashed] .= "$bc\tTL:Z:$tile\tCO:Z:$coord\n"; $retained ++; next;}
	$amb{$frag.$bc}++;
	$ambiguous +=2;
	substr($bc_array[$hashed],$pos,55)="";
	$retained -= 1;
	}                              


print $rep_out "$retained\t$ambiguous\n";
close $rep_out;

my @dinucleotide= split("",$ARGV[0]);
open(my $barcodes, ">>", "$ARGV[1]/$dinucleotide[0]_barcodes");
$barcodes->autoflush(1);
foreach my $i (0..$#bc_array) {
	if (exists($bc_array[$i]) ) {
		foreach my $bc (split("\n",$bc_array[$i])) {
			my $reverted= join("",@pam{unpack("A2" x $nts,sprintf ("%0${bits}b", $i))});
			print $barcodes $dinucleotide[1].$reverted.$bc."\n";
			
		}
	}
}	

exit;






