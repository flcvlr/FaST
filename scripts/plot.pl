use warnings;
use strict;
use GD;

if (substr($ARGV[0],-1,1) eq "/") {chop $ARGV[0];}

my $x_min=40000000;
my $x_max=0;
my $y_min=40000000;
my $y_max=0;
my @intensities;
my @logs=`ls $ARGV[0]/logs/*.tile.log`;

foreach my $log (@logs) {
	chomp $log;
	open (my $log_file, "<", "$log"); 
	my $bla=<$log_file>;
	$bla=<$log_file>;
	my $coord_info=<$log_file>;
	chomp $coord_info;
	my @d=split("\t",$coord_info);
	if ($d[1] < $x_min) {$x_min=$d[1];}
	if ($d[2] > $x_max) {$x_max=$d[2];}
	if ($d[3] < $y_min) {$y_min=$d[3];}
	if ($d[4] > $y_max) {$y_max=$d[4];}
	close $log_file;
	}

$x_max-=$x_min;
$y_max-=$y_min;

open (my $scores, "<", "$ARGV[0]/dge/whole_dataset.txt"); 
my $header=<$scores>;
while (<$scores>) {
#	if (substr($_,0,1) eq "x") {next;}
	my @d=split("\t",$_);
#	my $x=$d[0]-$x_min;
#	my $y=$d[1]-$y_min;
	$intensities[(($d[1]-$y_min)*$x_max)+($d[0]-$x_min)]+=$d[3];
}
close $scores;
my $scale;
my @distribution;
my $existing;
my $observed=0;
foreach my $score (@intensities) {
	if ($score) {$distribution[$score]+=1;$existing +=1;}
}
foreach my $i (reverse (0..$#distribution)) {
	if ($distribution[$i]) {$observed+=$distribution[$i];}
	if (($observed/$existing) > 0.05) {
		$scale = 255/$i; 
		last;
	}
}
if ($scale < 1) {$scale=1;}
my $img = GD::Image->new($x_max,$y_max);
my $white = $img->colorAllocate(255,255,255);
my @colors;
foreach my $c (0..255) {
	$colors[$c] = $img->colorAllocate($c,$c,$c);
	}
foreach my $index (0..$#intensities) {
	if($intensities[$index]) {
		my $col = 255- ($intensities[$index]*$scale);
		if ($col < 0) {$col=0;}
		if ($col==255) {next;}		
		my $x = $index % $x_max;
		my $y = int($index/$x_max);
		$img->setPixel($x,$y,$colors[$col]);
	}
}
open (my $out, ">", "$ARGV[0]/$ARGV[0]_reads_density.png") or die "$!\n";
binmode $out;
print $out $img->png;
close $out;

