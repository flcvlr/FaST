#!/bin/perl

use warnings;
use strict;
$|=1;
my $tile;
my $linein2;
my $L;
my $x;
my $y;
my $bc;
my $counter;
open (my $data_in, "-|", "pigz -cd $ARGV[0]") ;
my $firstline=  <$data_in>;
($L,$tile)=(split(":",$firstline))[3,4];
my $file=$L."_".$tile;
my $current_tile=$tile;
open(my $out, ">", "$ARGV[1]/$file.txt");
chomp $firstline;
my $secondline=<$data_in>;
print $out "$firstline\t$secondline";
my $bla=<$data_in>;
$bla=<$data_in>;
while (my $linein = <$data_in>) {
	if ((split(":",$linein))[4] ne $current_tile) {
		close $out;
		print "$file\n";
		($tile)=(split(":",$linein))[4];
		$L=(split(":",$linein))[3];
		$file=$L."_".$tile;
		$current_tile=$tile;
		open ($out, ">", "$ARGV[1]/$file.txt");
		$counter ++;
	}
	chomp $linein;
	$linein2 = <$data_in>;
	print $out "$linein\t$linein2";
	$bla = <$data_in>;	
	$bla = <$data_in>;
}

close $out;
print "$current_tile\n";
print STDERR "processed $counter tiles\n";
exit;

