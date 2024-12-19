#!/bin/perl

use warnings;
use strict;
use File::Basename;

my $L;
my $tile =basename($ARGV[0],".txt");
my $bc;
my $x;
my $y;
my %tile_map;
my %dup_bc;
my $counter=0;

open(my $in, "<", "$ARGV[1]/$tile.txt");

while (my $linein =<$in>) {
	chomp $linein;
	($x,$y,$bc)=(split("\t",$linein));
	$bc = reverse($bc); 
	$x = $x % 8000;
	$y = $y % 8000;
	$tile_map{$bc} = "$x\t$y\n";
}	


close $in;
open(my $out_idx, ">", "$ARGV[1]/$tile.idx");
open(my $out_gz, "|-", "pigz > $ARGV[1]/$tile.txt.gz");
while ((my $b,my $c)= each(%tile_map)){
	if ($counter < 650000) {
		print $out_idx "$b\n";
		}
	$counter ++;
	print $out_gz "$b\t$c";
}
my $ambigous =0;# keys(%dup_bc);

my $removed = 0;
#foreach my $k (keys(%dup_bc)) {
#	$removed += $dup_bc{$k};
#	}
print STDERR "removed $ambigous barcode sequences (corresponding to $removed pucks) from $tile\n";
close $out_gz;
close $out_idx;
