#!/bin/perl

use warnings;
use strict;


my $L;
my $tile =$ARGV[0];
my $bc;
my $x;
my $y;
my %tile_map;
my %dup_bc;
my $counter=0;

open(my $in, "<", "$ARGV[1]/$ARGV[0].txt");

while (my $linein =<$in>) {
	($x,$y)=(split(":",$linein))[5,6];
	$bc = reverse(substr((split("\t",$linein))[1],5,25)); 
	if ($bc !~ /N/) {
		$bc =~ tr/ACGT/TGCA/;
		if (exists $tile_map{$bc}) {
			$dup_bc{$bc}+=2;
			delete $tile_map{$bc};
		}
		else {
			if (!exists $dup_bc{$bc}) {
				$y = (split(" ",$y))[0];
				$tile_map{$bc} = "$x\t$y\n";
			}
			else {
				$dup_bc{$bc}++;
			}
		}
	}	
}

close $in;
unlink "$ARGV[1]/$ARGV[0].txt";
open(my $out_idx, ">", "$ARGV[1]/$tile.idx");
open(my $out_gz, "|-", "pigz > $ARGV[1]/$tile.txt.gz");
while ((my $b,my $c)= each(%tile_map)){
	if ($counter < 20000) {
		print $out_idx "$b\n";
		}
	$counter ++;
	print $out_gz "$b\t$c";
}
my $ambigous = keys(%dup_bc);

my $removed = 0;
foreach my $k (keys(%dup_bc)) {
	$removed += $dup_bc{$k};
	}
print STDERR "removed $ambigous barcode sequences (corresponding to $removed pucks) from $tile\n";
close $out_gz;
close $out_idx;
