use warnings;
use strict;

#### this script reads input from STAR in sam format, collects BC statistics and splits the output (after conversion to bed format) based on tile to multiple parallel processes.

my %H_of_tile_fh;

while (<STDIN>) {
	$_ =~/RG:Z:(.+$)/;
	if (! exists ($H_of_tile_fh{$1})) {
		open (${H_of_tile_fh{$1}}, "|-", "perl $ARGV[1]/bam_to_bed_w_tags.pl $ARGV[0] ${1} $ARGV[1]" );
		}
#	$H_of_BC_counts{$1}++;
	print {$H_of_tile_fh{$1}} $_;
	}

foreach my $tile (keys %H_of_tile_fh) {
	print  {$H_of_tile_fh{$tile}} "END_OF_RECORDS\n";
	}
	
exit;	
	
	
