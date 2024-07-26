use warnings;
use strict;
use File::Basename;
$|=1;
### ARGV: 0=R1_reads; 1=R2_reads; 2=sample; 3=script_directory; 4=tiles_map_directory; 5= ratio

### requires samtools, zcat, sed, cut, grep 

my $counter=0;  ### num of R1 records processed
my $retained=0;  ### num of retained reads ( w/ non ambiguous barcodes in "abundant" tiles)
my $discarded=0; ### discarded barcodes
my @barcode_files = `ls $ARGV[4]/*.idx`; 
my $ambiguous=0; ### number of ambiguous barcodes found
my $unmapped;  ### number of unmapped barcodes
my $tile; 
my %H_of_A_of_amb_bc;
my %H_of_BC;  ### main hash (dictionary) containing barcode info 
my %H_of_BC_counts;  ### hash containing barcode counts 
my $barcode_linein;  ### input line during barcode preparsing
my $R1_linein;  ### input line during bam file preparation
my @data;  ### array of input data after splitting input line
my $bc;  
my $x;
my $y;
my $t;
my $umi;
my $counts;
my $found=0;  ### counter to use in intermediate step, to count how many barcodes were found for each tile
my $all;    ### number of barcodes read from tile info file
my $ratio;
my %H_of_tiles;  ### hash (dictionary) containing tile information
my $candidate;
my %H_of_corrected_bc;
my $HD_recovered_bc;
my $HD_recovered_entries;
my $matched;
#### read barcodes from file R1 and store in hash

$ARGV[0] =~ y/,/ /;

$t = localtime;
print STDERR ((split(" ",$t))[3]."\tStarting preprocessing\n\n");
print ((split(" ",$t))[3]."\tStarting preprocessing\n\n");

open(my $R1, "-|","gzip -cd $ARGV[0] | sed -n 2~4p | cut -c 3-27 " ) || die "cannot open read 1 file\n"; 

while (<$R1>) {
	$H_of_BC_counts{$_}++;
}
$counter = $.;	
close $R1;
$t = localtime;
print STDERR ((split(" ",$t))[3]."\tCollected $counter R1 reads\n\n");

	

####### pick tile info for each of the barcodes found in R1 (and discard barcodes mapping to different tiles)

foreach my $file (@barcode_files) {   
	chomp $file;	
	$tile = basename($file, ".idx");   
	open(my $bc_fh, "<", "$file") || die "cannot open $file\n";	
	$found=0;  
	while ($bc = <$bc_fh>) {  
		if (exists($H_of_BC_counts{$bc})) {  
			$found ++;  
		}
	}
	$all=$.;   
	close $bc_fh; 
	$ratio = $found/$all; 
		
	if ($ratio > $ARGV[5]) {
		$H_of_tiles{$file}="$all\t$found\t$ratio\t$tile\n";   #### use the tile name to save in a hash the name of the retained tile
	}
}

open(my $puck_stats, ">","$ARGV[2]/barcodes_per_tile.txt");
while ((my $t,my $stat) = each %H_of_tiles) {
	print $puck_stats $t."\t".$stat;
}
close $puck_stats;

my $num_tiles = keys(%H_of_tiles);
$t = localtime;
print STDERR ((split(" ",$t))[3]."\tSelected $num_tiles tiles.Now starting to check for ambigous barcodes\n");
print ((split(" ",$t))[3]."\tSelected $num_tiles tiles.Now starting to check for ambigous barcodes\n");

foreach my $file (keys %H_of_tiles) {
	$tile = basename($file, ".idx");   
	open(my $tile_fh, "-|", "zcat $ARGV[4]/$tile.txt.gz ") || die "cannot open $tile.txt.gz\n";
	$barcode_linein = <$tile_fh>;
	while ($barcode_linein = <$tile_fh>) {
		chomp $barcode_linein;   
		($bc,$x,$y) = split("\t",$barcode_linein);
		if (exists $H_of_BC_counts{"$bc\n"}) {
			if ((!exists $H_of_A_of_amb_bc{"$bc\n"})  && (!exists $H_of_BC{"$bc\n"})) {
				$H_of_BC{"$bc\n"}=$tile.",".$x.",".$y;
			}
			else {
				$H_of_A_of_amb_bc{"$bc\n"}++;
			}
		}
	}	
}	

my $ambiguous_count = keys %H_of_A_of_amb_bc;

$t = localtime;
print STDERR ((split(" ",$t))[3]."\tRemoved $ambiguous_count ambiguous barcodes.\n");
print ((split(" ",$t))[3]."\tRemoved $ambiguous_count ambiguous barcodes.\n");



open (my $exchange, ">", "$ARGV[2]/bc_info.tmp");

while ((my $bc, my $info) = each %H_of_BC) {
	chomp $bc;
	print $exchange "$bc\t$info\n";
}
close $exchange;

exit;


