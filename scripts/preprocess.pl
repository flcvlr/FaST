use warnings;
use strict;
use File::Basename;

### ARGV: 0=R1_reads; 1=R2_reads; 2=sample; 3=script_directory; 4=tiles_map_directory 5=HD_recover(0/1)

### requires samtools, zcat, sed, cut, grep 
my $counter=0;  ### num of R1 records processed
my $retained=0;  ### num of retained reads ( w/ non ambiguous barcodes in "abundant" tiles)
my $discarded=0; ### discarded barcodes
my @barcode_files = `ls $ARGV[4]/*`; 
my $ambiguous=0; ### number of ambiguous barcodes found
my $unmapped;  ### number of unmapped barcodes
my $tile; 
my %H_of_A_of_amb_bc;
my %H_of_BC;  ### main hash (dictionary) containing barcode info 
my %H_of_BC_counts;  ### hash containing barcode counts (for spacemake-like statistics output)
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


$t = localtime;
print STDERR ((split(" ",$t))[3]."\tStarting preprocessing\n\n");


open(my $R1, "-|","pigz -cd $ARGV[0] | sed -n 2~4p | cut -c3-27 | fgrep -v N " ) || die "cannot open read 1 file\n"; 

while ( <$R1>) {  ##### read one line of input
	chomp;    ##### remove newline at the end of line
	$H_of_BC_counts{$_}++;   #### increment by one (or set if not set yet) the counter for that barcode. This will keep track of all barcodes
	}
$counter = $.;	#### before closing the file, copy from the special variable $. the number of lines read
close $R1;	

$t = localtime;
print STDERR ((split(" ",$t))[3]."\tRead $counter barcodes. Now starting to check for ambigous barcodes\n");

####### pick tile info for each of the barcodes found in R1 (and discard barcodes mapping to different tiles)

open(my $puck_stats, ">","$ARGV[2]/barcodes_per_tile.txt");   
foreach my $file (@barcode_files) {   
	chomp $file;	
	$tile = basename($file, ".txt.gz");   
	open(my $bc_fh, "-|", "zcat $file | sed 1,1d ") || die "cannot open $file\n";	
	$found=0;  
	
	while ($barcode_linein = <$bc_fh>) {  
		chomp $barcode_linein;   
		($bc,$x,$y) = split("\t",$barcode_linein);  
		$x=int($x/9);
		$y=int($y/9);
		if (exists($H_of_A_of_amb_bc{$bc})) {
			push @{$H_of_A_of_amb_bc{$bc}},join("\t",$tile,$x,$y);
			next;
		}
		if (exists($H_of_BC_counts{$bc})) {  
			if (!exists($H_of_BC{$bc})) {  
				$H_of_BC{$bc}=join("\t",$tile,$x,$y);
				$found ++;  
				$matched += $H_of_BC_counts{$bc};
			}
			else { ### if the same barcode was already met
				push @{$H_of_A_of_amb_bc{$bc}}, $H_of_BC{$bc};
				push @{$H_of_A_of_amb_bc{$bc}}, join("\t",$tile,$x,$y);   ### label as ambiguous
				delete $H_of_BC{$bc};
			}
		}				
	}

	$all=$.;   #### again, before closing the file (and resetting $.) save its value in a variable
	close $bc_fh;  ### close the tile info file
	$ratio = $found/$all;  ### pretty obvious....
	if ($ratio > 0.1) {    #### retain only the tiles in which we found more than 10% of barcodes
		print $puck_stats "$tile\t$all\t$found\t$ratio\t$tile\n";   #### print a line with tile info to puck count file
		$H_of_tiles{$tile}=1;   #### use the tile name to save in a hash the name of the retained tile
	}
}

close $puck_stats;
	
#### clean up: remove ambiguous barcodes from tile info and bc counts

### remove tile info for barcodes outside of selected tile (but keep the barcodes in %H_of_BC_counts to search for a HD==1 match in the selected tiles)

foreach my $a (keys %H_of_BC) {
	$tile = (split("\t",$H_of_BC{$a}))[0];
	if (!exists $H_of_tiles{$tile}) {
		delete $H_of_BC{$a};
	}
}
my $num_bc_in_tiles= keys %H_of_BC;

my $ambiguous_count = keys %H_of_A_of_amb_bc;
$t = localtime;
print STDERR ((split(" ",$t))[3]."\tFound $ambiguous_count ambiguous barcodes. Now starting to fix ambigous barcodes prioritizing selected tiles.\n");
### fix ambiguous barcodes by picking the most likely tile (if only one)

foreach $bc (keys %H_of_A_of_amb_bc) {
	$candidate="";
	foreach my $info (@{$H_of_A_of_amb_bc{$bc}}) {
		$tile=(split("\t",$info))[0];
		if (exists $H_of_tiles{$tile}) {
			if ($candidate ne "") {$candidate =""; last;}
			else {$candidate=$info;}
		}
	}
	if ($candidate ne "") {$H_of_BC{$bc}=$candidate;}
	if ($candidate eq "") {
		$discarded +=$H_of_BC_counts{$bc};
		delete $H_of_BC_counts{$bc};
		$ambiguous++;
	}
}
my $amb_discarded=$discarded;
undef %H_of_A_of_amb_bc;

my $num_tiles= keys %H_of_tiles;  
my $size = keys %H_of_BC_counts;   
my $counter1=0;
my $unrecovered_because_amb=0;
$t = localtime;
print ((split(" ",$t))[3]."\tCollected $size different barcodes from $counter R1 records.\nSelected $num_bc_in_tiles (corresponding to $matched out of $size) in $num_tiles tiles.\n");


if ($ARGV[5] eq "1") {
print "Now looking for HD==1 matches to the remaining barcodes in the $num_tiles selected tiles..\n";   #### print short output to terminal
foreach my $a (keys %H_of_BC_counts) {   #### iterate over the barcodes we found in the initial round
	if (!exists $H_of_BC{$a}) {  
		$counter1++; 
		$candidate ="";
		for my $i (0..24) {
			$bc=$a;
			substr($bc,$i,1)= "A"; if (exists($H_of_BC{$bc})) {
				if ($candidate ne "") {$candidate =""; $unrecovered_because_amb++; last;}
				else {$candidate = $bc;}
				}
			$bc=$a;
			substr($bc,$i,1)= "C"; if (exists($H_of_BC{$bc})) {
				if ($candidate ne "") {$candidate =""; $unrecovered_because_amb++; last;}
				else {$candidate = $bc;}
				}
			$bc=$a;
			substr($bc,$i,1)= "G"; if (exists($H_of_BC{$bc})) {
				if ($candidate ne "") {$candidate =""; $unrecovered_because_amb++; last;}
				else {$candidate = $bc;}
				}
			$bc=$a;
			substr($bc,$i,1)= "T"; if (exists($H_of_BC{$bc})) {
				if ($candidate ne "") {$candidate =""; $unrecovered_because_amb++; last;}
				else {$candidate = $bc;}
				}
			}
		
		if ($candidate ne "") {
			$H_of_corrected_bc{$a}=$candidate;
			$H_of_BC_counts{$candidate}+=$H_of_BC_counts{$a};
			$HD_recovered_entries += $H_of_BC_counts{$a};
			delete $H_of_BC_counts{$a};
			$HD_recovered_bc++;
#			if ($counter1 % 1000000==0) {
#				print "\rtested $counter1 unmapped bc, remapped $HD_recovered_bc;";
#			}
#			
		}
		if ($candidate eq "") {
			$discarded +=$H_of_BC_counts{$a};  #### increase the number of reads discarded by the number of times we found that orphan barcode
			delete $H_of_BC_counts{$a};   #### remove it
			$unmapped ++;   #### increment by 1 the number of unmapped (orphan) barcodes
			next;
		}
	}
}
}
else {$unmapped=$size - $num_bc_in_tiles;}


$size = keys %H_of_BC;   #### collect statistics (as above, we use here the "keys" funtion in scalar context)
$retained=$counter - $discarded;   #### obvious calculations

#$t = localtime;
#print ((split(" ",$t))[3]."\tFound $size barcodes in $num_tiles tiles;\n $discarded entries were discarded as corresponding to $ambiguous ambiguous ($amb_discarded)  and $unmapped unmapped barcodes;\nretained $retained records (including $HD_recovered_bc unmapped barcodes (corresponding to $HD_recovered_entries entries) with HD == 1 from mapped barcodes) out of $counter\n $unrecovered_because_amb bc could not be remapped because there were > 1 possible valid bc at HD==1\n");  #### print message to screen

#### prepare bam output
open(my $bam, "|-", "sed s/^\@// | samtools view -Sb -@ 2 -1 - > $ARGV[2]/input.bam");
print $bam "\@\@HD	VN:1.6\n"; 

#### add valid tiles as read groups
foreach $tile (keys %H_of_tiles) {print $bam "\@\@RG\tID:$tile\n";} 

### open fastq files
open(my $R1_R2, "-|","$ARGV[3]/prepare_input.sh $ARGV[0] $ARGV[1] $ARGV[3]" ) || die "cannot open read 1 file\n";   

### process fastq reads and write to bam, collecting Barcodes sequences into hash
$counter = 0; 
while ($R1_linein = <$R1_R2>) { 
	chomp $R1_linein;
	@data=split("\t",$R1_linein);  
	if (exists $H_of_corrected_bc{$data[0]}) {$data[0]=$H_of_corrected_bc{$data[0]};}
	if (exists $H_of_BC{$data[0]}) {
		($tile, $x,$y) =split("\t",$H_of_BC{$data[0]}); 
		print $bam "$data[1]\t4\t*\t0\t0\t*\t*\t0\t0\t$data[3]\t$data[4]\tCX:i:$x\tCY:i:$y\tCB:Z:$data[0]\tMI:Z:$data[2]\tRG:Z:$tile\n";  
		$counter ++;
	}
}
		
close $bam;
close $R1_R2;
$t = localtime;
print STDERR ((split(" ",$t))[3]."\tWritten $counter entries to input bam file. The number is slightly lower than retained entries in previous step because entries with Ns in the UMI were discarded in this latter step\n");

### write tile stats to file

exit;

