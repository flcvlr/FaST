use warnings;
use strict;

########### this script will parse a bed format input. It expects alla alignements of the same read to appear in consecutive lines of the input. It will look for the overlapping gene in field #11 of the bed input. It will output a signle line for each read in the format GENENAME\tLOC\tBC\tUMI\tXCOORD\tYCOORD\tTILE. This output is to be piped to preprocess_2_4.pl

my @data;
my $read_name;
my $gene;
my $loc;
my $bc;
my $umi;
my $x;
my $y;
my $l= <STDIN>;
my %anno_data;
my @anno;

@data=split("\t",$l);
($loc,$x,$y,$bc,$umi) = split("\\|",$data[6]);
if ($data[10] =~ /(.*?)#_INTRON$/) {$data[10] = $1;$anno_data{"INTRON"}=1; }
$anno_data{$data[10]}=1;
$read_name=$data[3];
open(my $pipe, "|-", "LC_ALL=C sort | uniq -c | sed 's/^ \\+//;s/ /\t/'  | pigz > $ARGV[0]/dge/$ARGV[1].DGE_input.txt.gz"); 

while ($l =<STDIN>) {
	@data=split("\t",$l);
	if ($read_name eq $data[3]) {
		if ($data[10] =~ /(.*?)#_INTRON$/) {$data[10] = $1;$anno_data{"INTRON"}=1; }
		$anno_data{$data[10]}=1;
		}
	else {  
		@anno = keys %anno_data;
		if (($#anno == 0) && ($anno[0] ne ".")) {
			if ($anno[0] =~ /^M[Tt]-/) {$loc="CY";}
			print $pipe "$anno[0]\t$loc\t$bc\t$umi\t$x\t$y\n";
			}
		if ( ($#anno == 1) and  (exists $anno_data{"INTRON"}) && (!exists $anno_data{"."}) ) {
			delete $anno_data{"INTRON"}; 
			@anno = keys %anno_data; 
			print $pipe "$anno[0]\tNU\t$bc\t$umi\t$x\t$y\n";
		}
		undef %anno_data;
		($loc,$x,$y,$bc,$umi) = split("\\|",$data[6]);
		$read_name=$data[3];
		if ($data[10] =~ /(.*?)#_INTRON$/) {$data[10] = $1;$anno_data{"INTRON"}=1; }
		$anno_data{$data[10]}=1;
	}
}

@anno = keys %anno_data;
if (($#anno == 0) && ($anno[0] ne ".")) {
	print $pipe "$anno[0]\t$loc\t$bc\t$umi\t$x\t$y\n";
}
if (($#anno == 1) and (exists $anno_data{"INTRON"}) && (!exists $anno_data{"."}) ) {delete $anno_data{"INTRON"}; @anno = keys %anno_data; print $pipe "$anno[0]\tNU\t$bc\t$umi\t$x\t$y\n";}
close $pipe;
exit;
