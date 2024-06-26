use warnings;
use strict;

#### this script reads input from STAR in sam format, collects BC statistics and splits the output (after conversion to bed format) based on tile to multiple parallel processes.

my @bam;
my @bed;
my @cigar;
my $c;
my @tags;
my $op;
my %H_of_tile_fh;
my @in_data;
my $tile;
my $bc;
my $out;
my %H_of_BC_counts;

open (my $pipe ,"|-","intersectBed -a - -b $ARGV[0]/GENE_annotation_intron_exon.bed -s -wao | perl $ARGV[2]/preparse_bed.pl $ARGV[0] $ARGV[1]"); 
while (my $linein =<STDIN>) {
	if ($linein eq "END_OF_RECORDS\n") {last;}
	@bam=split("\t",$linein);
	if ($bam[2] !~/^chr/) {next;}	
	@tags = ($linein =~ /CX:i:([\d]+).+?CY:i:([\d]+).+?CB:Z:([ACTG]+).+?MI:Z:([ACTG]+)/g);
	if ($bam[5] =~/N/) {unshift @tags, "CY";} else {unshift @tags, "NA";}
	$bed[0]=$bam[2];
	$bed[1]=$bam[3]-1;
	$bed[3]=$bam[0];
	$bed[4]=$bam[4];
	$bed[2]=$bed[1];
	if ($bam[1] & 16) {$bed[5]="-";} else {$bed[5]="+";}
	@cigar = ($bam[5] =~ /(\d+[MIN])/g);
		foreach $c (@cigar) {
		$op = chop($c);
		if ($op eq "M") {$bed[2]+=$c;}
		if ($op eq "I") {$bed[2]+=$c;}
		if ($op eq "N") {$out=join("\t",@bed)."\t".join("|",@tags); print $pipe $out."\n"; 
		$bed[1]=$bed[2]+$c; $bed[2]=$bed[1];}
		}
	$out=join("\t",@bed)."\t".join("|",@tags);
	print $pipe $out."\n";
	}
close $pipe;
close STDIN;
exit;



	
	
	
