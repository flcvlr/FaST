use warnings;
use strict;

#### this script reads input from STAR in sam format, intersects alignments with annotation and reports a single entry for each read name (only if all alignments of that read align to a single gene). The script also adds a nuclear or cy

my @bam;
my @bed;
my @cigar;
my $c;
my $tags;
my $op;
my %H_of_tile_fh;
my @in_data;
my $bc;
my $out="";
my %H_of_BC_counts;
my $Phix=0;
my $rRNA=0;
my $references = "\t";
my %chroms;
my $chroms_index = 0;
my $genes_index =0;
my @maps;
my $binsize = 500;
my %genes;
my @genes_id;
my %output;
my $tile_offset;
my $tile_len;
my $umi_offset;
my $umi_len;

if ($ARGV[3] eq "Illumina") {
	$tile_offset=-21;
	$tile_len=6;
	$umi_offset=-67;
	$umi_len=9;
	}

if ($ARGV[3] eq "Stereo-seq") {
	$tile_offset=-19;
	$tile_len=4;
	$umi_offset=-66;
	$umi_len=10;
	}

open (my $bed, "<", "$ARGV[0]/GENE_annotation_intron_exon.bed");
while (my $l=<$bed>) {
	chomp $l;	
	my @ref_data=split("\t",$l);
	if (! exists $chroms{$ref_data[0]}) {$chroms{$ref_data[0]} = $chroms_index; $chroms_index++;}
	my $index=$chroms{$ref_data[0]};
	if (!exists($genes{$ref_data[3]})) {$genes{$ref_data[3]}=$genes_index;push @genes_id , $ref_data[3]; $genes_index +=1;}
	$ref_data[3] = $genes{$ref_data[3]};
	if ($ref_data[5] eq "-") {$ref_data[5] = 1;} else {$ref_data[5]=0;}
	foreach my $kb (int($ref_data[1]/$binsize)..int($ref_data[2]/$binsize)) {
		$maps[$index][$kb] .= pack("VVVA",$ref_data[1],$ref_data[2],$ref_data[3],$ref_data[5]);
	}
}

my @list;
my %cand;
my $strand;
my @tiles = `cut -f 5 $ARGV[0]/barcodes_per_tile.txt `;
my %tile_fh;

foreach my $n (@tiles) {
	chomp $n;
	open($tile_fh{$n} ,">>", "$ARGV[0]/$n"); 
	$tile_fh{$n}->autoflush(1);
}
my $target;
my $name="";
my %anno_data;
my @anno;
my $linein;
my $tile;

while ($linein =<STDIN>) {
	@bam=split("\t",$linein,7);
	if ($bam[0] ne $name ) { 
		@anno=keys %anno_data;
		if ($#anno == 0) {
			my $incipit =substr($anno[0],0,3);
			if (($incipit eq "MT-") || ($incipit eq "Mt-") || ($incipit eq "mt-") ) {substr($tags,0,1)="Y";}
			if (!exists $anno_data{"rRNA"}) {
				$output{$tile} .=  "$anno[0]\t$tags\n";
				if (length($output{$tile}) > 2000) { print { $tile_fh{$tile} } $output{$tile};$output{$tile}="";}
			}
		}
		elsif ( ($#anno == 1) &&  (exists $anno_data{"INTRON"}) ) {
			delete $anno_data{"INTRON"}; 
			substr($tags,0,1) ="U";
			@anno = keys %anno_data; 
			$output{$tile} .=  "$anno[0]\t$tags\n";
			if (length($output{$tile}) > 2000) { print { $tile_fh{$tile} } $output{$tile};$output{$tile}="";}
		}
		elsif ( $#anno > 1) {
			if (!exists $anno_data{"rRNA"}) {
				$output{$tile} .="Multigene\t$tags\n";
				if (length($output{$tile}) > 2000) { print { $tile_fh{$tile} } $output{$tile};$output{$tile}="";}
			}
		}
		elsif ( ($#anno == -1) && ($name ne "")) {
			$output{$tile} .="Intergenic\t$tags\n";
			if (length($output{$tile}) > 2000) { print { $tile_fh{$tile} } $output{$tile};$output{$tile}="";}
		}
		%anno_data= ();
		$tile=substr($linein,$tile_offset,$tile_len);
		$tags="N\t".substr($linein,$umi_offset,$umi_len)."\t".substr($linein,-9,8);
		$name=$bam[0];
		}
	if (! exists($chroms{$bam[2]})) {			
		if (substr($bam[2],0,11) eq  "NC_001422.1") {$Phix ++;next;}
		if (substr($bam[2],0,6) eq  "URS000") {$rRNA ++;$anno_data{"rRNA"}=1;next;}
		next;
	}
	$strand = ($bam[1] & 16)/16;
	my $index=$chroms{$bam[2]};
	$bed[1]=$bam[3]-1;
	$bed[2]=$bed[1];
	@cigar = ($bam[5] =~ /(\d+[MIN])/g);
		foreach $c (@cigar) {
			$op = chop($c);
			if ($op eq "N") {
				substr($tags,0,1) ="Y";
				foreach my $kb (int($bed[1]/$binsize)..int($bed[2]/$binsize)) {
					if ($maps[$index][$kb]) {
					@list= unpack("(A13)*",$maps[$index][$kb]); 
					@cand{@list} = undef;
					}
				}
							
				foreach my $i (keys %cand) {
					my @cdata=unpack("VVVA",$i);
					if ( ($strand != $cdata[3]) || ($cdata[0] >= $bed[2]) || ($cdata[1] <= $bed[1])) {next;}
					@anno_data{split("#_",$genes_id[$cdata[2]])}=1;
				}
				%cand=();
				$bed[1]=$bed[2]+$c; $bed[2]=$bed[1];
				next;
			}	
			$bed[2]+=$c;
		}
	
		foreach my $kb (int($bed[1]/$binsize)..int($bed[2]/$binsize)) {
			if ($maps[$index][$kb]) {
			@list= unpack("(A13)*",$maps[$index][$kb]); 
			@cand{@list} = undef;
			}
		}			
	foreach my $i (keys %cand) {
		my @cdata=unpack("VVVA",$i);
		if (($strand != $cdata[3]) ||($cdata[0] >= $bed[2]) || ($cdata[1] <= $bed[1])) {next;}
		@anno_data{split("#_",$genes_id[$cdata[2]])}=1;
		
	}
	%cand=();
}

@anno = keys %anno_data;
if ($#anno == 0) {
	if ($anno[0] =~ /^[Mm][Tt]-/) {substr($tags,0,1)="Y";}
	$output{$tile} .= "$anno[0]\t$tags\n";
}
if ( ($#anno == 1) &&  (exists $anno_data{"INTRON"}) ) {
	delete $anno_data{"INTRON"}; 
	substr($tags,0,1) ="U";
	@anno = keys %anno_data; 
	$output{$tile} .= "$anno[0]\t$tags\n";
}

while (my($k, $v) =each %output) {
	print {$tile_fh{$k}} $v;
}

my $t = (split(" ",localtime))[3];
if (($rRNA + $Phix) > 0) {
	my $t = (split(" ",localtime))[3];
	my $thread =$ARGV[4]+1;
	print "$t\tThread $thread finished processing ${.} alignments.\tPhiX mapped reads: $Phix\trRNA mapped reads:\t$rRNA\n";
}

exit;




	
	
	
