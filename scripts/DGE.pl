use warnings;
use strict;
use File::Basename;



my $inline;
my @A_of_s;  
my $gene_list=" ";
my $non_zero=0;
my @map;
$map[ord("A")]="00";
$map[ord("C")]="01";
$map[ord("G")]="10";
$map[ord("T")]="11";
my %pam=("00"=>"A","01"=>"C","10"=>"G","11"=>"T");
my $N = 0;
my $M = 0;
my @cyto_scores;
my @nuc_scores;
my %genes;
my $tile = $ARGV[1];
my $umi_len;
my $offset_x;
my $offset_y;
my $tile_w;
if ($ARGV[3] eq "Illumina") {
	$tile_w=3000;
	$umi_len=9;
	my $swath= substr($tile,3,1);
	my $tile_id = substr($tile,4,2);
	#$tile_id =~ s/^0//;
	$offset_x= int(($swath*2164.421)+0.5);
	$offset_y = int((2327.559 * $tile_id)+0.5);
	if (($swath % 2) == 0) {
		if (substr($tile,2,1) eq "1") {$offset_y -= 387;}
		else {$offset_y += 387;}
	}
}

if ($ARGV[3] eq "Stereo-seq") {
	$tile_w=8000;
	$umi_len=10;
	$offset_x=substr($tile,0,2)*8000;
	$offset_y=substr($tile,2,2)*8000;
	}

my $umi_shift = $umi_len+1;


open(my $in, "<","$ARGV[0]/$tile") || die "cannot open $tile fifo\n$!";
my %counts;
my %umis;
while ($inline = <$in>) {
	chomp $inline;
	my ($gene,$loc,$umi,$coord) = split("\t",$inline);
	if (! exists $genes{$gene}) { $genes{$gene} = ++$N;}
	$counts{$gene}+=1;
	if (exists($A_of_s[$coord])) {
		my $idx = index($A_of_s[$coord], "\t$gene ");
		if ($idx == -1) {
			$non_zero++;
			$A_of_s[$coord] .= "$gene ${umi}$loc\t";
			$umis{$gene}+=1;
		}
		else {
			$idx += length($gene)+2;
			my $info = substr($A_of_s[$coord],$idx,(index($A_of_s[$coord], "\t", $idx)-($idx)));
			my $position=0;
			foreach my $rec_umi (unpack("(A$umi_shift)*",$info)) {
				if ($rec_umi eq $umi."N") { 
					substr($info,$position+$umi_len,1)=$loc;
					last;
				}
				if ($rec_umi eq $umi."Y") {
					if ($loc eq "U") {
						substr($info,$position+$umi_len,1)=$loc;
					}
				last;
				}
				if ($rec_umi eq $umi."U") {
				last;
				}
			$position+=$umi_shift;
			}				
			if ($position == length($info)) {
				$info = $umi.$loc.$info;
				$umis{$gene}+=1;
			}
			substr($A_of_s[$coord],$idx,(index($A_of_s[$coord], "\t", $idx)-($idx))) = $info;
		}
	}
	else {
		$M ++;
		$non_zero++;
		$A_of_s[$coord] = "\t$gene ".$umi.$loc."\t";
		$umis{$gene}+=1;
	}
}
my $t = (split(" ",localtime))[3];
print "$t\tCollected ${.} reads for tile $tile\n";
### writes the sparse matrix stored in %H_of_H to txt matrixmarket format
open (my $spateo, ">", "$ARGV[0]/dge/$tile.spateo.txt");
print $spateo "x\ty\tgeneID\tMIDCounts\tEXONIC\tINTRONIC\n";

open(my $out_fh,">", "$ARGV[0]/dge/dge.$tile.txt");
my $bc_seq=0;
my @barcodes_list;
my @scores;

my $xmin=1000000;
my $xmax=0;
my $ymin=1000000;
my $ymax=0;
print $out_fh "%\%MatrixMarket  matrix coordinate integer general\n$M $N $non_zero\n";
foreach my $coord (0..$#A_of_s) {
	if ($A_of_s[$coord]) {
		$bc_seq++;
		$barcodes_list[$bc_seq]=$coord;
		substr($A_of_s[$coord],0,1)="";
		my $bc_count=0;
		my $cyto_score= 0;
		my $nuc_score=0;
		my $x = ($coord % $tile_w)+$offset_x;
		if ($xmin > $x) {$xmin=$x;}
		if ($xmax < $x) {$xmax=$x;}
		my $y=int($coord/$tile_w)+$offset_y;
		if ($ymin > $y) {$ymin=$y;}
		if ($ymax < $y) {$ymax=$y;}
		foreach my $info (split("\t",$A_of_s[$coord])) {
			my @data = split(" ",$info);
			my $cyto_score = 0;
			my $nuc_score= 0;
			my $gene_idx = $genes{$data[0]};
			my $intronic = $data[1] =~ tr/U/N/;
			my $spliced = $data[1] =~ tr/Y/N/;
			my $num = $data[1] =~ tr/N/N/;
			my $exonic = $num - $intronic;
			print $out_fh "$bc_seq $gene_idx $num\n";
			print $spateo "$x\t$y\t$data[0]\t$num\t$exonic\t$intronic\n";
			$bc_count += $num;
			$cyto_score +=$spliced;
			$nuc_score += $intronic;
		}
		if ($bc_count > 4) {
			$cyto_score = $cyto_score/$bc_count;
			$nuc_score = $nuc_score/$bc_count;
			if ($cyto_score > 0) {push @cyto_scores, $cyto_score;}
			if ($nuc_score > 0) {push @nuc_scores, $nuc_score;}
			push @scores, "$bc_seq\t$cyto_score\t$nuc_score";
		}
	}		
}			
	

close $out_fh;
close $spateo;

### save ordered lists of genes and barcodes per tile (to be used when opening the sparse matrix file).
my @genes;
while (my ($k, $n) = each %genes) {
	$genes[$n]=$k;
}
open(my $out_gene_list,">","$ARGV[0]/dge/$tile.gene_list");
open(my $out_gene_stats,">","$ARGV[0]/dge/$tile.counts_stats");
my $tile_counts=0;
my $tile_umis=0;
foreach my $n (1..$#genes) {
	print $out_gene_list "$genes[$n]\n";
	print $out_gene_stats "$genes[$n]\t$counts{$genes[$n]}\t$umis{$genes[$n]}\n";
	$tile_counts += $counts{$genes[$n]};
	$tile_umis += $umis{$genes[$n]};
	}
close $out_gene_list;
open(my $tile_log, ">", "$ARGV[0]/logs/$tile.tile.log");
print $tile_log "UMIs:\t$tile_umis\nCOUNTS\t$tile_counts\nCOORDS:\t$xmin\t$xmax\t$ymin\t$ymax\n";
close $tile_log;
$t = (split(" ",localtime))[3];
print "$t\tcollected $tile_umis UMIs ($tile_counts counts) for tile $tile\n";

open(my $out_bc_list,">","$ARGV[0]/dge/$tile.bc_list");
foreach my $b (@barcodes_list[1..$bc_seq])  {print $out_bc_list "$b\n";}
close $out_bc_list;

unlink "$ARGV[0]/$tile";
		
### prepare nucleo/cytoplasm data
my @cyto_bc;
my @nuc_bc;
my @tile_nuc_counts;
my @tile_cyto_counts;
if ($ARGV[2] ne "apex") {
	my @sorted = sort { $b <=> $a } @cyto_scores;
	my $cyto_threshold=$sorted[int($#sorted*0.25)];
	@sorted = sort { $b <=> $a } values @nuc_scores;
	my $nuc_threshold=$sorted[int($#sorted*0.25)];
	
	foreach my $item (@scores) {
		my ($id, $nuc_score, $cyto_score) =split("\t",$item);
		if ($cyto_score > $cyto_threshold) {
			push  @cyto_bc, $id;
			}
		if ($nuc_score > $nuc_threshold) {
			push @nuc_bc, $id;
		}
	}
	open(my $MM_in, "<","$ARGV[0]/dge/dge.$tile.txt");
		my $header= <$MM_in>;
		$header= <$MM_in>;
		my $nuc = shift @nuc_bc;
		while (my $entry=<$MM_in>) {
			my ($bc_id, $gene_id, $counts)=split(" ", $entry);
			while ($nuc < $bc_id) {$nuc = shift @nuc_bc;}
			if ($nuc == $bc_id) {$tile_nuc_counts[$gene_id] += $counts;}
		}
	close $MM_in;
	open (my $nuclear, ">", "$ARGV[0]/dge/$tile.nuclear");
	my $i =1;
	foreach my $gene (keys %genes) {
		if (exists $tile_nuc_counts[$i]) {
			print $nuclear "$gene\t$tile_nuc_counts[$i]\n";
		}
		$i++;
	}
	close $nuclear;
	
	open($MM_in, "<","$ARGV[0]/dge/dge.$tile.txt");
		$header= <$MM_in>;
		$header= <$MM_in>;
		my $cyto = shift @cyto_bc;
		while (my $entry=<$MM_in>) {
			my ($bc_id, $gene_id, $counts)=split(" ", $entry);
			while ($cyto < $bc_id) {$cyto = shift @cyto_bc;}
			if ($cyto == $bc_id) {$tile_cyto_counts[$gene_id] += $counts;}
		}
	close $MM_in;
	open (my $cytoplasmic, ">", "$ARGV[0]/dge/$tile.cyto");
	$i =1;
	foreach my $gene (keys %genes) {
		if (exists $tile_cyto_counts[$i]) {
			print $cytoplasmic "$gene\t$tile_cyto_counts[$i]\n";
		}
		$i++;
	}
	close $cytoplasmic;
}
exit;
