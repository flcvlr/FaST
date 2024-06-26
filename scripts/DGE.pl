use warnings;
use strict;
use File::Basename;

##### input is the output preparse_bed.pl (piped then through sort | uniq -c) in the format:  Number\tGENENAME\tLOC\tBC\tUMI\tx_coord\ty_coord where Number is the # of reads found for each molecule of transcript.
##### this script will:
#1)  output a sparse matrix containing puck level counts for each gene (along with an ordered list of genes and Barcodes (one barcode = one puck)
#2)   print a baysor input file to be used for segmentation of cells.

my $inline;
my %H_o_H;  ### hash of hashes sparse matrix: $H_of_H{barcode}{gene}
my %H_o_H_intronic;  ### hash of hashes sparse matrix: $H_of_H_intronic{barcode}{gene}
my $exonic;
my $intronic;
my %gene_list;
my %bc_list;
my $num;
my $gene;
my %H_of_counts;
my $loc;
my $bc;
my $umi;
my $x;
my $y;
my $k;
my $v;
my $non_zero=0;
my $gene_in_tile=1;
my $bc_in_tile=1;
my %nuc_score;
my %cyto_score;

my @cyto_scores;
my @nuc_scores;

my $tile = basename($ARGV[1],".DGE_input.txt.gz");
my %tile_nuc_counts;
my %tile_cyto_counts;

### open bam tag summary file (#reads\tGENE\tBARCODE\tUMI\txcoord\tycoord) and adds the gene name and barcode to the %gene_list and %bc_list (if not found yet), then adds 1 count to $H_of_H{bc}{gene} sparse matrix and updates the %H_of_counts hash

open (my $input_DGE, "-|","gzip -cd $ARGV[1]")  || die "cannot open $ARGV[2]\n";

while ($inline = <$input_DGE>) {
	chomp $inline;
	($num,$gene,$loc,$bc,$umi,$x,$y) = split("\t",$inline);
	if ($loc eq "NU") {$nuc_score{$bc}++;$H_o_H_intronic{$bc}{$gene}++;}
	if ($loc eq "CY") { $cyto_score{$bc} ++;}
	if (!exists($gene_list{$gene})) {$gene_list{$gene} =$gene_in_tile;$gene_in_tile +=1} ### so each slot of %gene_list will contain a new integer
	if (!exists($bc_list{$bc})) {$bc_list{$bc} =$bc_in_tile;$bc_in_tile +=1}  ### so each slot of %bc_list will contain a new integer
	if (!exists($H_o_H{$bc}{$gene})) { $non_zero++;}
	$H_o_H{$bc}{$gene}++;
	$H_of_counts{$bc}{"counts"}++;
	$H_of_counts{$bc}{"reads"}+= $num; 
	$H_of_counts{$bc}{"coord"}="$x,$y";
	}

close $input_DGE;
### writes the sparse matrix stored in %H_of_H to txt matrixmarket format

open(my $out_fh,">", "$ARGV[0]/dge/dge.$tile.txt");

my $N = keys %gene_list;
my $M = keys %bc_list;
print $out_fh "%\%MatrixMarket  matrix coordinate integer general\n$M $N $non_zero\n";
foreach $bc (keys %H_o_H) {
	foreach $gene (keys %{$H_o_H{$bc}}) {
		print $out_fh "$bc_list{$bc} $gene_list{$gene} $H_o_H{$bc}{$gene}\n";
	}
}
close $out_fh;

### save ordered lists of genes and barcodes per tile (to be used when opening the sparse matrix file).

open(my $out_gene_list,">","$ARGV[0]/dge/$tile.gene_list");
foreach  $gene (sort { $gene_list{$a} <=> $gene_list{$b} } keys %gene_list) {print $out_gene_list "$gene\n";}
close $out_gene_list;

open(my $out_bc_list,">","$ARGV[0]/dge/$tile.bc_list");
foreach  $bc (sort { $bc_list{$a} <=> $bc_list{$b} } keys %bc_list) {print $out_bc_list "$bc\n";}
close $out_bc_list;

### save spateo input

open (my $spateo, ">", "$ARGV[0]/dge/$tile.spateo.txt");
print $spateo "x\ty\tgeneID\tMIDCounts\tEXONIC\tINTRONIC\n";
foreach $bc (keys %H_o_H) {
	foreach $gene (keys %{$H_o_H{$bc}}) {
		if (exists $H_o_H_intronic{$bc}{$gene}) {
			$exonic= $H_o_H{$bc}{$gene} - $H_o_H_intronic{$bc}{$gene}; 
			$intronic = $H_o_H_intronic{$bc}{$gene};
		}
		else {
			$exonic = $H_o_H{$bc}{$gene};$intronic=0;
		}
		($x,$y)=split(",",$H_of_counts{$bc}{'coord'});
		print $spateo "$x\t$y\t$gene\t$H_o_H{$bc}{$gene}\t$exonic\t$intronic\n";
	}
	if ($H_of_counts{$bc}{'counts'} > 4) {
		if (!exists $nuc_score{$bc}) {push  @nuc_scores, 0;} else {push @nuc_scores, $nuc_score{$bc}/$H_of_counts{$bc}{'counts'};}
		if (!exists $cyto_score{$bc}) {push  @cyto_scores, 0;} else {push @cyto_scores, $cyto_score{$bc}/$H_of_counts{$bc}{'counts'};}
	}
}

close $spateo;
		
### prepare nucleo/cytoplasm data

if ($ARGV[2] ne "apex") {
	my @sorted = sort { $b <=> $a } @cyto_scores;
	my $cyto_threshold=$sorted[int($#sorted*0.25)];
	@sorted = sort { $b <=> $a } values @nuc_scores;
	my $nuc_threshold=$sorted[int($#sorted*0.25)];
	foreach $bc (keys(%H_o_H)) {
		if ($H_of_counts{$bc}{'counts'} < 5) {next;}
		if (exists $cyto_score{$bc}) {
			if ($cyto_score{$bc}/$H_of_counts{$bc}{'counts'} > $cyto_threshold) {
				foreach $gene (keys %{$H_o_H{$bc}}) {
					$tile_cyto_counts{$gene} += $H_o_H{$bc}{$gene}
				}
			}
		}
		if (exists $nuc_score{$bc}) {
			if ($nuc_score{$bc}/$H_of_counts{$bc}{'counts'} > $nuc_threshold) {
				foreach $gene (keys %{$H_o_H{$bc}}) {
					$tile_nuc_counts{$gene} += $H_o_H{$bc}{$gene}
				}
			}
		}
	}

	open (my $nuclear, ">", "$ARGV[0]/$ARGV[1]/complete_data/dge/$tile.nuclear");
	while (($k, $v)=each %tile_nuc_counts) {
		print $nuclear "$k\t$v\n";
	}
	close $nuclear;

	open (my $cytoplasmic, ">", "$ARGV[0]/$ARGV[1]/complete_data/dge/$tile.cytoplasmic");
	while (($k, $v)=each %tile_cyto_counts) {
		print $cytoplasmic "$k\t$v\n";
	}
	close $cytoplasmic;
}
exit;
