use warnings;
use strict;
use File::Basename;

my %map;
my %kmers;
my @barcode_files = `ls $ARGV[4]/*.idx`;
my %nuc = ("AA" => 0,"AC" => 1,"AG" => 2,"AT" => 3,"CA" => 4,"CC" => 5,"CG" => 6,"CT" => 7,"GA" => 8,"GC" => 9,"GG" => 10,"GT" => 11,"TA" => 12,"TC" => 13,"TG" => 14,"TT" => 15);
my @cun=reverse %nuc;
my @fhs;
my @fhs_fifo;
my @pids;

system("mkfifo $ARGV[2]/reports");

my $t = (split(" ",localtime))[3];
print "$t\tStarting preprocessing\n";

### launch parallel workers to build R1 database

while (my ($i,$v) = each(%nuc)) {
	system("mkfifo $ARGV[2]/${i}_barcodes");
	$pids[$v] = open($fhs[$v], "|-", "perl $ARGV[3]/collector_1.pl $i $ARGV[2] $ARGV[9]") || die "cannot open collector\n$!";
	$fhs[$v]->autoflush(1);
	open($fhs_fifo[$v], ">>", "$ARGV[2]/${i}_barcodes");
	}


### read R1 and pass data to workers

my @zips_1= split(",",$ARGV[0]);
my @unzipper_pid;
my @unzipper_fhs;
foreach my $k (0..$#zips_1) {
	$unzipper_pid[$k] = open($unzipper_fhs[$k], "-|", "perl $ARGV[3]/unzipper.pl $zips_1[$k] $ARGV[2] $ARGV[9]") || die "cannot open unzipper\n$!";
}
foreach my $id (0..$#zips_1) {
	close $unzipper_fhs[$id];
	waitpid($unzipper_pid[$id],0);
}


while (my ($i,$v) = each(%nuc)) {
	$fhs_fifo[$v]->autoflush(1);
	print {$fhs_fifo[$v]} "END OF RECORDS\n";
	close $fhs_fifo[$v];
	}
	

### read tile indexes, pass to workers and collect stats to identify tiles

my @all;
my $found;
my %H_of_tiles;
my $tile;
my $bc;
my $rep_fh;
open($rep_fh, "+<", "$ARGV[2]/reports") || die "cannot open reports\n$!";

system("mkfifo $ARGV[2]/A_barcodes");
system("mkfifo $ARGV[2]/C_barcodes");
system("mkfifo $ARGV[2]/G_barcodes");
system("mkfifo $ARGV[2]/T_barcodes");

my $start;
my $all;
foreach my $file (@barcode_files) {
	chomp $file;
	$tile = basename($file, ".idx");   
	open(my $bc_fh, "<", "$file") || die "cannot open $file\n";	
	$found=0;  
	while ($bc = <$bc_fh>) {
		$start=substr($bc,0,2,"");
		print {$fhs[$nuc{$start}]} $bc;
		}
		$all= $.;
	close $bc_fh; 
	
	for my $i (0..$#fhs) {
		print {$fhs[$i]} "END OF FILE\n";
		my $rep=readline($rep_fh);
		chomp $rep;
		$found += $rep;
		}
	my $ratio = $found/$all;
	if ($ratio > $ARGV[5]) {
		$H_of_tiles{$file}="$all\t$found\t$ratio\t$tile\n";  
	}
}

for my $i (0..$#fhs) {
	print {$fhs[$i]} "END OF RECORDS\n";
	$fhs[$i]->autoflush(0);
	}

$t = (split(" ",localtime))[3];
print "$t\tBased on the collected barcodes, the RNA capture was done using the following tiles:\n";	
open(my $puck_stats, ">","$ARGV[2]/barcodes_per_tile.txt");
while ((my $t,my $stat) = each %H_of_tiles) {
	my $tile_found=basename($t,".idx");
	print "\t\t\t$tile_found\n";
	print $puck_stats $t."\t".$stat;
	}
close $puck_stats;


###  focus on selected tiles
system("mkfifo $ARGV[2]/input.sam");

foreach my $file (keys %H_of_tiles) {
	$tile = basename($file, ".idx");
	open(my $tile_fh, "-|", "pigz -cd $ARGV[4]/$tile.txt.gz ") || die "cannot open file $ARGV[4]/$tile.txt.gz\n";
	foreach my $i (0..$#fhs) {print {$fhs[$i]} "NEW_TILE_\t$tile\t0\n";}
	while (my $line = <$tile_fh>) {
		$start=substr($line,0,2,"");
		print {$fhs[$nuc{$start}]} $line;
	}
}


my $num_tiles = keys %H_of_tiles;
my $STARpid = open(my $merge, "-|", " STAR --runThreadN $ARGV[6] --genomeDir $ARGV[7]/reference/$ARGV[8]/$ARGV[8].star --clip3pAdapterSeq polyA --readFilesType SAM SE --readFilesIn $ARGV[2]/input.sam --outFileNamePrefix  $ARGV[2]/Aligned.$ARGV[8]. --outStd SAM --outTmpDir $ARGV[2]/star_tmp | samtools view -Sb - -1 -@ 2 > $ARGV[2]/Aligned.bam");


my $ambiguous=0;
my $rep;
my $ret;
my $retained=0;
while( my ($k, $i) = each %nuc) {
	$fhs[$i]->autoflush(1);	
	print {$fhs[$i]} "END OF RECORDS\n";
	$fhs[$i]->autoflush(0);
	($ret,$rep)=split("\t",readline($rep_fh));
	chomp $rep;
	$ambiguous +=$rep;
	$retained += $ret;
	}
	
my @filter_fh;
my @filter_pid;
my @reads_fh;

open (my $in_sam, ">>", "$ARGV[2]/input.sam");
foreach my $file (keys %H_of_tiles) {
	$tile = basename($file, ".idx");
	print $in_sam "\@RG\tID:$tile\n";
	}


my %nuc_bc = ("A" => 0,"C" => 1,"G" => 2,"T" => 3);

while( my($N, $i) = each %nuc_bc) {
	system("mkfifo $ARGV[2]/${N}_reads");
	$filter_pid[ord($N)] = open ($filter_fh[ord($N)], "|-","perl $ARGV[3]/filter.pl $N $ARGV[2] $ARGV[9]") || die "cannot launch filter.pl\n$!";
	}
	
my $percent;
$t =(split(" ",localtime))[3];
if ($ARGV[9] eq "Illumina") {
	$percent= sprintf("%.2f" , 100*($ambiguous/$retained));
	print "$t\tRetained $retained different barcodes. Found $ambiguous ambiguous barcodes ($percent%).\n$t\tStarting alignment with STAR...\n";
}

if ($ARGV[9] eq "Stereo-seq") {
	print "$t\tRetained $retained different barcodes.\n$t\tStarting alignment with STAR...\n";
}


while (my ($i,$v) = each(%nuc)) {
	close $fhs[$v];
	waitpid $pids[$v],0;
}
	

my @zips_2= split(",",$ARGV[1]);
my @preprocesser_pid;
my @preprocesser_fhs;
foreach my $k (0..$#zips_1) {
 	$preprocesser_pid[$k] = open($preprocesser_fhs[$k], "-|", "perl $ARGV[3]/preprocess.pl $zips_1[$k] $zips_2[$k] $ARGV[2] $ARGV[9]") || die "cannot open preprocesser\n$!";
}
foreach my $id (0..$#zips_1) {
	close $preprocesser_fhs[$id];
	waitpid($preprocesser_pid[$id],0);
}

foreach my $k (keys %nuc_bc) {
	close $filter_fh[ord($k)];
	waitpid($filter_pid[ord($k)],0)
}
close $in_sam;
close $rep_fh;
close $merge;

unlink "$ARGV[2]/reports";
while (my $file= glob("$ARGV[2]/*_barcodes") ) {unlink $file;}
while (my $file= glob("$ARGV[2]/*_reads") ) {unlink $file;}
unlink "$ARGV[2]/input.sam";
#print "passed $aligned records to STAR for genome mapping\n";
$t =(split(" ",localtime))[3];
print "$t\tFinished alignment with STAR, now building DGE matrix\n";
exit;


