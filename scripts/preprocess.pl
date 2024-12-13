use warnings;
use strict;

$ARGV[0] =~ y/,/ /;
$ARGV[1] =~ y/,/ /;

my $R1;
if ($ARGV[3] eq "Illumina") {
open($R1, "-|","pigz -cd $ARGV[0] | sed -n 2~4p | cut -c 3-27 " ) || die "cannot open read R1 file\n";
}

if ($ARGV[3] eq "Stereo-seq") {
open($R1, "-|","pigz -cd $ARGV[0] | sed -n 2~4p " ) || die "cannot open read R1 file\n";
}


open(my $R2, "-|","pigz -cd $ARGV[1]  " ) || die "cannot open read R2 file\n";

my @reads_fh;
my @buffer;
my @nuc = ("A","C","G" ,"T");
foreach my $N (@nuc) {
	open ($reads_fh[ord($N)], ">>","$ARGV[2]/${N}_reads") || die "cannot open reads fifo\n$!";
	$reads_fh[ord($N)]->autoflush(1);
	$buffer[ord($N)] = "";
	}

my $score;
my $umi;
my $seq;

my $bc;
my $start;
if ($ARGV[3] eq "Stereo-seq") {
	while ($bc= <$R1>) {
		my $name = <$R2>;
		$seq= <$R2>;
		<$R2>;
		$score=<$R2>;
		$umi=substr($bc,25,10,"");
		if (index($umi.$bc, "N") != -1) { next;}
		chomp $seq;
		chomp $score;
		chomp $name;
		substr($name,0,1,"");
		my $N=ord(substr($bc,0,1));
		$buffer[$N] .= "$name\t4\t*\t0\t0\t*\t*\t0\t0\t$seq\t$score\tMI:Z:$umi\tCB:Z:$bc";
		if (length($buffer[$N]) > 2000 ) {
			print {$reads_fh[$N]} $buffer[$N];
			$buffer[$N] ="";
		}
	}
}
if ($ARGV[3] eq "Illumina") {
	while ($bc= <$R1>) {
		my $name = <$R2>;
		$seq= <$R2>;
		<$R2>;
		$score=<$R2>;
		$umi=substr($seq,0,9,"");
		if (index($umi.$bc, "N") != -1) { next;}
		substr($score,0,9)="";
		chomp $seq;
		chomp $score;
		$name=substr($name,1,(index($name, " ")-1));
		my $N=ord(substr($bc,0,1));
		$buffer[$N] .= "$name\t4\t*\t0\t0\t*\t*\t0\t0\t$seq\t$score\tMI:Z:$umi\tCB:Z:$bc";
		if (length($buffer[$N]) > 2000 ) {
			print {$reads_fh[$N]} $buffer[$N];
			$buffer[$N] ="";
		}
	}
}

foreach my $N (@nuc) {
	if ($buffer[ord($N)] ne "") {print {$reads_fh[ord($N)]} $buffer[ord($N)];}
	}

exit;
