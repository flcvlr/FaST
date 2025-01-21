use warnings;
use strict;

my %nuc = ("AA" => 0,"AC" => 1,"AG" => 2,"AT" => 3,"CA" => 4,"CC" => 5,"CG" => 6,"CT" => 7,"GA" => 8,"GC" => 9,"GG" => 10,"GT" => 11,"TA" => 12,"TC" => 13,"TG" => 14,"TT" => 15);
my @fhs;
my @output;
while (my ($i,$v) = each(%nuc)) {
	open($fhs[$v], "+>", "$ARGV[1]/${i}_barcodes") || die "cannot write to collector\n$!";
	$output[$v] = "";
	$fhs[$v]->autoflush(1);
	}

my $in;
my $index;
if ($ARGV[2] eq "Illumina") {
	open($in, "-|", "pigz -cd $ARGV[0] | sed -n 2~4p | cut -c 3-27 | fgrep -v N") || die "cannot open input\n$!";
}

if ($ARGV[2] eq "Stereo-seq") {
	open($in, "-|", "pigz -cd $ARGV[0] | sed -n 2~4p | cut -c 1-25 | fgrep -v N") || die "cannot open input\n$!";
}


while (my $line=<$in>) {
	$index =$nuc{substr($line,0,2,"")};
	$output[$index] .= $line;
	if (length($output[$index]) > 1000) {
		print {$fhs[$index]} $output[$index];
		$output[$index]="";
	}	
}
my $collected=$.;
my $t = (split(" ",localtime))[3];
print STDERR "$t\tcollected $collected barcodes from file $ARGV[0]\n";
system("echo 'collected $collected barcodes from file $ARGV[0]' >> $ARGV[1]/logs/run.log");

foreach my $i (0..15) {
	if ($output[$i] ne "") {print {$fhs[$i]} $output[$i];}
	close $fhs[$i];
}



close $in;

exit;
