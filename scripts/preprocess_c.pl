
use warnings;
use strict;
use File::Basename;
#$|=1;
my $barcode_linein;
my %H_of_BC;
my $x;
my $y;
my $bc;
my $tile;
my @data;
my $t;
my $info;
my $bam;
open(my $exchange,"<","$ARGV[2]/bc_info.tmp");

while (<$exchange>) {
	chomp;   
	($bc,$info) = split("\t",$_);
	$H_of_BC{$bc}=$info;	
}
	
close $exchange;

#### prepare STAR input
if ($ARGV[4] eq "bam") {
	open($bam, "|-", "samtools view -Sb - -@ 2 > $ARGV[2]/input.bam");
	}
else{
	open($bam, ">", "$ARGV[2]/input.bam");
	}
print $bam "\@HD	VN:1.6\n"; 

### open fastq files
open(my $R1_R2, "-|","$ARGV[3]/prepare_input.sh $ARGV[0] $ARGV[1] $ARGV[3]" ) || die "cannot open read 1 file\n";   

### process fastq reads and write to bam
my $counter = 0; 
while (my $R1_linein = <$R1_R2>) { 
	chomp $R1_linein;
	@data=split("\t",$R1_linein);  
	$bc=substr($data[11],5);
	if (exists $H_of_BC{$bc}) {
		($tile, $x,$y) = split(",",$H_of_BC{$bc}); 
		print $bam "$R1_linein\tCX:i:$x\tCY:i:$y\tRG:Z:$tile\n";  
		$counter ++;
	}
}
		
close $bam;
close $R1_R2;
$t = localtime;
print STDERR ((split(" ",$t))[3]."\tWritten $counter entries to input bam file.\n");
print ((split(" ",$t))[3]."\tWritten $counter entries to input bam file.\n");


### write tile stats to file

exit;

