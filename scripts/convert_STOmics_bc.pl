use warnings;
use strict;
my %pam=("00"=>"A","01"=>"C","11"=>"G","10"=>"T");
open(my $in, "<", "$ARGV[0]");
my %tile_fh;
while (<$in>) {
	if ($_ =~ /^ {0,6}\(\d+.*/ ) {
		$_=~ y/ (),/    /;
		my @d=(split(" ",$_))[0,1,4];
		if ($d[2] != 0) {
			my $bc=join("",@pam{unpack("A2" x 25 ,sprintf ("%050b", $d[2]))});
			my $tile = sprintf("%02d", int($d[0]/8000)).sprintf("%02d",int($d[1]/8000));
			my $x = $d[0] % 8000;
			my $y = $d[1] % 8000;
			if (!exists $tile_fh{$tile}) {
				open($tile_fh{$tile}, ">", "$ARGV[1]/$tile.txt");
				}
			print {$tile_fh{$tile}} "$x\t$y\t$bc\n";
			}
		}
	}
close $in;
foreach my $tile (keys %tile_fh){
	close $tile_fh{$tile};
	}
exit;
