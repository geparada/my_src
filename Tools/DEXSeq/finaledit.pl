use strict;
open(FILE, $ARGV[0]) or die "I did not read any file";
while(<FILE>){
#    chomp($_);
#	chomp($_);
	$_ =~ s///g;
	$_ =~ s/\r//g;
	if($_ =~ /(\S+)\t(\S+)\t(exon)\t(\d+)\t(\d+)\t\t(\+|\-)\t\t(.+)/){
	#    $7 =~ s///g;
		print "$1\t$2\t$3\t$4\t$5\t\.\t$6\t0\t$7;\n";
	}
}
close(FILE);
