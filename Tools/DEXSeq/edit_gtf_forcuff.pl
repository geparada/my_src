############### CHANGE THE NAMES OF CHROMOSOMES, from 2R to chr2R, eliminates the Het annotations ########
open(FILE, $ARGV[0]) or die "I could not open the gtf file";
while(<FILE>){
   $_ =~ /^(\S+).*/;
   $chrname = $1;
   if($chrname !~ /(2L$|2R$|3L$|3R$|4$|X$|YHet$)/){
      next;
   }
   print "chr$_";
}
