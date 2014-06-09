#!/usr/bin/perl
use strict;
open(IN,"/home/geparada/Desktop/refseq_tags_human.tabular") or die "Cannot open file gi"; 
open(OUT, "> out") or die "Cannot create file!";

while(<IN>) {
        chomp $_;
        my @v = split(/\t/,$_);
        get_pid(@v);

}
sub get_pid {
         my @line = @_;
         my $pid = (100.0 - (&pslCalcMilliBad(@line) * 0.1));
         print "@line $pid\n";
         return $pid;
}

sub pslCalcMilliBad {
         my @cols = @_;

         # sizeMul depens of dna/Prot
         my $sizeMul = 1;                  #核酸等于1，蛋白质等于3
        
     # cols[0]  matches
         # cols[1]  misMatches
         # cols[2]  repMaches
         # cols[4]  qNumInsert
         # cols[6]  tNumInsert
         # cols[11] qStart
         # cols[12] qEnd
         # cols[15] tStart
         # cols[16] tEnd

         my $qAliSize = $sizeMul * ($cols[12] - $cols[11]);
         my $tAliSize = $cols[16] - $cols[15];

         # I want the minimum of qAliSize and tAliSize
         my $aliSize;
         $qAliSize < $tAliSize ? $aliSize = $qAliSize : $aliSize =  
$tAliSize;

         # return 0 is AliSize == 0
         return 0 if ($aliSize <= 0);

         # size diff
         my $sizeDiff = $qAliSize - $tAliSize;
         if ($sizeDiff < 0) {
             $sizeDiff = 0;  
         }

         # insert Factor
         my $insertFactor = $cols[4]; 

# $insertFactor += $cols[6]; # 若为蛋白质

         my $milliBad = (1000 * ($cols[1]*$sizeMul + $insertFactor +  
&round(3*log( 1 + $sizeDiff)))) / ($sizeMul * ($cols[0] + $cols[2]
+ $cols[1]));
                return $milliBad;

}

sub round {
         my $number = shift;
         return int($number + .5);
}
