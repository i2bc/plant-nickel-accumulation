#!/usr/bin/perl 
use strict;
use warnings;
if ($#ARGV <4) { print "myLancer_BRH.pl  <dossier GIT > <repertoire des genomes fasta> <repertoire des blasts> <nombre de coeurs> <run tag>\n"; exit();}

my $work_dir = $ARGV[0];
my $genomes=$ARGV[1]; # directory where the genomes are stored (one file per genome suffixed .fa)
my $blastp=$ARGV[2]; # directory where the BLASTP results are stored
my $cpu=$ARGV[3]; # number of cpu to use
my $tag = $ARGV[4];

system("date");
print "STEP 1\n";
system("perl my1_BRH.pl $work_dir $blastp $cpu $tag");
system("date");
print "STEP 2\n";
system("perl my2_BRH.pl $work_dir $genomes $cpu $tag");
system("date");
print "STEP 3\n";
system("perl my3_BRH.pl $work_dir $genomes $tag");
system("date");
print "FIN";
