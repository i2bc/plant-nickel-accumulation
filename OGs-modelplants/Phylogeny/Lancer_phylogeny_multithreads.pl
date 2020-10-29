#!/usr/bin/perl 
#===============================================================================
#
#         FILE:  Lancer_phylogeny.pl
#
#        USAGE:  ./Lancer_phylogeny.pl  
#
#  DESCRIPTION:  lance l'ensemble des scripts de phylogeny
#
#      OPTIONS:  ---
# REQUIREMENTS:  ---
#         BUGS:  ---
#        NOTES:  ---
#       AUTHOR:  YOUR NAME (), 
#      COMPANY:  
#      VERSION:  1.0
#      CREATED:  10/06/2014 15:45:27
#     REVISION:  ---
#===============================================================================

#use strict;
use warnings;
use lib "/home/chris/fungipath/lib/";

$out=$ARGV[0]; # dossier resultats 
$cpu=$ARGV[1];

$work= "/home/chris/fungipath/Phylogeny/";
$fasta= "/home/chris/fungipath/DATA/Proteomes/";
$blast= "/home/chris/fungipath/blastp";
$evalue_cutoff = 0.001;
$alignment_cutoff = 0;

mkdir($out);

#Usage?
@nbfasta=glob("$fasta/*"); 
$nb=@nbfasta;
$nbmax=$nb*5;
#die;
system("date");
print "1_phylo\n";
print "system perl 1_phylogeny_launch_jobs.pl $work $fasta $blast $out $evalue_cutoff $alignment_cutoff $cpu";
system("perl 1_phylogeny_launch_jobs.pl $work $fasta $blast $out $evalue_cutoff $alignment_cutoff $cpu");
system("date");
#die;
print "2_phylo\n";
system("perl 2_phylogeny_launch_jobs.pl $work $out $fasta ");
system("date");

print "3_phylo\n";
system("perl 3_phylogeny_launch_jobs.pl $work $nbmax $out $cpu");
system("date");
#die;
print "4_phylo\n";
system("perl 4_phylogeny_launch_jobs.pl $work $fasta $out $cpu");
system("date");
print "5_phylo\n";
system("perl 5_phylogeny_modif.pl $out $fasta ");
system("date");
