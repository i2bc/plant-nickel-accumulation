#!/usr/bin/perl -w
use strict;
use threads;
# Execution en local de Inparanoid entre deux genomes. 
# Inparanoid n'execute pas les BLASTs mais utilise les BLASTs déjà calculés 
# mais qu'il utilise les resultats BLAST executes au prealable sur le cluster

my $sp1 = $ARGV[0];
my $sp2 = $ARGV[1];
my $blast_dir = $ARGV[2]; # INPUT : repertoire contenant les resultats BLAST realises au prealable
my $out_dir = $ARGV[3]; # OUTPUT 
my $max_threads = $ARGV[4]; # allowed multithreading

sleep(1) while(threads->list(threads::running) >= $max_threads);

sub launch_inparanoid 
{
  my ($sp1,$sp2) = @_;
  system("perl inparanoid_4.1/inparanoid_blastp.pl $sp1.fa $sp2.fa $blast_dir $out_dir");
}
