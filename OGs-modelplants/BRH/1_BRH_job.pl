#!/usr/bin/perl -w
use strict;
use File::Basename;


my $blastfile = $ARGV[0];
my $out_dir = $ARGV[1]; 

# Lecture du fichier genere par le BLAST
my $current_query_seq = 'none';
open (IN, $blastfile) or die("erreur blastfile");
open (OUT, '>'.$out_dir.'/'.basename($blastfile).'.besthit') or die("erreur outfile");

# Ecriture des Meilleurs hits
print $out_dir.'/'.basename($blastfile).'.besthit'."\n";
while (my $line = (<IN>)) 
{
    unless($line =~ /^#/)
	{
    	chomp($line);
    	my @Line = split(/\t/,$line);
		unless ($current_query_seq eq $Line[0]) 
		{
			#ecriture dans le fichier de sortie de l'ID du depart ainsi que de son meilleur hit, des longueurs alignees et du score
			print OUT "$Line[0]\t$Line[1]\t".($Line[7] - $Line[6] + 1)."\t".($Line[9] - $Line[8] + 1)."\t$Line[11]\t$Line[12]\t$Line[13]\n" ;
			$current_query_seq = $Line[0]; # ID de la sequence utilisee pour le BLAST
		}
	}
}
close(IN);
close(OUT);
