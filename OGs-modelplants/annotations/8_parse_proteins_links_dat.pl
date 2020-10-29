#!/usr/bin/perl -w
use strict;
# Usage: Extraire des fonctions des proteines de MetaCyc

my $dirWork = $ARGV[0]; # Repertoire de travail (git clone)
my $version = $ARGV[1];  # Version of Metacyc (meta.tar.gz) dans le dossier annot_source
my $tag = '_'.$ARGV[2]; # tag du run (facultatif)
my $dirIn = $dirWork.'/annot'.$tag.'/annot_source/'.$version; # dossier metacyc 
my $dirOut = $dirWork.'/annot'.$tag.'/annot_files/'; # Repertoire contenant les fichiers utilisés lors de l'annotation des groupes
mkdir($dirOut) unless (-e $dirOut);
my $filelinks = $dirIn.'/data/protein-links.dat'; #inputfile

# Lecture du fichier decrivant toutes les proteines
if (-e $filelinks) 
{
	# pour chaque proteine du fichier, on lit son ID de metacyc (ATTRIBUTE=UNIQUE-ID) et on lit les identifiants de la proteine dans d'autres bases de donnees (ATTRIBUTE=DBLINKS)
	open IN,$filelinks or die ("pb avec l'ouverture du fichier $filelinks\n");	
	while (my $l = <IN>) 
	{
		chomp($l);		
		unless ($l=~ /^#/)
		{
			my @l = split (/\t/, $l);
			if(!$l[1]) { $l[1] = "NR1";}
			if(!$l[2]) { $l[2] = "NR2";}
			if(!$l[3]) { $l[3] = "NR3";}			
			print $l[0]." ; ".$l[1]." ; ".$l[2]." ; ".$l[3]."\n";
			
		}
		
	}
}

		
				
				
