#!/usr/bin/perl -w
#script de test 
# le script teste si tous les monomers des complexes figurent dans le fichier general proteins.dat
use strict;

my $dirWork = $ARGV[0]; # Repertoire de travail (git clone)
my $version = $ARGV[1];  # Version of Metacyc (meta.tar.gz) dans le dossier annot_source
my $tag = '_'.$ARGV[2]; # tag du run (facultatif)
my $dirIn = $dirWork.'/annot'.$tag.'/annot_source/'.$version; # dossier metacyc 
my $filelinks = $dirIn.'/data/protcplxs.col'; #inputfile

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
			#print $l[0]." ; ".$l[$#l]."\n";
			my @monomers = split(/\*|,/, $l[$#l]);		
			for (my $i =1; $i<$#monomers; $i= $i+2) 
			{ 
				my $test = `grep "UNIQUE-ID - $monomers[$i]" $dirWork/annot$tag/annot_source/20.1/data/proteins.dat `; 
				print $test;
				if (!$test ) 
				{ 
					#print "no $monomers[$i] in proteins.dat";
				}
				#else { print "$test=ok";}
			}
		}
		
	}
}

		
				
				
