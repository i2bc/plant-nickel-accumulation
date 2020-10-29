#!/usr/bin/perl -w
use strict;

# chaque job va creer un fichier par paire de genomes (*tri) 
# chaque job permet de parser les fichiers blast generes pour extraire les resultats repondant aux criteres fixes
# les fichiers *tri seront utilises ensuite pour construire les familles de protéines
# script adapté pour pouvoir etre lancé sur le cluster du lix
# a besoin du fichier 1_phylogeny_job_clusterLIX.pl

#[cluster]perl 1_phylogeny_launch_jobs.pl /home/cpereira/ 0.001 0 

#############
# Declaration
#############

# Repertoires
#my $dirWork = '/projet/extern/sgrossetete/'; # INPUT : Repertoire de travail 
my $dirDONNEES=$ARGV[0];
my $dirWork ="$dirDONNEES/Phylogeny/";
my $dirIn = $dirDONNEES.'data/';  # INPUT : repertoire contenant donnees
my $dirGenomeNew = $dirIn; # INPUT : repertoire contenant les nouveax genomes au format fasta
my $dirBlast = $dirDONNEES.'blastp/'; # INPUT : Repertoire contenant les fichiers BLAST generes lors de la methode BRH 
my $dirOut = $dirWork.'result/dataPhylo_1e3_0/'; # OUTPUT : Repertoire ou sera enregistre les fichiers generes par le programme
my $dirTMP=$dirWork.'result/TMP_phylo/';
# Seuil a respecter sandrine avait donné 0.001 et 0
die ("Erreur : Trois parametres doivent etre donnes en argument, le dossier de travail, le seuil sur la E-value et celui sur le pourcentage d'alignement\n") unless ($#ARGV == 2);
my $evalue = $ARGV[1]; #seuil sur la evalue pour former familles d'homologues
my $palign = $ARGV[2]; #seuil sur le pourcentage d'alignement pour former familles d'homologues

# Autres variables
my @File = ();
my @AllGenome = ();
my $esp1 = ''; #nom de l'espece1
my $esp2 = ''; #nom de l'espece2
my $file = '';
my $diresp1 = ''; #nom du directory dans lequel on place les blasts de esp1
my $diresp2 = '';
my $nvnomesp1 = '';
my $nvnomesp2nv= '';
my $nvnomesp2ac='';
my $nomdecompresse = '';
my @File2=();
my @File3=();
my @File4=();
my $dossierTMP="";
my $dossierTMPanc="";


###########
# Programme
###########

# affichage seuil
print "Resultat enregistre dans :".$dirOut;
print "\nEvalue seuil : ".$evalue;
print "\npourcentage alignement seuil : ".$palign."\n";

# Creation du repertoire de sortie si non cree
print "creation des repertoires de sortie\n";
system("mkdir $dirWork/result/")unless (-e $dirWork.'/result');
system("mkdir $dirOut") unless (-e $dirOut);
system("mkdir $dirTMP") unless(-e $dirTMP);

# Recup du nom des genomes a analyser
print "recuperation des noms de genomes à analyser\n";
if (-e $dirGenomeNew){
	chdir($dirGenomeNew);     # Devient le repertoire de travail
	@File = glob("*.fa");             # Etabli liste du nom des genomes
	foreach $esp1 (@File){
		  $esp1 =~ s/\.fa$//;
		  push(@AllGenome,$esp1);
	}
}

#Execution des jobs entre un nouveau genome et tous les autres genomes
print "execution des jobs\n";
for (my $i = 0; $i <= $#AllGenome; $i ++) {
  $esp1 = $AllGenome[$i];
  unless((-e "$dirTMP/$esp1\_tri.err") && (-z "$dirTMP/$esp1\_tri.err")){#a moins que le fichier erreur de l'essait pressédent soit vide
	chdir($dirBlast);
	open(OUT,">".$dirTMP."/$esp1\_tri.qsub");
	unlink("$dirTMP/$esp1\_tri\.out");
	unlink("$dirTMP/$esp1\_tri\.err");
	print OUT "#!/bin/bash\n";
	print OUT '#$ -N '."$esp1\_tri\n";
	print OUT '#$ -cwd'."\n";
	print OUT '#$ -o '.$dirTMP."/$esp1\_tri\.out\n";
	print OUT '#$ -e '.$dirTMP."/$esp1\_tri\.err\n\n";
	# Comparaison de $esp1 avec tous les autres genomes
	for (my $j = $i; $j <= $#AllGenome; $j ++) {
		$esp2 =  $AllGenome[$j];
		#a moins que le fichier a deja été genere
		unless (((-e $dirOut.$esp2.'_'.$esp1.'.tri') && !(-z $dirOut.$esp2.'_'.$esp1.'.tri' )) || ((-e $dirOut.$esp1.'_'.$esp2.'.tri') && !(-z $dirOut.$esp1.'_'.$esp2.'.tri' ))) {
			#Lancement du job : extraction des resultats repondant aux criters entre ces deux genomes
			print OUT "perl $dirWork/1_phylogeny_job_clusterLIX.pl $dirOut $evalue $palign $esp1 $esp2 $dirBlast\n";
		}
	}
	close OUT;
	system("qsub $dirTMP/$esp1\_tri.qsub");
  	#suppression du script du job une fois le job lancé
  	unlink("$dirTMP/$esp1\_tri.qsub");
  }
}
