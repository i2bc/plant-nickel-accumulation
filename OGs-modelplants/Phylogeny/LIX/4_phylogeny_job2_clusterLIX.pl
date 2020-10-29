#!/usr/bin/perl -w

#script recupere

# Argument 1: dossier contenant les arbres! (chemin absolu!)
# Argument 2: dossier contenant la d�finition des groupes! (chemin absolu)
# Argument 3: le dossier contenant les fichiers tris (pour rechercher la dPam)
# Ecrit sur la sortie standard toutes les relations d'orthologie entre les prot�ines de chaques groupes/arbres

# Ce programme d�t�cte les orthologues � partir d'un arbre
# au format Newick donn� en argument. Cet arbre doit �tre enracin�
# comme il faut pour que le programme fonctionne. Par exemple, point milieu.
# Ensemble de fonctions pour parcourir des arbres... -> R�cursivit� �
# mort!
# Les arbres sont impl�ment�s avec des tableaux de hashage -> Un pour
# les enfants, et un pour les parents.


use strict;
# use lib '/mnt/data2/phylo';
use lib::Arbre;

my $reparbr = $ARGV[0];
my $repgrp = $ARGV[1];
my $filearbr = $ARGV[2];
my @ortho;
my $fichier;


my $i = 0;

print STDERR "Analyse ... ".$reparbr.$filearbr."\n";
@ortho = ();

if(-f $reparbr.$filearbr){
	my @temp = split(/\_/,$filearbr);
	my $groupe = $temp[$#temp];
	$groupe=~s/\..*$//g;
	# On initialise un nouvel arbre
	my $arbre = newDPam Arbre($reparbr.$filearbr,$repgrp."groupe_".$groupe);
	# On regarde ses orthologues
	$arbre->parcourOrtho(Arbre::RACINE,\@ortho);

	my $s = $arbre->simplicite();

	foreach my $o (@ortho){
	    print $s."\t"; #0 si que esp diff, 1 si inpara seulement, 2 si inpara et/ou outpara
	    if(exists $o->[2]){
			print "$o->[0]\t$o->[1]\t$o->[2]\n"; #perso: id1 id2 ?
	    }
	    else{
			# Si la distance pam n'existe pas dans l'arbre, on la recherche dans les fichiers
			print "$o->[0]\t$o->[1]\t0\n";
	    }
	}
}
else{
	print "No file ".$reparbr.$filearbr."\n";
}

