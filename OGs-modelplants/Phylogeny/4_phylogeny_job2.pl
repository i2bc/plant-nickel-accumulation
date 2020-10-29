#!/usr/bin/perl -w

#script recupere (??)

# Argument 1: dossier contenant les arbres! (chemin absolu!)
# Argument 2: dossier contenant la definition des groupes! (chemin absolu)
# Argument 3: le dossier contenant les fichiers tris (pour rechercher la dPam)
# Ecrit sur la sortie standard toutes les relations d'orthologie entre les proteines de chaques groupes/arbres

# Ce programme detecte les orthologues a partir d'un arbre au format newick donne en argument. 
# Cet arbre doit etre enracine comme il faut (??) pour que le programme fonctionne. 
# Par exemple, point milieu.
# Ensemble de fonctions pour parcourir des arbres... -> Recursivite e
# mort! (??)
# Les arbres sont implementes avec des tableaux de hashage:
# Un pour les enfants, et un pour les parents.


use strict;
use FindBin;
use lib "$FindBin::Bin";
use lib::MyArbre;

my $reparbr = $ARGV[0]; # dossier de sortie STEP4-1_<tag> = repertoire des arbres phyML apres retree (arbre_<numero du groupe>.txt)
my $repgrp = $ARGV[1]; # dossier STEP2_<tag> repertoires des fichiers groupe_<numero du groupe> = graphes decrit par des lignes nodes et edge
my $filearbr = $ARGV[2]; # fichier à traiter (jobs elementaires parallelisables)
my @ortho; 
my $fichier;


my $i = 0;

print STDERR "Analyse ... ".$reparbr.$filearbr."\n";
@ortho = ();

if(-f $reparbr.$filearbr)
{
	# extraction du numero du groupe
	my @temp = split(/\_/,$filearbr); 
	my $groupe = $temp[$#temp];
	$groupe=~s/\..*$//g;
	
	# On initialise un nouvel arbre 
	# utilise les methodes analyse (acquisition d'un arbre PhyML dans STEP4-1) et analyse_dpam (acquisition des donnee des fichiers groupe dans STEP2)
	my $arbre = newDPam MyArbre($reparbr.$filearbr,$repgrp."groupe_".$groupe);

	# On identifie les orthologues
	# parcourt de l'arbre à partir de la racine en utilisant deux hash (parent et enfants)
	#la fonction construit un tableau
	$arbre->parcourOrtho(MyArbre::RACINE,\@ortho);
	
	# On identifie les especes representees dans l'arbre PhyML par un code
	# 0 si seulement des orthologues (sp differentes), 1 si seulement des inparalogues, 2 si inpara et/ou outpara
	# Cette donnée figurera dans le fichier de sortie (1er champ de chaque ligne)
	my $s = $arbre->simplicite();
	#print "simplicite " .$s."\n";
	#exit();
	foreach my $o (@ortho)
	{
		# pour toutes les lignes la valeur sera la meme dans le fichier
		print $s."\t"; #0 si que esp diff, 1 si inpara seulement, 2 si inpara et/ou outpara
	    if(exists $o->[2]) # si il y a une distance (un ratio de scores) enregistré pour les 2 proteines
	    {
			print "$o->[0]\t$o->[1]\t$o->[2]\n"; #perso: id1 id2 ?
	    }
	    else
	    {
			# Si aucune distance n'est enregistre 0
			print "$o->[0]\t$o->[1]\t0\n";
	    }
	}
}
else
{
	print "No file ".$reparbr.$filearbr."\n";
}

