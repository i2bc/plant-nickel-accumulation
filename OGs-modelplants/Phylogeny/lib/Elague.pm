#!/usr/bin/perl -w
package Elague;
use Graph;
use Graph::Undirected;
use strict;

sub elague{
    # Fichier d'entr�e (contenant le graphe de d�part)
    my ($fichier) = @_;
    # Nouveau graphe
    my $g = Graph::Undirected->new;
    # Sauvegarder le label des arcs -> distance PAM
    my %label;

    print "Mise en m�moire du Graphe\n";
    
    # On parcourt le fichier d'entr�e
    open(IN,$fichier);
    my $compt = 0;
    while(<IN>){
	chomp;
	s/\(|\)//g;
	# Pour chaque arc
	if(/edge/){
	    my @temp = split(/\s/,$_);
	    # On l'ajoute au graphe construit
	    $g->add_edges($temp[1],$temp[2]);
	    $label{$temp[1].".".$temp[2]} = $temp[3];
	}
	$compt++;
	if($compt%100000 eq 0){
		print $compt."\n";
	}
    }
    close(IN);

    print " ... Mise en m�moire Finie\n";

    print "D�tection des ar�tes \"Bridge\"\n";

    # On r�cup�re les ar�tes "bridges": Celles qui font du graphe un
    # ensemble non connexe
    # => Qui ne font pas partie d'un cycle
    my @b  = $g->bridges();
    
    print " ... D�tection Finie\n";

    print "Suppression des bridges du graphe\n";

    # On parcourt les arcs "bridges"
    foreach my $bridge (@{b}){
	# on les enl�ves du graphe
	$g->delete_edge($bridge->[0],$bridge->[1]);    
    }
    
    print " ... Suppression finie\n";

    print "R�cup�ration des composantes connexes\n";

    # On r�cup�re alors les composantes connexes du graphe ainsi obtenu
    my @cc = $g->connected_components();

    print " ... R�cup�ration effectu�e\n";
    
    # Nombre de sommets du graphe complet
    my $nbVertices = $g->vertices();
    print "Nombre de sommets total: $nbVertices\n";
    print "Apr�s �lagation:\n";
    my $noGraphe = 0;
    
    print "Parcour des composantes connexes pour les sauvegarder\n";

    # On parcourt les composantes connexes
    foreach my $c (@cc){
	my @edges;
	my %EDGE;
	# Si la composante n'est pas compos�e d'un seul sommet, on la
	# garde
	print $#{$c}."\n";
	if($#{$c}>0){
	    # On cr� un fichier contenant cette composante connexe
	    open(OUT,">".$fichier."-".$noGraphe);
	    print OUT "(nodes ";
	    # On imprime les sommets
	    foreach my $v (@$c){
		my @e = $g->edges_at($v);
		push @edges,@e;
		print OUT "$v ";
	    }
	    print OUT ")\n";
	    # On imprime les ar�tes et leur label (Distance PAM)
	    foreach my $e (@edges){
		if(!defined $EDGE{$e->[0]}{$e->[1]} && !defined $EDGE{$e->[1]}{$e->[0]}){
		    if(!exists $label{$e->[0].".".$e->[1]}){
			print OUT "(edge $e->[0] $e->[1] ".$label{$e->[1].".".$e->[0]}.")\n";	
			$EDGE{$e->[1]}{$e->[0]}=1;
		    }
		    else{
			print OUT "(edge $e->[0] $e->[1] ".$label{$e->[0].".".$e->[1]}.")\n";
			$EDGE{$e->[0]}{$e->[1]}=1;
		    }
		}
	    }
	    close(OUT);
	    $noGraphe++;
	}
    }
    unlink($fichier);

    print " ... Parcour et sauvegarde finis\n";

}
1;
