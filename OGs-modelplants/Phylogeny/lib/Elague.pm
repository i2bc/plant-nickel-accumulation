#!/usr/bin/perl -w
package Elague;
use Graph;
use Graph::Undirected;
use strict;

sub elague{
    # Fichier d'entrée (contenant le graphe de départ)
    my ($fichier) = @_;
    # Nouveau graphe
    my $g = Graph::Undirected->new;
    # Sauvegarder le label des arcs -> distance PAM
    my %label;

    print "Mise en mémoire du Graphe\n";
    
    # On parcourt le fichier d'entrée
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

    print " ... Mise en mémoire Finie\n";

    print "Détection des arêtes \"Bridge\"\n";

    # On récupère les arêtes "bridges": Celles qui font du graphe un
    # ensemble non connexe
    # => Qui ne font pas partie d'un cycle
    my @b  = $g->bridges();
    
    print " ... Détection Finie\n";

    print "Suppression des bridges du graphe\n";

    # On parcourt les arcs "bridges"
    foreach my $bridge (@{b}){
	# on les enlèves du graphe
	$g->delete_edge($bridge->[0],$bridge->[1]);    
    }
    
    print " ... Suppression finie\n";

    print "Récupération des composantes connexes\n";

    # On récupère alors les composantes connexes du graphe ainsi obtenu
    my @cc = $g->connected_components();

    print " ... Récupération effectuée\n";
    
    # Nombre de sommets du graphe complet
    my $nbVertices = $g->vertices();
    print "Nombre de sommets total: $nbVertices\n";
    print "Après élagation:\n";
    my $noGraphe = 0;
    
    print "Parcour des composantes connexes pour les sauvegarder\n";

    # On parcourt les composantes connexes
    foreach my $c (@cc){
	my @edges;
	my %EDGE;
	# Si la composante n'est pas composée d'un seul sommet, on la
	# garde
	print $#{$c}."\n";
	if($#{$c}>0){
	    # On cré un fichier contenant cette composante connexe
	    open(OUT,">".$fichier."-".$noGraphe);
	    print OUT "(nodes ";
	    # On imprime les sommets
	    foreach my $v (@$c){
		my @e = $g->edges_at($v);
		push @edges,@e;
		print OUT "$v ";
	    }
	    print OUT ")\n";
	    # On imprime les arêtes et leur label (Distance PAM)
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
