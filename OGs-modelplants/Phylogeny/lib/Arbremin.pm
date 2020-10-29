 #!/usr/bin/perl -w

# Argument 1: fichier contenant l'arbre e analyser
# Ce programme detecte les orthologues e partir d'un arbre
# au format Newick donne en argument. Cet arbre doit etre enracine
# comme il faut pour que le programme fonctionne. Par exemple, point milieu.
# Ensemble de fonctions pour parcourir des arbres... -> Recursivite e
# mort!
# Les arbres sont implementes avec des tableaux de hashage -> Un pour
# les enfants, et un pour les parents.


# Reste e reprendre la distance PAM pour pouvoir decider d'une paire
# d'Orthologues lorsque on a une proteine face e groupe de paralogues
# et lorsqu'on a un groupe de paralogues face e un groupe de paralogues.

package Arbre;

use strict;
# La racine de l'arbre
use constant RACINE => 0;

# constructeur prenant en parametre un fichier de graphe et un fichier contenant la distance pam entre toutes les paires de proteines
sub newDPam
{
    my ($classe,$fichier,$fichier_dPam) = @_;
    my $this = {};

    unless(defined $fichier){die "You have to give the file name in Arbre->new(file)\n";}

    bless $this,$classe;

    #Nom de l'arbre
    $this->{nom} = $fichier;
    # Arbre avec les fleches orientees vers les parents
    # pour chaque cle (id du noeud) on a une valeur: le parent du noeud en cle
    $this->{parent} = {};
    # Arbre avec les fleches orientees vers les enfants
    # Pour chaque cle (id du noeud) on a une valeur: tableau d'enfants;
    $this->{enfants} = {};
    # Garde les proteines qui appartiennent a des groupes de paralogues
    $this->{para}={};
    # Hash contenant toutes les distances PAM entre les proteines du groupes
    # correspondant \`a l'arbre
    $this->{dpam}={};

    # Lecture du fichier et analyse de l'arbre
    $this->analyse($fichier);
    $this->analyse_dPam($fichier_dPam);

    # on retourne la reference vers l'objet lui meme
    return $this;
}


##################################################
#                                                #
#         Definition des Methodes                #
#                 de la classe                   #
#                                                #
##################################################

# Analyse la chaine de caractere passee en argument, et representant
# l'arbre. Construit alors l'arbre en memoire, avec le hash parent, et
# le hash enfants.
sub analyse
{
    my ($this,$fichier) = @_;

    # On essaie d'ouvrir le fichier d'arbre
    open(IN,$fichier) or die "The tree file $fichier doesn't exists";
    # lecture du fichier
    my $s = "";

    my $taille = 0;

    while(<IN>)
    {
		chomp;
		$s .= $_;
    }
    close(IN);

    $s =~ s/\;//g;
    $s =~ s/:.*?([,\)])/$1/g; #suppression des valeurs des distances dans l'arbre PhyML
    #Ex: (((8_490,((((2_3593,4_1109),5_4032),1_3802),7_1285)),9_2887),(6_3486,(3_4661,0_1028)))
    
    # On transforme la chaine en un tableau de caractere
    my @t = split(//,$s);

    # Initialisation des variables de la boucle
    my $nbNoeud = -1; # pas de noeud
    my $noeudCourant; # pas de noeud courant   
    my $species = ""; # pas d'espece courante non plus
    
    # On parcourt le tableau representant la chaine de caractere
    for (my $i=0;$i<=$#t;$i++)
    {
		# $v est le ieme +1 caractere de l'arbre PhyML simplifié
		# soit un des 3 caracteres specifiques du format Newick [(,)] soit un caractere appartenant a une nom de proteine	
		my $v = $t[$i];
		# Si on rencontre une parenthese ouvrante, on cree un noeud
		if($v eq "(")
		{
			if($nbNoeud>=0) # il ne s'agit pas de la premiere '(' de l'arbre
			{
				my $temp = $nbNoeud+1; #valeur du noeud actuel (= variable numerique incrementale)
				$this->{parent}{$temp} = $noeudCourant; # valeur du noeud precedent/parent 
				push @{$this->{enfants}{$noeudCourant}},$temp; # enregistre la valeur du noeud actuel dans la table des enfants du noeud precedent
				
				#print "test liste enfants apres ( \n";
				#print "Noeud". $noeudCourant."\n";
				#foreach $v (@{$this->{enfants}{$noeudCourant}}) { print $v."\n";}
				#print "//\n";
			}
			else
			{
				$this->{parent}{0} = "nil"; # premiere parenthese rencontrée, pas de parents
			}
			$nbNoeud++; 
			$noeudCourant=$nbNoeud;	# actualisation de la valeur du noeud precedent a la valeur actuelle		
		}
		# Si on rencontre une virgule, on peut ajouter l'espece (ou plutot la proteine????, CD)
		# courante comme feuille
		elsif($v eq ",") #
		{
			# 2 cas 
			# 1) la virgule relie 2 proteines enfants du même noeud 			
			if($t[$i-1] ne ")")
			{
				$this->{parent}{$species} = $noeudCourant;
				# on enregistre la proteine (contenu de la variable Species) comme enfant (feuille de l'arbre) du noeud precedent
				push @{$this->{enfants}{$noeudCourant}},$species; 
				
				#print "test liste enfants apres virgule\n";
				#print "Noeud". $noeudCourant."\n";
				#foreach $v (@{$this->{enfants}{$noeudCourant}}) { print $v."\n";}
				#print "//\n";
			}
			# cas 1 (= premiere proteine d'un couple) 
			# cas 2 = la virgule est precedee par une parenthese fermante
			# re-initialisation de la variable species/proteine apres ou non ecriture d'une proteine
			$species=""; 
			
		}
		# Si on rencontre une parenthese fermante, on remonte d'un
		# niveau dans l'arbre et on ajoute l'especes courante comme
		# enfant du noeud courant.
		elsif($v eq ")")
		{
			if($noeudCourant eq "nil"){die "bad formed tree file\n";}
			if($noeudCourant>=0)
			{
				if($t[$i-1] ne ")") # la parenthese precede un nom de proteine
				{
					$this->{parent}{$species} = $noeudCourant; # enregistrement du noeud precedent comme parent du noeud actuel
					push @{$this->{enfants}{$noeudCourant}},$species; # enregistrement de la proteine comme enfant du noeud precedent
					
					#print "test liste enfants apres ) \n";
					#print "Noeud". $noeudCourant."\n";
					#foreach $v (@{$this->{enfants}{$noeudCourant}}) { print $v."\n";}
					#print "//\n";
				}
				$noeudCourant = $this->{parent}{$noeudCourant}; # on remonte d'un niveau dans l'arbre
			}
			$species="";
		}
		# tous autres caracteres que (), on reconstruit le nom d'une proteine
		else 
		{ 
			$species.=$v; 
		} 
    }
}

# Analyse du fichier contenant les distances PAM (fichiers groupe_<numero du groupe> dans STEP2_<tag> posterieurement au run de mcl (STEP3)
# recupere la valeur de distance pour les 2 proteines de chaque edge dans un hash
sub analyse_dPam
{
    my ($this,$fichier) = @_;
    if (-e $fichier)
    {
        open(IN,$fichier); 
        while(<IN>)
        {
            if(/\(edge/)
            {
                s/\(|\)//g;
                my @temp = split(/\s/);
                $this->{dpam}{$temp[1]}{$temp[2]}=$temp[$#temp];
                $this->{dpam}{$temp[2]}{$temp[1]}=$temp[$#temp];
            }
        }
        close(IN);
    }
    else
    {
        print "No file $fichier\n";
    }
}

# Analyse l'arbre en question pour savoir si c'est un arbre:
#   - simple (que especes differentes):    retourne 1
#   - avec seulement des inparalogs:       retourne 2
#   - avec des out-paralogs ET/OU in-paralogs: retourne 3

sub simplicite
{
    my ($this) = @_;
    
    # On liste toutes les proteines de l'arbre tousEnfants(RACINE);
    if($this->especesDifferentes(RACINE))
    {
		# especesDifferentes return a protein list;
		return "1";
    }
    #elsif($this->outpara(RACINE))
    elsif ($this->para(RACINE)) # une seule espece pour toutes les proteines
    {
		
		#para retourne le genome commun aux proteines\n";
		return "2";
    }
    else
    {
		# print "other case simplicite return 2\n";
		return "3";
    }
}

# Comme la methode plus bas, mais prend comme 2 enfants tout groupe PARA
# ? (CD)
# le script est dependant de la structure des nom de proteines! (CD)
# initialement les 5 premieres lettres correspondent a un code de l'espece
# code modifié pour l'adapter a la numerotation préalable des genomes (<numerotation des genomes>_<numerotation des proteines du fasta du genome)
sub especesDifferentes
{
    my ($this,$noeud) = @_;

    if(!(defined $noeud) || (!exists( $this->{enfants}{$noeud}) && !(exists $this->{parent}{$noeud})))
    {
		die "The node $noeud doesn't exists in the tree ".$this->{nom}."";
    }

    # On prend tous les enfants du noeud (CD: la liste ne comprends que les proteines/feuilles de l'arbre)
    my $t = $this->tousEnfants($noeud);
    my %o; # hash qui enregistre la presence d'une espece dans l'arbre vaut 1 si on rencontre une nouvelle espece 
	
	# exception: le noeud est une feuille
    if($#{$t} == 0)
    {
		return 0;
    }
	# parcourt de la liste des proteines
    foreach my $p (@$t) 
    {
		if ($p =~ /_/)
		{
			#print "protein is $p\n";		
			#my $s = substr($p,0,5);
			my @s = split(/_/,$p); 
			my $s = $s[0]; # numero du genome de la proteine $p
			#print "genome is $s\n";
			if(exists $o{$s})
			{
				# l'espece a deja ete vue pour une precedente proteine, au moins 2 proteines sont de la même espece
				return 0;
			}
			else
			{
				$o{$s} = 1; # enregistrement de l'espece
			}
		}
    }   
    return $t; #si toutes les especes sont differentes on renvoie leur liste 
}

# Dit si l'arbre a des duplications pas seulement aux feuilles
# Pour cela: parcour l'arbre, et si il y a intersection et que
# le noeud n'est pas aux feuilles, alors --> renvoie vrai
sub outpara
{
    my ($this,$noeud) = @_;

    if(!(defined $noeud) || (!exists( $this->{enfants}{$noeud}) && !(exists $this->{parent}{$noeud})))
    {
		die "The node doesn't exists";
		
    }

    if($this->duplication($noeud) && $this->profondeur($noeud)>=3)
    {
		return 1;
    }
    elsif($this->feuille($noeud))
    {
		return 0;
    }
    else
    {
		foreach my $e (@{$this->{enfants}{$noeud}})
		{
			#print "noeud:".$e."\n";	    
			if($this->outpara($e) && $this->profondeur($e)>=3)
			{
				return 1;
			}
		}
    }
    return 0;
}

# Retourne la profondeur du sous arbre partant du noeud en argument
sub profondeur
{
    my ($this,$noeud) = @_;
    my @profondeur;

    if(!(defined $noeud) || (!exists( $this->{enfants}{$noeud}) && !(exists $this->{parent}{$noeud})))
    {
		die "The node doesn't exists";
    }

    if ($this->feuille($noeud))
    {
		return 1;
    }
    else
    {
		foreach my $e (@{$this->{enfants}{$noeud}})
		{
			push @profondeur,($this->profondeur($e)+1);
		}
		return max(\@profondeur);
    }
}

# Retourne le nombre de Noeuds de l'Arbre
sub taille
{
    my ($this) = @_;
    return ($#{$this->tousEnfants(RACINE)}+1);
}

# Retourne vrai (1) si le noeud en argument est une feuille
# et faux sinon
sub feuille
{
    my ($this,$noeud) = @_;
    if(!(defined $noeud) || (!exists( $this->{enfants}{$noeud}) && !(exists $this->{parent}{$noeud})))
    {
		die "The node doesn't exists";
    }
    return (!exists($this->{enfants}{$noeud}));
}



# retourne les enfants directes du noeud en question
sub enfantsD
{
    my ($this,$noeud) = @_;
    if(!(defined $noeud) || (!exists( $this->{enfants}{$noeud}) && !(exists $this->{parent}{$noeud})))
    {
		die "The node doesn't exists";
    }

   return $this->{enfants}{$noeud};
}

# Affiche l'arbre
sub affiche
{
    my ($this,$noeud,$prefix) = @_;
    my @temp = split(/\//,$noeud);
    my $n = $temp[$#temp];
    foreach my $f (@{$this->{enfants}{$n}})
    {
        if($f eq $this->{enfants}{$n}->[$#{$this->{enfants}{$n}}])
        {
			# pour la derniere entree du repertoire, on met un `
			# plutet qu'un |
			print $prefix.'`--'.$f."\n";
			my $fullname = $noeud."/".$f;
			if (exists $this->{enfants}{$f})
			{
				$this->affiche($fullname,$prefix.'   ');
			}
		}
		else
		{
            print $prefix.'|--'.$f."\n";
			my $fullname = $noeud."/".$f;
			if (exists $this->{enfants}{$f})
			{
				$this->affiche($fullname,$prefix.'|  ');
			}
		}
    }
}


# affiche les enfants de tous les noeuds
sub enfants
{
    my ($this) = @_;
    foreach my $v (keys(%{$this->{enfants}}))
    {
		print "$v";
		foreach my $e (@{$this->{enfants}{$v}})
		{
			print " - $e";
		}
		print "\n";
    }
}


# Affiche les parents de tous les noeuds
sub parents
{
    my ($this) =@_;
    foreach my $v (keys(%{$this->{parent}}))
    {
		print "$v - ".$this->{parent}{$v}."\n";
    }
}

# Renvoie un tableau contenant tous les enfants d'un noeuds, en
# considerant les groupes PARA comme un seul noeud/enfant.
sub enfantsPara
{
    my ($this,$noeud) = @_;
    if(!(defined $noeud) || (!exists( $this->{enfants}{$noeud}) && !(exists $this->{parent}{$noeud})))
    {
		die "The node doesn't exists";
    }
    my @tab;
    if(exists($this->{enfants}{$noeud}))
    {
		my $sp = $this->para($noeud);
		if($sp ne "0")
		{
			my @t = ($sp."_para_".$noeud);
			return \@t;
		}
		else
		{
			foreach my $e (@{$this->{enfants}{$noeud}})
			{
				push @tab,@{$this->enfantsPara($e)};
			}
			return \@tab;
		}
    }
    my @t = ($noeud);
    return(\@t);
}

#renvoie un tableau contenant tous les enfants d'un noeud (indirectes)
sub tousEnfants
{
    my ($this,$noeud) = @_;
    my @tab;
    if(!(defined $noeud) || (!exists( $this->{enfants}{$noeud}) && !(exists $this->{parent}{$noeud}))){
	die "The node doesn't exists";
    }

    if(exists($this->{enfants}{$noeud}))
    {
		foreach my $e (@{$this->{enfants}{$noeud}})
		{
			push @tab,@{$this->tousEnfants($e)};
		}
		return \@tab;
    }
    my @t = ($noeud);
    return(\@t);
}

# Regarde si tous les enfants appartiennent e la meme espece. Si non:
# renvoie 0, si oui: renvoie l'especes en question
sub para
{
    my ($this,$noeud) = @_;

    if(!(defined $noeud) || (!exists( $this->{enfants}{$noeud}) && !(exists $this->{parent}{$noeud}))){
	die "The node doesn't exists";
    }

    my $t = $this->tousEnfants($noeud);
    
    # CD test
    #foreach my $test (@$t)
    #{
		#print "liste in para tousenfants:".$test."\n";
    #}
    # fin
    
    #my $species = substr($t->[0],0,5);
    my @species = split(/_/, $t->[0]); # numero du genome de la premiere proteine
    #print "first species is $species[0]\n";
    foreach my $p (@$t)
    {
		my @other_sp = split (/_/, $p);
		if ($other_sp[0] ne $species[0])
		#if(substr($p,0,5) ne $species)
		{
			#print $species[0].','.$other_sp[0]."\n";
			return "0"; 
		}
    }
    # On ajoute au hash contenant les groupes de paralogue, le groupe
    # en question
    if(!exists $this->{para}{$species[0]."_para_".$noeud}){push @{$this->{para}{$species[0]."_para_".$noeud}},@$t};
    #print "toutes les proteines appartiennent au genome $species[0] \n"; 
    return "sp".$species[0];
}


# Renvoie vrai si le noeud en question provient d'une duplication
# et faux sinon (faux: speciation)
sub duplication
{
    my ($this,$noeud) = @_;
	
	# exceptions
    if(!(defined $noeud) || (!exists( $this->{enfants}{$noeud}) && !(exists $this->{parent}{$noeud})))
    {
		die "The node doesn't exists";
    }
    if($this->feuille($noeud))
    {
		return 0;
    }
	# extraction de la liste de toutes les proteines des noeuds
    my @tabEnfants = @{$this->{enfants}{$noeud}};
    print "In duplication/tabEnfants\n";
    foreach my $enfant(@tabEnfants) { print "$enfant\t";}
	
    
    # On considere qu'il n'y a pas d'intersection entre les enfants
    my $intersect = 0;
    for(my $i=0;$i<=$#tabEnfants && !$intersect;$i++)
    {
		for(my $j=$i+1;$j<=$#tabEnfants && !$intersect;$j++)
		{
			if(intersect($this->tousEnfants($tabEnfants[$i]),$this->tousEnfants($tabEnfants[$j])))
			{
				$intersect = 1 ;
			}
		}
    }
    print "intersect value is $intersect\n";
    print "//\n";
    return $intersect;
}

# Renvoie vrai s'il existe une intersection non nulle 
# entre les 2 multi-ensembles en argument ( tableaux)
# C'est une intersectin au niveau des noms d'especes
# Methode Statique: pour l'appeler, ne pas la ratacher e
# un obje de type arbre mais juste intersect(\@dskd,\@jlj);
sub intersect
{
    my ($e1,$e2) = @_;
    my $intersect = 0;

    # On parcour la premiere liste
    for(my $i=0;$i<=$#{$e1} && !$intersect;$i++)
    {
		# On parcour la deuxieme liste
		for(my $j=0;$j<=$#{$e2} && !$intersect;$j++)
		{
			#print "$e1->[$i] $e2->[$j]\n";
			my @e1 = split (/_/, $e1->[$i]); # $e1[0] est le numero de genome de la proteine $e1->[$i]
			my @e2 = split (/_/, $e2->[$j]); # $e2[0] est le numero de genome de la proteine $e2->[$i]
			#if(substr($e1->[$i],0,5) eq substr($e2->[$j],0,5)){
			if($e1[0] eq $e2[0]) # même espece pour 2 proteines de l'arbre
			{
				#print " yes";
				$intersect = 1;
			}
			#print "\n";
		}
    }
    return $intersect;
}

# Regarde si tous les enfants sont d'especes differentes
# Compte comme un seul enfants tout groupe PARA
# Renvoie null si le noeud ne constitue pas un groupe d'orthologues,
# et un tableau contenant le groupe d'orthologue sinon
sub ortho
{
    my ($this,$noeud) = @_;

    if(!(defined $noeud) || (!exists( $this->{enfants}{$noeud}) && !(exists $this->{parent}{$noeud}))){
	die "The node doesn't exists";
    }

    # On prend tous les enfants du noeud, enconsiderant les groupes
    # PARA comme un seul enfant
    my $t = $this->enfantsPara($noeud);
    my %o;

    if($#{$t} == 0){
	return 0;
    }

    foreach my $p (@$t)
    {
		my @s= split (/_/,$p);
		my $s = $s[0];
		#my $s = substr($p,0,5);
		if(exists $o{$s})
		{
			return 0;
		}
		else
		{
			$o{$s} = 1;
		}
    }
    return $t;
}
# Appel dans my4_phylogeny_job2.pl : $arbre->parcourOrtho(Arbre::RACINE,\@ortho);
sub parcourOrtho
{
    my ($this,$noeud,$ortho) = @_;
    print "\nparcourOrtho $noeud\n";

    if(!(defined $noeud) || (!exists($this->{enfants}{$noeud}) && !(exists $this->{parent}{$noeud})))
    {
	die "The node $noeud doesn't exists";
    }
    
    if($this->duplication($noeud)) # les enfants du noeud sont 2 proteines de la même espece
    {
		print "after duplication\n";
		foreach my $n (@{$this->enfantsD($noeud)})
		{
			print "$n\t";
			$this->parcourOrtho($n,$ortho);
		}
		print "//\n";
    }
    else
    {
		my $tab = $this->enfantsD($noeud);
		print "tab/enfantsD: \t";
		for ( my $k = 0 ; $k <= $#{$tab} ; $k++)
		{
			print $tab->[$k],"\t";
		}
		print "//\n";
				
		for (my $i = 0 ; $i <= $#{$tab} ; $i++)
		{
			for(my $j = $i+1 ; $j <= $#{$tab} ; $j++)
			{
				foreach my $e1 ( @{$this->tousEnfants($tab->[$i])})
				{
					foreach my $e2 (@{$this->tousEnfants($tab->[$j])})
					{
						if(defined $this->{dpam}{$e1}{$e2})
						{
							push(@$ortho,[$e1,$e2,$this->{dpam}{$e1}{$e2}]);
							#print STDERR"$e1\t$e2\t".$this->{dpam}{$e1}{$e2}."\n";
						}
						else
						{
							push(@$ortho,[$e1,$e2]);
							#print STDERR "$e1\t$e2\t\n";
						}
			#else{
			#    print "$e1 - $e2 - ?\n";
			#}
					}
				}
			}
		}
		foreach my $e (@$tab)
		{
			#print $e."\n";
			$this->parcourOrtho($e,$ortho);
		}	
    }
}

# Renvoie le maximum du tableau en argument
# Methode Statique
sub max
{
    my ($tab) = @_;
    my $temp = $tab->[0];
    foreach my $v (@$tab)
    {
		if($v > $temp)
		{
			$temp = $v;
		}
    }
    return $temp;
}
1;
