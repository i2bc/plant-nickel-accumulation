 #!/usr/bin/perl -w

# Argument 1: fichier contenant l'arbre � analyser
# Ce programme d�t�cte les orthologues � partir d'un arbre
# au format Newick donn� en argument. Cet arbre doit �tre enracin�
# comme il faut pour que le programme fonctionne. Par exemple, point milieu.
# Ensemble de fonctions pour parcourir des arbres... -> R�cursivit� �
# mort!
# Les arbres sont impl�ment�s avec des tableaux de hashage -> Un pour
# les enfants, et un pour les parents.


# Reste � reprendre la distance PAM pour pouvoir d�cider d'une paire
# d'Orthologues lorsque on a une prot�ine face � groupe de paralogues
# et lorsqu'on a un groupe de paralogues face � un groupe de paralogues.

package Arbre;

use strict;
# La racine de l'arbre
use constant RACINE => 0;


# Le constructeur prend en param�tre: le fichier d'entree
sub new
{
    my ($classe,$fichier) = @_;
    my $this = {};

    unless(defined $fichier){die "You have to give the file name in Arbre->new(file)\n";}

    bless $this,$classe;

    #Nom de l'arbre
    $this->{nom} = $fichier;
    # Arbre avec les fl�ches orient�es vers les parents
    # pour chaque cl� (id du noeud) on a une valeur: le parent du noeud en cl�
    $this->{parent} = {};
    # Arbre avec les fl�ches orient�es vers les enfants
    # Pour chaque cl� (id du noeud) on a une valeur: tableau d'enfants;
    $this->{enfants} = {};
    # Garde les prot�ines qui appartiennent � des groupes de paralogues
    $this->{para}={}; 
    # Lecture du fichier et analyse de l'arbre
    $this->analyse($fichier);

    # on retourne la r�f�rence vers l'objet lui m�me
    return $this;
}

# Deuxieme constructeur prenant en parametre en plus le fichier de graphe
# contenant la distance pam entre toutes les paires de proteines
sub newDPam
{
    my ($classe,$fichier,$fichier_dPam) = @_;
    my $this = {};

    unless(defined $fichier){die "You have to give the file name in Arbre->new(file)\n";}

    bless $this,$classe;

    #Nom de l'arbre
    $this->{nom} = $fichier;
    # Arbre avec les fl�ches orient�es vers les parents
    # pour chaque cl� (id du noeud) on a une valeur: le parent du noeud en cl�
    $this->{parent} = {};
    # Arbre avec les fl�ches orient�es vers les enfants
    # Pour chaque cl� (id du noeud) on a une valeur: tableau d'enfants;
    $this->{enfants} = {};
    # Garde les prot�ines qui appartiennent � des groupes de paralogues
    $this->{para}={};
    # Hash contenant toutes les distances PAM entre les proteines du groupes
    # correspondant \`a l'arbre
    $this->{dpam}={};

    # Lecture du fichier et analyse de l'arbre
    $this->analyse($fichier);

    $this->analyse_dPam($fichier_dPam);

    # on retourne la r�f�rence vers l'objet lui m�me
    return $this;
}


##################################################
#                                                #
#         D�finition des M�thodes                #
#                 de la classe                   #
#                                                #
##################################################

# Analyse la cha�ne de caract�re pass�e en argument, et repr�sentant
# l'arbre. Construit alors l'arbre en m�moire, avec le hash parent, et
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
	$s.=$_;
    }
    close(IN);

    $s=~s/\;//g;
    $s=~s/:.*?([,\)])/$1/g;
    # On met la cha�ne dans un tableau
    my @t = split(//,$s);

    # Il n'y a aucun noeuds au d�part
    my $nbNoeud = -1;
    # Pas de noeud courant
    my $noeudCourant;
    # Pas d'esp�ce courante non plus
    my $species = "";
    
    # On parcourt tous les tableau repr�sentant la chapine de caract�re
    for (my $i=0;$i<=$#t;$i++){
	# $v est la valeur du tableau � l'indice $i
	my $v = $t[$i];
	# Si on rencontre une parenth�se ouvrante, on cr� un nouveau
	# noeud � l'arbrbe
	if($v eq "("){
	    if($nbNoeud>=0){
		my $temp = $nbNoeud+1;
		$this->{parent}{$temp} = $noeudCourant;
		push @{$this->{enfants}{$noeudCourant}},$temp;
	    }
	    else{
		$this->{parent}{0} = "nil";
	    }
	    $nbNoeud++;
	    $noeudCourant=$nbNoeud;
	}
	# Si on rencontre une virgule, on peut ajouter l'esp�ce
	# courante comme feuille
	elsif($v eq ","){
	    if($t[$i-1] ne ")"){
		$this->{parent}{$species} = $noeudCourant;
		push @{$this->{enfants}{$noeudCourant}},$species;
	    }
	    $species="";
	}
	# Si on rencontre une parenth�se fermante, on remonte d'un
	# niveau dans l'arbre et on ajoute l'esp�ces courante comme
	# enfant du noeud courant.
	elsif($v eq ")"){
	    if($noeudCourant eq "nil"){die "bad formed tree file\n";}
	    if($noeudCourant>=0){
		if($t[$i-1] ne ")"){
		    $this->{parent}{$species} = $noeudCourant;
		    push @{$this->{enfants}{$noeudCourant}},$species;
		}
		$noeudCourant = $this->{parent}{$noeudCourant};
	    }
		$species="";
	}
	else{
	    $species.=$v;
	}
    }
}

# Analyse du fichier contenant les distances PAM
sub analyse_dPam{
    my ($this,$fichier) = @_;

    if (-e $fichier)
    {
        open(IN,$fichier); 
        while(<IN>){
            if(/\(edge/){
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
#   - simple (que esp�ces diff�rentes):    retourne 0
#   - avec seulement des inparalogs:       retourne 1
#   - avec des out-paralogs ET/OU in-paralogs: retourne 2
sub simplicite{
    my ($this) = @_;
    
    if($this->especesDifferentes(RACINE)){
	#print "0\n";
	return "0";
    }
    elsif($this->outpara(RACINE)){
	#print "2\n";
	return "2";
    }
    else{
	#print "1\n";
	return "1";
    }
}

# Dit si l'arbre a des duplicaions pas seulement aux feuilles
# Pour cela: parcour l'arbre, et si il y a intersection et que
# le noeud n'est pas aux feuilles, alors --> renvoie vrai
sub outpara{
    my ($this,$noeud) = @_;

    if(!(defined $noeud) || (!exists( $this->{enfants}{$noeud}) && !(exists $this->{parent}{$noeud}))){
	die "The node doesn't exists";
    }

    if($this->duplication($noeud) && $this->profondeur($noeud)>=3){
	return 1;
    }
    elsif($this->feuille($noeud)){
	return 0;
    }
    else{
	foreach my $e (@{$this->{enfants}{$noeud}}){
	    #print "noeud:".$e."\n";
	    
	    if($this->outpara($e) && $this->profondeur($e)>=3){
		return 1;
	    }
	}
    }
    return 0;
}

# Retourne la profondeur du sous arbre partant du noeud en argument
sub profondeur{
    my ($this,$noeud) = @_;
    my @profondeur;

    if(!(defined $noeud) || (!exists( $this->{enfants}{$noeud}) && !(exists $this->{parent}{$noeud}))){
	die "The node doesn't exists";
    }

    if ($this->feuille($noeud)){
	return 1;
    }
    else{
	foreach my $e (@{$this->{enfants}{$noeud}}){
	    push @profondeur,($this->profondeur($e)+1);
	}
	return max(\@profondeur);
    }
}

# Retourne le nombre de Noeuds de l'Arbre
sub taille{
    my ($this) = @_;

    return ($#{$this->tousEnfants(RACINE)}+1);

}

# Retourne vrai (1) si le neud en argument est une feuille
# et Faux sinon
sub feuille{
    my ($this,$noeud) = @_;

    if(!(defined $noeud) || (!exists( $this->{enfants}{$noeud}) && !(exists $this->{parent}{$noeud}))){
	die "The node doesn't exists";
    }

    return (!exists($this->{enfants}{$noeud}));
}

# retourne le parent du noeud en question
sub parentD{
    my ($this,$noeud) = @_;
    

    if(!(defined $noeud) || (!exists( $this->{enfants}{$noeud}) && !(exists $this->{parent}{$noeud}))){
	die "The node doesn't exists";
    }

    return $this->{parent}{$noeud};
}

# retourne les enfants directes du noeud en question
sub enfantsD{
    my ($this,$noeud) = @_;

    if(!(defined $noeud) || (!exists( $this->{enfants}{$noeud}) && !(exists $this->{parent}{$noeud}))){
	die "The node doesn't exists";
    }

   return $this->{enfants}{$noeud};
}

# Enracine l'arbre entre les deux noeuds en arguments
sub enracine{
    my ($this,$noeud1,$noeud2);

    print "Non implemente\n";
}

# Affiche l'arbre
sub affiche{
    my ($this,$noeud,$prefix) = @_;

    my @temp = split(/\//,$noeud);
    my $n = $temp[$#temp];

    foreach my $f (@{$this->{enfants}{$n}}){
        if($f eq $this->{enfants}{$n}->[$#{$this->{enfants}{$n}}]){
	    # pour la derni�re entr�e du r�pertoire, on met un `
	    # plut�t qu'un |
	    print $prefix.'`--'.$f."\n";
	    my $fullname = $noeud."/".$f;
	    if (exists $this->{enfants}{$f}){
		$this->affiche($fullname,$prefix.'   ');
	    }
	}
	else{
            print $prefix.'|--'.$f."\n";
	    my $fullname = $noeud."/".$f;
	    if (exists $this->{enfants}{$f}){
		$this->affiche($fullname,$prefix.'|  ');
	    }
	}
    }
}


# affiche les enfants de tous les noeuds
sub enfants{
    my ($this) = @_;
    foreach my $v (keys(%{$this->{enfants}})){
	print "$v";
	foreach my $e (@{$this->{enfants}{$v}}){
	    print " - $e";
	}
	print "\n";
    }
}


# Affiche les parents de tous les noeuds
sub parents
{
    my ($this) =@_;
    foreach my $v (keys(%{$this->{parent}})){
	print "$v - ".$this->{parent}{$v}."\n";
    }
}

# Renvoie un tableau contenant tous les enfants d'un noeuds, en
# consid�rant les groupes PARA comme un seul noeud/enfant.
sub enfantsPara
{
    my ($this,$noeud) = @_;


    if(!(defined $noeud) || (!exists( $this->{enfants}{$noeud}) && !(exists $this->{parent}{$noeud}))){
	die "The node doesn't exists";
    }

    my @tab;
    if(exists($this->{enfants}{$noeud})){
	my $sp = $this->para($noeud);
	if($sp ne "0"){
	    my @t = ($sp."_para_".$noeud);
	    return \@t;
	}
	else{
	    foreach my $e (@{$this->{enfants}{$noeud}}){
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

    if(exists($this->{enfants}{$noeud})){
	foreach my $e (@{$this->{enfants}{$noeud}}){
	    push @tab,@{$this->tousEnfants($e)};
	}
	return \@tab;
    }
    my @t = ($noeud);
    return(\@t);
}

# Regarde si tous les enfants appartiennent � la m�me esp�ce. Si non:
# renvoie 0, si oui: renvoie l'esp�ces en question
sub para
{
    my ($this,$noeud) = @_;

    if(!(defined $noeud) || (!exists( $this->{enfants}{$noeud}) && !(exists $this->{parent}{$noeud}))){
	die "The node doesn't exists";
    }

    my $t = $this->tousEnfants($noeud);

    my $species = substr($t->[0],0,5);
    foreach my $p (@$t)
    {
	if(substr($p,0,5) ne $species)
	{
	    return "0";
	}
    }
    # On ajoute au hash contenant les groupes de paralogue, le groupe
    # en question
    if(!exists $this->{para}{$species."_para_".$noeud}){push @{$this->{para}{$species."_para_".$noeud}},@$t};
    return $species;
}


# Renvoie vrai si le noeud en question provient d'une duplication
# et faux sinon (faux: sp�ciation)
sub duplication{
    my ($this,$noeud) = @_;

    if(!(defined $noeud) || (!exists( $this->{enfants}{$noeud}) && !(exists $this->{parent}{$noeud}))){
	die "The node doesn't exists";
    }

    if($this->feuille($noeud)){
	return 0;
    }

    my @tabEnfants = @{$this->{enfants}{$noeud}};

    # On consid�re qu'il n'y a pas d'intersection entre les enfants
    my $intersect = 0;
    for(my $i=0;$i<=$#tabEnfants && !$intersect;$i++){
 	for(my $j=$i+1;$j<=$#tabEnfants && !$intersect;$j++){
	    if(intersect($this->tousEnfants($tabEnfants[$i]),$this->tousEnfants($tabEnfants[$j]))){
		$intersect = 1 ;
	    }
	}
    }
    return $intersect;
}

# Renvoie vrai s'il existe une intersection non nulle 
# entre les 2 multi-ensembles en argument ( tableaux)
# C'est une intersectin au niveau des noms d'esp�ces
# M�thode Statique: pour l'appeler, ne pas la ratacher �
# un obje de type arbre mais juste intersect(\@dskd,\@jlj);
sub intersect{
    my ($e1,$e2) = @_;

    my $intersect = 0;

    # On parcour la premi�re liste
    for(my $i=0;$i<=$#{$e1} && !$intersect;$i++){
	# On parcour la deuxi�me liste
	for(my $j=0;$j<=$#{$e2} && !$intersect;$j++){
	    #print "$e1->[$i] $e2->[$j]\n";
	    if(substr($e1->[$i],0,5) eq substr($e2->[$j],0,5)){
		#print " yes";
		$intersect = 1;
	    }
	    #print "\n";
	}
    }
    return $intersect;
}

# Comme la m�thode plus bas, mais prend comme 2 enfants tout groupe PARA
sub especesDifferentes{
    my ($this,$noeud) = @_;

    if(!(defined $noeud) || (!exists( $this->{enfants}{$noeud}) && !(exists $this->{parent}{$noeud}))){
	die "The node $noeud doesn't exists in the tree ".$this->{nom}."";
    }

    # On prend tous les enfants du noeud
    my $t = $this->tousEnfants($noeud);
    my %o;

    if($#{$t} == 0){
	return 0;
    }

    foreach my $p (@$t)
    {
	my $s = substr($p,0,5);
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

# Regarde si tous les enfants sont d'esp�ces diff�rentes
# Compte comme un seul enfants tout groupe PARA
# Renvoie null si le noeud ne constitue pas un groupe d'orthologues,
# et un tableau contenant le groupe d'orthologue sinon
sub ortho
{
    my ($this,$noeud) = @_;

    if(!(defined $noeud) || (!exists( $this->{enfants}{$noeud}) && !(exists $this->{parent}{$noeud}))){
	die "The node doesn't exists";
    }

    # On prend tous les enfants du noeud, enconsid�rant les groupes
    # PARA comme un seul enfant
    my $t = $this->enfantsPara($noeud);
    my %o;

    if($#{$t} == 0){
	return 0;
    }

    foreach my $p (@$t)
    {
	my $s = substr($p,0,5);
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


sub parcourOrtho{
    my ($this,$noeud,$ortho) = @_;

    if(!(defined $noeud) || (!exists($this->{enfants}{$noeud}) && !(exists $this->{parent}{$noeud}))){
	die "The node $noeud doesn't exists";
    }
    
    if($this->duplication($noeud)){
	foreach my $n (@{$this->enfantsD($noeud)}){
	    $this->parcourOrtho($n,$ortho);
	}
    }
    else{
	my $tab = $this->enfantsD($noeud);
	for(my $i = 0 ; $i <= $#{$tab} ; $i++){
	    for(my $j = $i+1 ; $j <= $#{$tab} ; $j++){
		foreach my $e1 ( @{$this->tousEnfants($tab->[$i])}){
		    foreach my $e2 (@{$this->tousEnfants($tab->[$j])}){
			if(defined $this->{dpam}{$e1}{$e2}){
			    push(@$ortho,[$e1,$e2,$this->{dpam}{$e1}{$e2}]);
			    #print STDERR"$e1\t$e2\t".$this->{dpam}{$e1}{$e2}."\n";
			}
			else{
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
	foreach my $e (@$tab){
	    #print $e."\n";
	    $this->parcourOrtho($e,$ortho);
	}
	
    }
}

# Parcour tout l'arbre et regarde les noeuds qui amm�nent � des
# prot�ines orthologues. Pour cela, utilise la fonction "ortho"
sub parcourOrtho_OLD
{
    # On part du noeud en argument;
    my ($this,$noeud) = @_;
    
    if(!(defined $noeud) || (!exists( $this->{enfants}{$noeud}) && !(exists $this->{parent}{$noeud}))){
	die "The node doesn't exists";
    }

    if((my $t = $this->ortho($noeud)) != 0){
	print "\nNoeud $noeud:\n------------\n";
	# On va parcourir toutes les paires de prot�ines/para, pour
	# afficher les paires d'Orthologues
	for(my $i = 0; $i<=$#{$t};$i++){
	   #for(my $j=$i+1;$j<$#{$t};$j++){
	   #	# Plusieurs possibilit�s:
	   #	# - Soit $t->[$i] est une prot�ine et $t->[$j] aussi
	   #	# - Soit $t->[$i] est une prot�ine et $t->[$j] un groupe
	   #	#   de paralogue: Dans ce cas on prend la prot�ine la
	   #	#   plus proche de $t->[$i] dans $t->[$j]
	   #	# - Soit le cas inverse (m�me traitement)
	   #	# - Soit les deux sont des groupes de paralogues: Dans
	   #	#   ce cas, on regarde pour chaque prot�ine quelle est
	   #	#   sa plus proche.
	   #
	   #	# Cas o� aucun n'est un groupe de paralogue
	   #	if(!exists($para{$t->[$i]}) && !exists($para{$t->[$j]})){
	   #	    print $t->[$i]."-".$t->[$j]."\n";
	   #	}
	   #	
	   #	# Cas o� un des deux est un groupe de paralogue
	   #	elsif(!exists($para{$t->[$i]}) && exists($para{$t->[$j]}))
	   #	{
	   #	    # print $t->[$i]."-".$t->[$j].": ";
	   #	    foreach my $v (@{$para{$t->[$j]}}){
	   #		print $t->[$i]."-".$v."\n";
	   #	    }
	   #	}
	   #
	   #	# Cas o� l'autre est un groupe de paralogue
	   #	elsif(!exists($para{$t->[$j]}) && exists($para{$t->[$i]}))
	   #	{
	   #	    # print $t->[$i]."-".$t->[$j].": ";
	   #	    foreach my $v (@{$para{$t->[$i]}}){
	   #		print $v."-".$t->[$j]."\n";
	   #	    }
	   #	}
	   #	# Cas o� les deux sont des groupes de paralogues
	   #	else{
	   #	    # print $t->[$i]."-".$t->[$j].": ";
	   #	    foreach my $v (@{$para{$t->[$i]}}){
	   #		# print $t->[$i]."-".$t->[$j].": ";
	   #		foreach my $w (@{$para{$t->[$j]}}){
	   #		    print $v."-".$w."\n";
	   #		}
	   #	    }
	   #	}
	   #}
	    print $t->[$i]."\n";
	}
    }
    elsif(exists $this->{enfants}{$noeud}){
	foreach my $e( @{$this->{enfants}{$noeud}}){
	    $this->parcourOrtho($e);
	}
    }
}


# Renvoie pour l'instant tous les noeuds dont tous les enfants 
# sont paralogues c'est � dire dont tous les enfants sont de la m�me esp�ce
sub fusionnePara
{
    # Noeud initial
    my ($this,$noeud) = @_;
    my @tabTemp;
    
    if(!(defined $noeud) || (!exists( $this->{enfants}{$noeud}) && !(exists $this->{parent}{$noeud}))){
	    die "The node doesn't exists";
    }
    
    if($this->para($noeud) ne "0")
    {
	if(exists $this->{enfants}{$noeud})
	    {
		push @tabTemp,$noeud;
		return \@tabTemp;
	    }
    }
    else{
	foreach my $n (@{$this->{enfants}{$noeud}})
	{
	    push @tabTemp,@{$this->fusionnePara($n)};
	}
    }
    return \@tabTemp;
}

# Renvoie le maximum du tableau en argument
# Methode Statique
sub max{
    my ($tab) = @_;

    my $temp = $tab->[0];

    foreach my $v (@$tab){
	if($v > $temp){
	    $temp = $v;
	}
    }
    return $temp;
}


1;
