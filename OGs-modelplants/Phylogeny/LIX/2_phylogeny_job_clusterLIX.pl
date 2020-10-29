#!/usr/bin/perl -w
#script recupere
# Ce script permet de construire l'ensemble des graphes d'homologie
# entre les protéines. Chaque graphe correspond à un groupe de
# protéines, et est stocké dans un fichier (dans le dossier fourni en
# argument)

# Argument 1: dossier contenant les fichiers d'homologies
# Argument 2: stringence (%age de l'alignement)
# Argument 3: Dossier de sortie: dans lequel il y aura les fichiers 
#             contenant les infos sur chaque groupe: les proteines du groupe et 
#             les distances PAM

# Argument 4 (Ajouter par Sandrine): repertoire contenant les génomes ==> Pour récupérer nom des génomes
use strict;
use lib::Elague;
print "2 phylogeny job \n";
if($#ARGV<3){
    die "Pas le bon nombre d'arguments";
}

system("mkdir $ARGV[2]") unless (-e $ARGV[2]);

# Tableau d'espèces
my @species;
# clé: pid, valeur: groupe
my %groupe;
# Clé: groupe, valeur: liste de protéines
my %pid;
# Clé : pid : valeur: espece
my %pidSpecies;
# 2 clés: pid1-pd2, valeur: dPAM
# distance PAM pour chaque paire de protéine
my $dPam;

# Numéro provisoire des PIds (indice dans un tableau)
my @numProvPid;

#Nb total de groupes;
my $nbGroupes =0;

# Nb total de protéines
my $nbProt =0;

# Dossier d'entree
my $inDir = $ARGV[0];
# Pourcentage minimum d'identité entre les protéines
my $stringence = $ARGV[1];

if($#ARGV<3){
    exit(-1);
}

my $outdir = $ARGV[2];

my @File = ();
my $file = '';
print @ARGV;
if (-e $ARGV[3]){
	chdir($ARGV[3]);
	print "$ARGV[3]\n";
	@File = glob("*.fa");
	foreach $file (@File) {
		$file =~ s/.fa//;
		push(@species,$file);
		print $file."\n";
	}
}

# On trie le tableau d'espèces pour les prendre dans le bon ordre
@species = sort(@species);

# Va contenir les données de chaque lignes du fichier, analyse du
# fichier ligne par ligne
my @temp;
# Inutilisée...
my $prot_1;
# Inutilisée...
my $prot_2;
my $spec_1;
my $spec_2;
# Inutilisée...
#my $longueurProt_1;
# Inutilisée...
#my $longueurProt_2;
# Inutilisée...
#my $longueurAlign_1;
# Inutilisée...
#my $longueurAlign_2;
# Donne la proportion de l'alignement par rapport à la pls petite
# séquence
my $pourcentage;

# On parcour les espèces
for(my $i=0;$i<=$#species;$i++)
{
    # On parcourt chaque paires
    for(my $j=$i;$j<=$#species;$j++)
    {
	# Fichier de log
	open(OOO,">>".$outdir."log");
        # Fichier contenant les comparaisons génome-génome
        if (-e $inDir.$species[$i].'_'.$species[$j].".tri")
        {
	   print OOO "$i - $j Opening ".$inDir.$species[$i].'_'.$species[$j].".tri\n";
           open(IN,$inDir.$species[$i].'_'.$species[$j].".tri");
        }
        elsif($inDir.$species[$j].'_'.$species[$i].".tri")
        {
           print OOO "$i - $j Opening ".$inDir.$species[$j].'_'.$species[$i].".tri\n";
           open(IN,$inDir.$species[$j].'_'.$species[$i].".tri");
        }
        else
        {
            die("Erreur: no file for $species[$i] and $species[$j]\n");
        }
	close(OOO);


	# On le parcourt ligne par ligne.
	# Une ligne correspopnd à une comparaison prot-prot pour une
	# certaine comparaison génome-génome
	while(<IN>)
	{
	    chomp;
	    @temp = split(/\t/,$_);

            if ($temp[0] =~ /^([a-zA-Z]{5})2?_(.*)$/)
            {
                $spec_1 =  lcfirst($1);
	       $prot_1 = $2;
               $prot_1 =~ s/\_/%U/g;
               $prot_1 =~ s/\./%P/g;
            }
            else
            {
                die ("Erreur with ".$inDir.$species[$i].'_'.$species[$j].".tri ($temp[0])\n");
            }
            if ($temp[1] =~ /^([a-zA-Z]{5})2?_(.*)$/)
            {
                $spec_2 = lcfirst($1);
	       $prot_2 = $2;
               $prot_2 =~ s/\_/%U/g;
               $prot_2 =~ s/\./%P/g;
            }
            else
            {
                die ("Erreur with ".$inDir.$species[$i].'_'.$species[$j].".tri ($temp[1])\n");
            }
# 	    my @lignetemp = split(/\_/,$temp[0]);
# 	    $prot_1 = $lignetemp[1];
# 	    @lignetemp = split(/\_/,$temp[1]);
# 	    $prot_2 = $lignetemp[1];
	    #$longueurProt_1 = $temp[6];
	    #$longueurProt_2 = $temp[7];
	    #$longueurAlign_1 = $temp[11];
	    #$longueurAlign_2 = $temp[12];

	    #print $prot_1." ".$prot_2." ".$temp[11]." ".$temp[6]." ".$temp[12]." ".$temp[7]."\n";
	    # Calcul du pourcentage de l'alignement
	    $pourcentage = min($temp[11]/$temp[6],$temp[12]/$temp[7])*100;

	    # On verifie qu'on est bien au dessus du seuil
	    #if($pourcentage>$stringence and $temp[13]<250)
	    if($pourcentage>$stringence)
	    {
		# On stocke la distance PAM entre les deux protéines
		$dPam->{$prot_1}{$prot_2} = $temp[13];
		# On stocke l'espèce pour les pids en question
		$pidSpecies{$prot_1} = $spec_1;
		$pidSpecies{$prot_2} = $spec_2;
		# Si la première protéine forme déjà un groupe et la
		# deuxième aussi
		if(defined $groupe{$prot_1} && defined $groupe{$prot_2})
		{
		    # Et si ce sont bien deux groupes distinctes
		    if($groupe{$prot_1} ne $groupe{$prot_2}){
			# On fusionne les deux groupes
			fusionnerGroupes($groupe{$prot_1},$groupe{$prot_2},\%groupe,\%pid);
		    }
		}
		# Si une des deux protéines n'appartient pas à un groupe
		else
		{
		    # Si la première appartient à un groupe
		    if(defined $groupe{$prot_1} && !defined $groupe{$prot_2})
		    {
			# On ajoute la 2èmè protéine au groupe contenant
			# la 1ère protéine, on passe en argument le hash
			ajouterProt($prot_2,$groupe{$prot_1},\%groupe,\%pid);
		    }
		    # Si ce n'est pas la première qui appartientà un groupe
		    else {
			# mais que c'est la seconde
			if(!defined $groupe{$prot_1} && defined $groupe{$prot_2})
			{
			    # On ajoute la 1èrè protéine au groupe contenant 
			    # la 2ème protéine
			    ajouterProt($prot_1,$groupe{$prot_2},\%groupe,\%pid);
			}
			# Si ce n'est pas non plus la seconde prot qui
			#appartient à un groupe
			else
			{
			    # on cré un nouveau groupe contenant les 2 protéines
			    creerGroupe($prot_1,$prot_2,\%groupe,\%pid);
			    # l'ordre des conditions logiques est
			    # important, pour qu'on ne cré pas un
			    # groupe à chaque fois
		
			}
		    }
		}
	    }
	}
	close(IN);
    }
}
#afficher(\%pid,\%groupe);

# On ecrit l'ensemble des fichier de groupes ainsi calculés
ecrireFichiers(\%pid,\%groupe,$outdir);


#statFam(\%pid,\%groupe);

# Fonction qui fusionne deux groupes
sub fusionnerGroupes
{
    my ($groupe_1,$groupe_2,$groupe,$pid) = @_;

    foreach my $prot (@{$pid->{$groupe_2}}){
	$groupe->{$prot} = $groupe_1;
    }

    push @{$pid->{$groupe_1}},@{$pid->{$groupe_2}};
    delete $pid->{$groupe_2};
}

# Ajoute une protéine prot au groupe gr
sub ajouterProt
{
    my ($prot,$gr,$groupe,$pid) = @_;
    
    $groupe->{$prot} = $gr;
    push(@{$pid->{$gr}},$prot);
}

# Crée un nouveau groupe contenant les protéines p1 et p2
sub creerGroupe
{
    my ($p1,$p2,$groupe,$pid) = @_;
    $nbGroupes++;
    $groupe->{$p1} = $nbGroupes;
    $groupe->{$p2} = $nbGroupes;
    push @{$pid->{$nbGroupes}},$p1;
    push @{$pid->{$nbGroupes}},$p2;
}

# Affiche les groupes.... Inutilisable pour un très grand nombre de groupes.
sub afficher
{
    my ($pid,$groupe) = @_;
 
   # open(OUT,">/home/lemoine/detect_ortho/total.txt");
    foreach my $g (keys(%$pid))
    {
	print "<groupe id=\"$g\">\n";
	foreach my $p (@{$pid->{$g}})
	{
	    print "<protein pid=\"$p\" species=\"$pidSpecies{$p}\"></pid>\n";
	}
	print "</groupe>\n";
    }
    print "\n";
    
#   close(OUT);
    #foreach my $p( keys(%$groupe))
    #{
    #  print "$p: ".$groupe->{$p}."\n";
    #}
}

# Ecris les fichiers de groupes (format de graphe ressemblant à celui
# de Tulip)
sub ecrireFichiers{
    my ($pid,$groupe,$outdir) = @_;

    foreach my $g (keys(%$pid))
    {
	# On ne garde que les edges qui sont supérieures au seuil défini plus haut. 
	# Donc, le graphe n'est pas complet
	# Lorsque on analyse les arbres après, on cherche les distances pam qui nous manquent!
	# Voir analayseArbres_all!
	#my $dPam = dPam($pid,$groupe,$g);
	open(OUT,">".$outdir."/groupe_".$g) or die "Impossible d'ouvrir le fichier ".$outdir."/groupe_".$g;
	print OUT "(nodes ";
	# On ecrit les noeuds du graphe
	foreach my $p (@{$pid->{$g}}){
	    print OUT $pidSpecies{$p}."_".$p." ";
	}
	print OUT ")\n";
	# On ecrit les liens
	foreach my $p1 (@{$pid->{$g}}){
	    foreach my $p2 (@{$pid->{$g}}){
		if (defined $dPam->{$p1}{$p2}){
		    print OUT "(edge ".$pidSpecies{$p1}."_$p1 ".$pidSpecies{$p2}."_$p2 ".$dPam->{$p1}{$p2}.")\n";
		    delete $dPam->{$p1}{$p2};
		}
		if (defined $dPam->{$p2}{$p1}){
		    print OUT "(edge ".$pidSpecies{$p1}."_$p1 ".$pidSpecies{$p2}."_$p2 ".$dPam->{$p2}{$p1}.")\n";
		    delete $dPam->{$p2}{$p1};
		}
		delete $groupe{$p2};
	    }
	    delete $groupe{$p1};
	}
	close(OUT);

#Ajouter le 9 octobre par sandrine, a supprimer
if (-z 	$outdir."/groupe_".$g)
{
	print $outdir."/groupe_".$g." fichier vide:\n" ;
	foreach my $p1 (@{$pid->{$g}}){
		print $p1.', ';
	}
	print "\n";
}
	# On casse les ponts dans le groupe
	if($#{$pid->{$g}}>3 && $#{$pid->{$g}}<30000){
	    Elague::elague($outdir."/groupe_".$g);
	  }
	else{
		if($#{$pid->{$g}}>=30000){
		    print STDERR "Skipping Elaguage of $g (".$#{$pid->{$g}}."proteins)";
		}
	}

	delete $pid->{$g};
	#undef $dPam;
    }
}


# Va rechercher la distance pam entre la protéine $pid à l'intérieur
# du groupe $groupe en argument dans les fichiers tris. Fonction très
# lente!!! Pas utilisée pour le moment
sub dPam{
    my ($pid,$groupe,$g) = @_;
    my $fichier;
    # Dit si on a trouvé les pids dans le fichier d'homologie (pour avoir la distance)

    my %species;
    my $dpam = {};

    foreach my $p (@{$pid->{$g}}){
	$species{$pidSpecies{$p}} = 1;
    }
    
    foreach my $s1 (keys(%species)){
	foreach my $s2 (keys(%species)){
	    if($s1 lt $s2){
	    $fichier = $inDir."/".$s1.'_'.$s2.".tri";
	    open(IN,$fichier) or die "Impossible d'ouvrir le fichier $fichier";
		while(<IN>)
		{
		    /^(\d*?)\t(\d*?)\t/;
# print "prot1: ".$1." et prot2: ".$2."\n";
		    if(defined $groupe{$1} && defined $groupe{$2} && $groupe{$1} eq $g && $groupe{$2} eq $g){
			my @temp = split(/\t/);
			$dpam->{$s1."_".$1}{$s2."_".$2} = $temp[13];
			}
		}
		close($fichier);
	    }
	}
    }
    return $dpam;
}


sub statFam
{
    my ($pid,$groupe) = @_;
    my $nbG=0;
    foreach my $g (keys(%$pid))
    {
	$nbG++;
	my $n = $#{$pid->{$g}}+1;
	print $n.",";
	#foreach my $p (@{$pid->{$g}})
	#{
	#    print "\t $p\n";
	#}
    }
    print "\n";
}

# renvoie le minimum des deux valeurs en argument
sub min
{
    my ($v1,$v2) = @_;
    ($v1<$v2)?return $v1: return $v2;
}
