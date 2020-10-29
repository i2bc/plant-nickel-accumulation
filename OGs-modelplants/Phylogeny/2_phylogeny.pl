#!/usr/bin/perl -w
# Ce script permet de construire l'ensemble des graphes d'homologie entre les proteines. 
# Chaque graphe correspond a  un groupe de proteines, et est stocke dans un fichier
# Le dossier de sortie contient un fichier par groupe: les proteines du groupe et les distances PAM (?)

# Argument 1: stringence (%age de l'alignement)
# Argument 2: dossier "git clone"
# Argument 3: tag du run (etiquette les dossiers successifs)
# Argument 4 (Ajouter par Sandrine): repertoire contenant les genomes ==> Pour recuperer leur code
use strict;
use lib::Elague;
print "2 phylogeny job \n";

if($#ARGV<2)
{
	print "2_phylogeny_job.pl <stringence> <dossier de travail> <dossier des genomes fasta>\n";
	exit();
}

# Pourcentage minimum d'identite entre les protines
my $stringence = $ARGV[0];
# repertoire de travail (Phylogeny)
my $workdir = $ARGV[1].'/Phylogeny';
my $fastadir = $ARGV[2];
# Dossier d'entree
my $inDir = $workdir.'/STEP1/';
my $outdir = $workdir.'/STEP2/';
unless ( -e $outdir) { mkdir($outdir);}

# Tableau d'especes
my @species;
# cle: pid, valeur: groupe
my %groupe;
# Cle: groupe, valeur: liste de protines
my %pid;
# Cle: pid : valeur: espece
my %pidSpecies;
# 2 cles: pid1-pid2, valeur: dPAM
# distance PAM pour chaque paire de proteine
my $dPam;

# Numro provisoire des PIds (indice dans un tableau)
my @numProvPid;

#Nb total de groupes;
my $nbGroupes =0;

# Nb total de protines
my $nbProt =0;



my @File = ();
my $file = '';
if (-e $fastadir){
	chdir($fastadir);
	@File = glob("*.fa");
	foreach $file (@File) {
		$file =~ s/\.fa$//;
		push(@species,$file);
	}
}

# On trie le tableau d'especes pour les prendre dans le bon ordre
@species = sort(@species);

# Va contenir les donnes de chaque lignes du fichier, analyse du
# fichier ligne par ligne
my @temp;
my $prot_1;
my $prot_2;
my $spec_1;
my $spec_2;

chdir $workdir;
my $pourcentage;

# On parcour les especes
for(my $i=0;$i<=$#species;$i++)
{
    # On parcourt chaque paires
    for(my $j=$i;$j<=$#species;$j++)
    {
		print $inDir.$species[$i].'_vs_'.$species[$j].".tri\n";
		# Fichier de log
		open(OOO,">>".$outdir."/log");
        # Fichier contenant les comparaisons genome-genome
        if (-e $inDir.$species[$i].'_vs_'.$species[$j].".tri")
        {
			print OOO "$i - $j Opening ".$inDir.$species[$i].'_vs_'.$species[$j].".tri\n";
            open(IN,$inDir.$species[$i].'_vs_'.$species[$j].".tri") or die ("impossible d'ouvrir $inDir$species[$i]_vs_$species[$j].tri");
        }
        elsif( -e $inDir.$species[$j].'_vs_'.$species[$i].".tri")
        {
           print OOO "$i - $j Opening ".$inDir.$species[$j].'_vs_'.$species[$i].".tri\n";
           open(IN,$inDir.$species[$j].'_vs_'.$species[$i].".tri") or die ("impossible d'ouvrir $inDir$species[$j]_vs_$species[$i].tri");
        }
        else
        {
            die("Erreur: no file for $species[$i] and $species[$j]\n");
        }
		close(OOO);

		# On le parcourt ligne par ligne.
		# Une ligne correspopnd a  une comparaison prot-prot pour une
		# certaine comparaison genome-genome
		while(<IN>)
		{
			chomp;	    
			@temp = split(/\t/,$_);
			if($temp[0] =~ /^(.+?)_(.+)$/)
			{
				$spec_1 = $1; #print $spec_1."\n";
	       		$prot_1 = $temp[0];	#print $prot_1."\n";			
            }
            else
            {
                die ("Erreur with ".$inDir.$species[$i].'_vs_'.$species[$j].".tri ($temp[0])\n");
            }
            
			if($temp[1] =~ /(.+?)_(.+)$/)
			{
				$spec_2 = $1;
	       		$prot_2 = $temp[1];
            }
            else
            {
                die ("Erreur with ".$inDir.$species[$i].'_vs_'.$species[$j].".tri ($temp[1])\n");
            }
            
			# Calcul du pourcentage de l'alignement
			$pourcentage = min($temp[11]/$temp[6],$temp[12]/$temp[7])*100;

			# On verifie qu'on est bien au dessus du seuil
			#if($pourcentage>$stringence and $temp[13]<250)
			if($pourcentage>$stringence)
			{
				# On stocke la distance PAM entre les deux protines
				$dPam->{$prot_1}{$prot_2} = $temp[13];
				# On stocke l'espece pour les pids en question
				$pidSpecies{$prot_1} = $spec_1;
				$pidSpecies{$prot_2} = $spec_2;
				# Si la premiere proteine forme dja  un groupe et la
				# deuxieme aussi
				if(defined $groupe{$prot_1} && defined $groupe{$prot_2})
				{
					# Et si ce sont bien deux groupes distinctes
					if($groupe{$prot_1} ne $groupe{$prot_2})
					{
						# On fusionne les deux groupes
						fusionnerGroupes($groupe{$prot_1},$groupe{$prot_2},\%groupe,\%pid);
					}
				}
				# Si une des deux protines n'appartient pas a  un groupe
				else
				{
					# Si la premiere appartient a  un groupe
					if(defined $groupe{$prot_1} && !defined $groupe{$prot_2})
					{
						# On ajoute la 2eme proteine au groupe contenant
						# la 1ere protine, on passe en argument le hash
						ajouterProt($prot_2,$groupe{$prot_1},\%groupe,\%pid);
					}
					# Si ce n'est pas la premiere qui appartienta  un groupe
					else 
					{
						# mais que c'est la seconde
						if(!defined $groupe{$prot_1} && defined $groupe{$prot_2})
						{
							# On ajoute la 1ere protine au groupe contenant 
							# la 2eme protine
							ajouterProt($prot_1,$groupe{$prot_2},\%groupe,\%pid);
						}
						# Si ce n'est pas non plus la seconde prot qui
						#appartient a  un groupe
						else
						{
							# on cree un nouveau groupe contenant les 2 protines
							creerGroupe($prot_1,$prot_2,\%groupe,\%pid);
							# l'ordre des conditions logiques est
							# important, pour qu'on ne crepas un
							# groupe a  chaque fois
		
						}
					}
				}
			}
		}
		close(IN);
    }
}
#afficher(\%pid,\%groupe);

# On ecrit l'ensemble des fichier de groupes ainsi calculs
ecrireFichiers(\%pid,\%groupe,$outdir);


#statFam(\%pid,\%groupe);

# Fonction qui fusionne deux groupes
sub fusionnerGroupes
{
    my ($groupe_1,$groupe_2,$groupe,$pid) = @_;

    foreach my $prot (@{$pid->{$groupe_2}})
    {
		$groupe->{$prot} = $groupe_1;
    }
    push @{$pid->{$groupe_1}},@{$pid->{$groupe_2}};
    delete $pid->{$groupe_2};
}

# Ajoute une protine prot au groupe gr
sub ajouterProt
{
    my ($prot,$gr,$groupe,$pid) = @_;
    
    $groupe->{$prot} = $gr;
    push(@{$pid->{$gr}},$prot);
}

# Cre un nouveau groupe contenant les protines p1 et p2
sub creerGroupe
{
    my ($p1,$p2,$groupe,$pid) = @_;
    $nbGroupes++;
    $groupe->{$p1} = $nbGroupes;
    $groupe->{$p2} = $nbGroupes;
    push @{$pid->{$nbGroupes}},$p1;
    push @{$pid->{$nbGroupes}},$p2;
}

# Affiche les groupes.... Inutilisable pour un tres grand nombre de groupes.
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

# Ecris les fichiers de groupes (format de graphe ressemblant a  celui
# de Tulip)
sub ecrireFichiers
{
    my ($pid,$groupe,$outdir) = @_;
    foreach my $g (keys(%$pid))
    {
		# On ne garde que les edges qui sont suprieures au seuil dfini plus haut. 
		# Donc, le graphe n'est pas complet
		# Lorsque on analyse les arbres apres, on cherche les distances pam qui nous manquent!
		# Voir analayseArbres_all!
		#my $dPam = dPam($pid,$groupe,$g);
		open(OUT,">".$outdir."/groupe_".$g) or die "Impossible d'ouvrir le fichier ".$outdir."/groupe_".$g;
		print OUT "(nodes ";
		# On ecrit les noeuds du graphe
		foreach my $p (@{$pid->{$g}})
		{
			print OUT $p." ";
		}
		print OUT ")\n";
		# On ecrit les liens
		foreach my $p1 (@{$pid->{$g}})
		{
			foreach my $p2 (@{$pid->{$g}})
			{
				if (defined $dPam->{$p1}{$p2})
				{
					print OUT "(edge ".$p1." ".$p2." ".$dPam->{$p1}{$p2}.")\n";
					delete $dPam->{$p1}{$p2};
				}
				if (defined $dPam->{$p2}{$p1})
				{
					print OUT "(edge ".$p1." ".$p2." ".$dPam->{$p2}{$p1}.")\n";
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
			foreach my $p1 (@{$pid->{$g}})
			{
				print $p1.', ';
			}	
			print "\n";
		}
	
		# On casse les ponts dans le groupe
		if($#{$pid->{$g}}>3 && $#{$pid->{$g}}<30000) 
		{
			Elague::elague($outdir."/groupe_".$g);
		}
		else
		{
			if($#{$pid->{$g}}>=30000)
			{
		    	print STDERR "Skipping Elaguage of $g (".$#{$pid->{$g}}."proteins)";
			}
		}
		delete $pid->{$g};
		#undef $dPam;
    }
}

# Va rechercher la distance pam entre la protine $pid a l'interieur
# du groupe $groupe en argument dans les fichiers tris. Fonction tres
# lente!!! Pas utilise pour le moment
sub dPam
{
    my ($pid,$groupe,$g) = @_;
    my $fichier;
    # Dit si on a trouveles pids dans le fichier d'homologie (pour avoir la distance)

    my %species;
    my $dpam = {};

    foreach my $p (@{$pid->{$g}})
	{
		$species{$pidSpecies{$p}} = 1;
    }
    
    foreach my $s1 (keys(%species))
	{
		foreach my $s2 (keys(%species))
		{
	    	if($s1 lt $s2)
			{
	    		$fichier = $inDir."/".$s1.'_vs_'.$s2.".tri";
	    		open(IN,$fichier) or die "Impossible d'ouvrir le fichier $fichier";
				while(<IN>)
				{
		   			 /^(\d*?)\t(\d*?)\t/;
					# print "prot1: ".$1." et prot2: ".$2."\n";
		    		if(defined $groupe{$1} && defined $groupe{$2} && $groupe{$1} eq $g && $groupe{$2} eq $g)
					{
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
