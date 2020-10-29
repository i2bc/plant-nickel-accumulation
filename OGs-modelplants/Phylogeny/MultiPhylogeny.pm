package MultiPhylogeny;

# Threads de my_2_phylogeny_launch_jobs.pl 
# Argument 1: dossier contenant les fichiers d'homologies
# Argument 2: stringence (%age de l'alignement)
# Argument 3: Dossier de sortie: dans lequel il y aura les fichiers contenant les infos sur chaque groupe: les proteines du groupe et les distances PAM
# command line: system("perl -I $dirLib 2_phylogeny_job.pl $dirInComparison 70  $dirOutFamily $dirGenomeNew)
sub parse_trifile
{
use strict;
my @arg = @_;
my $workdir = $arg[0]."/Phylogeny/";
my $dirLib = $workdir.'/Phylogeny/lib/'; # INPUT : repertoire avec les librairies utilisees par Phylogeny
my $tag =$arg[2]; # tag des repertoires d'entrée et de sortie input=dataphylo output=groupDist
my $inDir = $workdir.'/Phylogeny/dataphylo_'.$tag; # INPUT : repertoire avec les fichier .tri 
#example intput : /dataPhylo_1e3_0/' (fichiers des résultats du blast triés) 
my $outdir= $workdir.'/groupDist_'.$tag.'/'; # OUTPUT fichiers de graphe
unless ( -e $outdir) { mkdir($outdir);}
my $stringence = $arg[1]; # Pourcentage minimum d'identité entre les protéines
my $targetfile = $arg[3]; # repertoire contenant les fichiers fasta de sequence

# clé: pid, valeur: groupe
my $groupe = $arg[4]; # reference %groupe
my $pid = $arg[5]; # reference %pid
# Clé: groupe, valeur: liste de protéines
my %pid;
# Clé : pid : valeur: espece
my %pidSpecies;
# 2 clés: pid1-pd2, valeur: dPAM
# distance PAM pour chaque paire de protéine
my $dPam;

my $prot_1;
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

chdir $workdir;
my $pourcentage;

open(IN, $targetfile) or die ("impossible d'ouvrir $targetfile");
# On parcourt un fichier .tri ligne par ligne.
# Une ligne correspond à une comparaison prot-prot pour une certaine comparaison génome-génome
while(<IN>)
{
	    chomp;
		# pour le champ query, recuperation de l abrev pour l espece (jusqu'au premier _ ) dans spec_1 et du nom de la proteine dans prot_1, remplacement des underscoreet des points des noms par des %U et %P 
	    my @temp = split(/\t/,$_);
		if($temp[0]=~/^(.+?)_(.+)$/)
		{
			$spec_1 = $1;
	       	$prot_1 = $2;
            $prot_1 =~ s/\_/%U/g;
            $prot_1 =~ s/\./%P/g;
        }
        else
        {
        	die ("Erreur with ".$targetfile."\n");
        }
        # même topo pour le champ hit 
	    if($temp[1]=~/(.+?)_(.+)$/)
		{      
				$spec_2 = $1;
				$prot_2 = $2;
				$prot_2 =~ s/\_/%U/g;
				$prot_2 =~ s/\./%P/g;
		}
		else
		{
			die ("Erreur with ".$targetfile."\n");
		}
	    # Calcul du pourcentage de l'alignement
	    $pourcentage = min($temp[11]/$temp[6],$temp[12]/$temp[7])*100;

	    # On verifie qu'on est bien au dessus du seuil
	    if($pourcentage>$stringence)
	    {
			$dPam->{$prot_1}{$prot_2} = $temp[13];
			# On stocke l'espèce pour les pids en question
			$pidSpecies{$prot_1} = $spec_1;
			$pidSpecies{$prot_2} = $spec_2;
			# Si la première protéine forme déjà un groupe et la deuxième aussi
			if(defined $groupe->{$prot_1} && defined $groupe->{$prot_2})
			{
				# Et si ce sont bien deux groupes distinctes
				if($groupe->{$prot_1} ne $groupe->{$prot_2})
				{
				# On fusionne les deux groupes
				fusionnerGroupes($groupe->{$prot_1},$groupe->{$prot_2},$groupe,$pid);
				}
			}
			# Si une des deux protéines n'appartient pas à un groupe
			else
			{
				# Si la première appartient à un groupe
				if(defined $groupe->{$prot_1} && !defined $groupe->{$prot_2})
				{
					# On ajoute la 2èmè protéine au groupe contenant la 1ère protéine, on passe en argument le hash
					ajouterProt($prot_2,$groupe->{$prot_1},$groupe,$pid);
				}
				# Si ce n'est pas la première qui appartientà un groupe
				else 
				{
					# mais que c'est la seconde
					if(!defined $groupe->{$prot_1} && defined $groupe->{$prot_2})
					{
			    		# On ajoute la 1èrè protéine au groupe contenant 
						# la 2ème protéine
						ajouterProt($prot_1,$groupe->{$prot_2},$groupe,$pid);
					}
					# Si ce n'est pas non plus la seconde prot qui
					#appartient à un groupe
					else
					{
			    		# on cré un nouveau groupe contenant les 2 protéines
			    		creerGroupe($prot_1,$prot_2,$groupe,$pid);
			    		# l'ordre des conditions logiques est important, pour qu'on ne cré pas un groupe à chaque fois		
					}
				}
			}
	    }
	}
	close(IN);
}   
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
    

# renvoie le minimum des deux valeurs en argument
sub min
{
    my ($v1,$v2) = @_;
    ($v1<$v2)?return $v1: return $v2;
}
1;
