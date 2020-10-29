#!/usr/bin/perl -w

# - Lecture des fichiers de groupe ( ?) 
# - Pour chaque groupe:
# --- extraction des sequences des proteines
# --- Alignement avec Muscle (version locale)
# --- Construction d'un arbre phylogenetique du groupe avec PhyML (version locale)
# --- Modification de la racine de l'arbre avec fretree du package embassy-Phylip (=dependance du script) [etape modifie par CD voir bug avec retree_modif]

# exec: repertoire contenant les executables phyml, muscle, et retree
# Argument 1 : ou sont les executables (a l'exeption de fretree installé via synaptic par CD)
# Argument 2: repertoire de sortie des binaires
# Argument 3: repertoire de sortie du script
# Argument 4: fichier contenant le graphe du groupe en input (tache courrante)
# argument 5: repertoire des proteomes fasta

use strict;
if($#ARGV ne 4)
{
    print "4_phylogeny_job1.pl <repertoire des binaires> <repertoire de sortie des binaires> <repertoire des genomes fasta> \n";
    exit();
}
# command: system('perl '.$dirGENERAL."/my4_phylogeny_job1.pl $dirExec $dirTmp $dirTree $dirGroup"."$file ".$dirData);

# Repertoires de travail
my $binDir  = $ARGV[0];   #print "bindir $binDir\n";#exec dans l'exemple de CP 
my $workDir = $ARGV[1];   #print "workdir $workDir\n";
my $outDir  = $ARGV[2];   #print "outdir $outDir\n";  #arbre_1e3_70 dans l'exemple de CP
my $string = $ARGV[3];    #print "string $string\n"; #groupDist_1e3_70/group_numdugroup
my $fasta = $ARGV[4];     #print "fasta $fasta\n"; #genome fasta

# ensemble des pids du groupe
my @pid; 
# les especes auxquelles appartiennent les proteines
my %species;
# Nom du fichier qui contiendra les sequences des proteines du groupe (format fasta)
my $fichierFastaTemp = $workDir."temp.fa";
# Nom du fichier qui contiendra l'alignement du groupe (output de muscle)
my $fichierMuTemp = $workDir."temp.mu";
# Les pids temporaires assignes a chaque proteine (-> pour eviter les pid de + de 10 caracteres incompatible avec Phylip)
my %pidTemp;
# nom du groupe de proteine qu'on analyse
my $groupe;

# parsing du graphe du groupe (fichier $string),  documentation de @pid et %species
parseGenomes($string,\@pid,\%species,\$groupe);

#exception: si il y a moins de 2 proteines dans groupe, on ecrit l'arbre sans calcul
if($#pid<2)
{
    open(OUT,">".$outDir."arbre_".$groupe.".txt ");
    print OUT "(".$pid[0].":0,".$pid[1].":0);\n";
    close(OUT);
    exit(0);
}
#exception: si il y a plus de 1000 proteines dans le groupe, pas d'analyse
if($#pid>1000)
{
    print STDERR "Groupe ".$groupe." non analyse (".$#pid." proteines)\n";
    exit(0);
}

$fichierFastaTemp .="_".$groupe;
$fichierMuTemp .="_".$groupe;

print "fasta\n";
# Ecriture des sequences des proteines du groupe dans un fichier FASTA et stockage des  pids Temporaires dans la variable %pidTemp
# Les proteines extraites sont numerotés a cette etape (a partir de 1, par ordre d'apparition dans la ligne nodes)
ecrireFichierFasta($fichierFastaTemp,\@pid,\%species);

print "align\n";
# execution de MUSCLE
my $align = "".$binDir."muscle3.6/muscle -quiet -stable -in $fichierFastaTemp -out $fichierMuTemp";
if($#pid>500){$align.=" -maxiters 2 ";}
system $align;

# Conversion du fichier de sortie de MUSCLE au format phylip (pour input dans PhyML)
print "Fichier: ".$fichierMuTemp."\n";
fasta2phylip($fichierMuTemp);

print "phylip\n";

# Execution de PhyML
system "".$binDir."phyml_v2.4.4/exe/phyml_linux $fichierMuTemp 1 i 1 0 WAG 0.0 4 1.0 BIONJ n n";

# Suppression des fichiers temporaires
unlink($fichierMuTemp."_phyml_stat.txt");
unlink($fichierMuTemp."_phyml_lk.txt");
#unlink($fichierMuTemp);
#unlink($fichierFastaTemp);
#unlink("log");

# Deplacement de la racine de l'arbre PhyML au midpoint (phylip retree modifie pour accepter les arguments de la ligne de commande = retree_modif)
# CD: remplacement de retree_modif par fretree de embassy-phylip (pas de sources pour retree_modif, recompilation impossible)
print "Retree\n ";
unlink($outDir."arbre_".$groupe.".txt") if (-e $outDir."arbre_".$groupe.".txt");
open(OUT,"| fretree -intreefile ".$fichierMuTemp."_phyml_tree.txt -outtreefile ".$outDir."arbre_".$groupe.".txt >/dev/null");
print OUT "\n";
print OUT "M\n";
print OUT "W\n";
print OUT "R\n";
print OUT "Q\n";
close(OUT);

# retablissement des noms initiaux des proteines
my $outfile = "STEP4-1/arbre_$groupe.txt";
my $filecontent = back2name ($outfile,\@pid);
# remplacement du fichier numeroté par le fichier avec les noms initiaux
print $filecontent;
open(OUT,">$outfile");
print OUT $filecontent;
close OUT;

# suppression du fichier arbre temporaire
# unlink($fichierMuTemp."_phyml_tree.txt");
# Suppression du fichier d'entree contenant les proteines du groupe
# unlink($workDir.$string);

##################################################
#                                                #
#         Definition des Fonctions utilisees     #
#                                                #
##################################################

# $file (= $string): arbre courrant a traiter
# $pid: ref vers une liste des proteines de l'arbre
sub parseGenomes
{
    # $file = $string: fichier à parser
    # $pid: ref vers @: oe stocker les donnees
    my ($file,$pid,$species,$groupe) = @_;
    open(IN,$file);
    # On recupere le No du groupe
    my @f = split(/\_|\..*/,$file);
    $$groupe = $f[$#f];
    
    while(<IN>)
    {
		chomp;
		if(/\(nodes/) # extraction du nom des proteines de la premier ligne du fichier 
		# ex: (nodes 7_3700 3_4270 8_4148 0_1843 6_2019 1_4828 9_4149 6_2898 5_3556 4_4060 2_2875 2_4622 4_4418 )
		{
			print $_."\n";
			s/nodes|\(|\)//g;
			s/^\s|\s$//g;
			my @temp = split(/\s/); # @temp est la liste des proteines
			foreach my $p (@temp)
			{ 
				#pour chaque prot (ID = <numero du genome>_<numero de la proteine>)
				my @p = split(/\_/,$p);
				push (@$pid,$p) ; # construction d'une liste de proteines				
				$species->{$p} = $p[0]; # association proteine->genome
			}
		}
    }
    close(IN);
}

# A partir des pids /especes
# Ecrit un fichier fasta contenant les sequences des proteines du graphe sachant leur espece et  pid 
sub ecrireFichierFasta
{
    my ($fastafile,$pid, $species) = @_;
    my @proteins = @{$pid};
    open(OUT,">".$fastafile) or die "Impossible d'ouvrir le fichier $fastafile";

    # Pour chaque pid on retrouve sa sequence et son espece dans la base de donnees
    # On met les infos dans un fichier Fasta pour donner e muscle->phyml
    for (my $i = 0; $i< scalar(@proteins) ; $i++)
    {	
		my $match = 0;
		my $seq;	
        my $species = $species->{$proteins[$i]}; print $species."\n";
		if (-e $fasta."/".$species.".fa")
		{
			open(IN,$fasta."/".$species.".fa");
		}
		else
		{
			die ("Impossible d'ouvrir ".$species.".fa");
		}
	
   	   	while (my $l = <IN>)
   	   	{	
			if ($l =~ /^\>(.*)$/ and $match == 1)
			{
				last;
			}
			elsif ($l =~ /^\>(.*)$/ and $match == 0)
			{
				chomp($l);
				my $id = $1;
				if ($id  eq $proteins[$i])
				{
					$match = 1;				
					$seq = ">".$species.'_'.($i+1)."\n";
				}
			}	
			else 
			{
				if ($match==1)
				{				
						$seq .= $l;
				}
			}								
		}
		$seq =~ s/\*$//;
		print OUT $seq;	
		close(IN);						
	}
	close(OUT);	
	
}

# Convertit le fichier infile du format fasta au format phylip
sub fasta2phylip
{
    my ($infile)=@_;
    
    my $nom;
    my @tab=();
    my %h1=();
    my $nbseq=0;
    
    open IN,"$infile";
    while(<IN>)
    {
	chomp;
	#recherche ">" comme le debut du nom d'une sequence
	if(/^>/)
	{
	    #enleve le signe ">" dans le nom
	    s/>//;
	    $nom=$_;
	    #le scalaire $nom2 devient la cle du hash
	    $h1{$nom}="";
	    #cree un tableau correspondant e la liste des noms
	    push(@tab,$nom);
	}
	else
	{
	    s/\s//g;
	    #concatene les lignes de sequence 
	    $h1{$nom}.=$_;
	}
    }
    close IN;
    $nbseq=@tab;

    my $outfile=$infile;    
    open OUT ,"> $outfile";	
    my $longueur=length($h1{$tab[0]});
    #indique au debut du fichier le nombre de sequences alignees et la longueur de ces sequences
    print OUT "\t$nbseq\t$longueur\n";
    #initialisation d'un curseur de position sur la sequence
    my $pos = 0;
    #execution d'une boucle tant que le curseur n'a pas atteint le dernier element de la sequence
    while($pos<=$longueur)
    {
	#parcourt les cles du hash
	foreach(@tab)
	{
	    if ($pos < 60)
	    {
		#affiche l'un en face de l'autre le nom de la sequence et la partie de sequence 
		#extraite et retourne e la ligne pour afficher les sequences l'une en dessous de l'autre
		print OUT substr($_,0,10)."\t\t".substr($h1{$_},$pos,60)."\n";
	    }
	    else
	    {
		#affiche l'un en face de l'autre le nom de la sequence et la partie de sequence
		#extraite et retourne e la ligne pour afficher les sequences l'une en dessous de l'autre
		print OUT "\t\t\t".substr($h1{$_},$pos,60)."\n";
	    }
	}
	#ajoute 60 au curseur pour que le foreach puisse afficher le reste de la sequence
	$pos+=60;
	#saute une ligne entre chaque groupe de sequence
	print OUT "\n";
    }
    close OUT;
    return($outfile);
}

# Retablissement des noms initiaux des proteines
sub back2name
{
	
	my ($outfile, $pid) = @_;
	my @names = @{$pid};
	for ( my $i = 0; $i<scalar(@names); $i++) { print "$i is $names[$i]\n";}
	
	# retablissement des noms initiaux des proteines	
	open (IN, $outfile);
	
	my $rev = `cat $outfile`;
	print $rev;
	$rev =~ s/\n//g;
	my @rev = split(//,$rev);
	for ( my $i = 0; $i<scalar(@rev); $i++) { print " char $i is $rev[$i]\n";}
	my $arbre = $rev[0]; # variable qui stocke le contenu du fichier apres reversion des noms des proteines, le premier caractere sera '('
	my $string; # variables qui strocke les nom des proteines ou les distances
	for (my $i = 1; $i <= $#rev; $i++)
	{
		#print $rev[$i]."\n";
		# case: suite de caractere specifique du fichier newick 
		# on ajoute le caractere a l'arbre
		if ($rev[$i] =~ /[(),:;]/ and $rev[$i-1] =~ /[(),:]/)
		{
			$arbre .= $rev[$i]; 
		}
		# caractere de transition entre les caracteres newick et les noms de proteines ou les distances
		# on ajoute la distance a l'arbre
		elsif ( $rev[$i] =~ /[(),:]/ and $rev[$i-1] =~ /[0-9_]/) 
		{	
			#print $string."\n";
			if($string !~ /_/) # ce n'est pas un nom de proteine (mais une distance), ecriture  du string dans l'arbre 
			{
				$arbre .= $string;	#ajout de la distance			
				$arbre .= $rev[$i];	#ajout du caractere special courrant
				$string = ''; # re-initialisation de la variable string 
			}
			else # c'est un nom de proteine 
			{	
				my @p = split ('_',$string); #remplacement du nom de la proteine par son nom initial
				my $rank = $p[1]; 
				#print "rank $rank \n";	
				my $j =  int($rank) -1; 
				my $name = $names[$j];
				#print "name is $name \n";
				$arbre .= $name;	#ecriture du nom corrigé dans le fichier			
				$arbre .= $rev[$i]; # ecriture du caractere special courrant
				$string = ''; # reinitialisation de la variable string
			}
		}
		else
		{
			$string .= $rev[$i];
		}
	}
	close IN;
	#print $arbre;
	return $arbre;	
}

# Cette fonction ouvre le fichier d'arbre ($in), et remet les veritables pids des proteines
# Donc de taille > 10 caracteres
# CD: Ne fonctionne pas avec la numerotation préalable des proteines mise en place pour OrthoFinder
# remplacement par back2name
sub retablirPID
{
    my ($in,$pidTemp) = @_;
    my $total="";
	print "## in $in\n";
    open(IN,$in) or die("impossible d'ouvrir $in");
    while(<IN>){
		print "## _ $_\n";
		#(4932_sacch:0.79326,(4952_yarro:0.62467,((5270_ustil:0.63695,5207_crypt:0.61792):0.09315,
		
		#s/[\d]+_([0-9]*):/$pidTemp->{$1}:/g; #pour genome non modifie, abrev = 5 lettres
		s/[a-zA-Z]+?_([0-9]*):/$pidTemp->{$1}:/g;
		
		#s/([a-zA-Z]{5}2_)([0-9]*)/$1$pidTemp->{$2}/g; #pour genome modifie, abrev = 5 lettres + "2"
		$total.=$_;
		print "\ntotal $total\n";
    }
    close(IN);

    open(OUT,">".$in);
    print OUT $total;
    close(OUT);
}
