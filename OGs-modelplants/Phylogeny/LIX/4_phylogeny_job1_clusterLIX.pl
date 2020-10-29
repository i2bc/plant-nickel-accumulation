#!/usr/bin/perl -w
#script recupere

# Script qui enchaine les traitements:
# - Lecture des fichiers de groupe
# - Pour chaque groupe:
# --- On récupère les séquences de toutes les protéines
# --- On les aligne (Muscle)
# --- On reconstruit l'arbre phylogénétique

# Argument 1: repertoire contenant les exec phyml, muscle, et retree
# Argument 2: repertoire courant (où seront mis les fichiers temp...) et où est situé le fichier entree_retree
# Argument 3: reprtoire de sortie!
# Argument 4: fichier contenant le groupe à analyser (graphe)
# argument 5: le path des fichiers sgml

use strict;
#use blib;
#use DBI;
#use lib::DSN qw(%DSN);

# S'il n'y a pas le bon nombre d'arguments, on
# s'arr�te
if($#ARGV ne 4)
{
    print "not the right number of arguments\n";
    exit();
}

# R�pertoires de travail
my $binDir  = $ARGV[0];
my $workDir = $ARGV[1];
my $outDir  = $ARGV[2];

# Le fichier d'entree: contenant les proteines d'un seul groupe
# format xml:
# <groupe id="123">
# <protein pid="2122"></protein>
# ...
# </groupe>
my $string = $ARGV[3];

# Dossier contenant les fichies sgml
## my $sgml = $ARGV[4];
my $fasta = $ARGV[4];

# ensemble des pids du groupe
my @pid; 
# les esp�ces auxquelles appartiennent les proteines
my %species;
# Nom du fichier qui contiendra les s�quences des prot�ines pr�c�dentes : format fasta
my $fichierFastaTemp = $workDir."temp.fa";
# Nom du fichier qui contiendra les squences ALIGNEES (muscle)
my $fichierMuTemp = $workDir."temp.mu";

# Les pids temporaires assign�s � chaque proteine
# -> pour �viter les pid de + de 10 caract�res
my %pidTemp;

# nom du groupe de proteine qu'on analyse
my $groupe;

# On analyse le fichier d'entree
# Ca met ce qu'il faut dans le tableau @pid et %species
parseGenomes($string,\@pid,\%species,\$groupe);

#si moins de 2 id dans groupe, on affiche directement l'arbre
if($#pid<2){
    open(OUT,">".$outDir."arbre_".$groupe.".txt ");
    print OUT "(".$species{$pid[0]}."_".$pid[0].":0,".$species{$pid[1]}."_".$pid[1].":0);\n";
    close(OUT);
    exit(0);
}
#si plus de 700 id dans groupe, pas d'anlyse
if($#pid>1000){
    print STDERR "Groupe ".$groupe." non analys� (".$#pid." proteines)\n";
    exit(0);
}

# On personnalise le nom des fichiers en fonction du groupe
# Car sinon quand on lance le script sur le cluster: il peut y avoir des 
# fichiers de m�me nom en m�me temps pour deux groupes diff�rents...
$fichierFastaTemp.="_".$groupe;
$fichierMuTemp.="_".$groupe;

print "fasta\n";
# Ecriture des donn�es dans un fichier FASTA (ecrit les s�quences des
# prot�ines)
# Et �a retient les pids Temporaires
ecrireFichierFasta($fichierFastaTemp,\@pid,\%pidTemp);

foreach my $tmptmp (keys(%pidTemp))
{
	print $tmptmp."=>".$pidTemp{$tmptmp}."\n";
}

print "align\n";

# Alignement multiple
# On met la ligne de commande dans un scalaire
my $align = "".$binDir."muscle3.6/muscle -quiet -stable -in $fichierFastaTemp -out $fichierMuTemp";
if($#pid>500){$align.=" -maxiters 2 ";}
system $align;


# Conversion du fichier au format pylip
print "Fichier: ".$fichierMuTemp."\n";
fasta2phylip($fichierMuTemp);

print "phylip\n";



# Construction de la phylog�nie
system "".$binDir."phyml_v2.4.4/exe/phyml_linux $fichierMuTemp 1 i 1 0 WAG 0.0 4 1.0 BIONJ n n";
print "".$binDir."phyml_v2.4.4/exe/phyml_linux $fichierMuTemp 1 i 1 0 WAG 0.0 4 1.0 BIONJ n n";
# Suppression des fichiers temporaires
#unlink($fichierMuTemp."_phyml_stat.txt");
#unlink($fichierMuTemp."_phyml_lk.txt");
#unlink($fichierMuTemp);
#unlink($fichierFastaTemp);
#unlink("log");

# On r�tablie les vrais PID � partir du fichier d'arbre
#retablirPID($fichierMuTemp."_phyml_tree.txt",\%pidTemp);

# On bouge la racine de l'arbre: midpoint (phylip retree modifi� pour 
# accepter les arguments de la ligne de commande)
print "Retree\n ";
unlink($outDir."arbre_".$groupe.".txt") if (-e $outDir."arbre_".$groupe.".txt");
open(OUT,"|".$binDir."phylip3.66/exe/retree_modif ".$fichierMuTemp."_phyml_tree.txt ".$outDir."arbre_".$groupe.".txt >/dev/null");
print OUT "Y\n";
print OUT "M\n";
print OUT "W\n";
print OUT "R\n";
print OUT "Q\n";
close(OUT);
print $binDir."phylip3.66/exe/retree_modif ".$fichierMuTemp."_phyml_tree.txt ".$outDir."arbre_".$groupe.".txt < ".$workDir."entree_retree >/dev/null";
retablirPID($outDir."arbre_".$groupe.".txt",\%pidTemp);

# Suppression du fichier arbre temporaire
#unlink($fichierMuTemp."_phyml_tree.txt");
# Suppression du fichier d'entree contenant les prot�ines du groupe
#unlink($workDir.$string);

##################################################
#                                                #
#         D�finition des Fonctions utilis�es     #
#                                                #
##################################################

# A partir des pids /esp�ces
# Ecrit un fichier fasta contenant les s�qences de chaques prot�ines
sub ecrireFichierFasta
{
    my ($fichierTemp,$pid,$pidTemp) = @_;
    my @Header = ();
    my $nbPid = 1;
    my $id = "";
    my $specie = '';

    open(OUT,">".$fichierTemp) or die "Impossible d'ouvrir le fichier $fichierTemp";

    # Pour chaque pid on retrouve sa s�quence et son espece dans la base de donn�es
    # On met les infos dans un fichier Fasta pour donner � muscle->phyml
    foreach my $p (@$pid)
    {
        my $idTmp = $p; #ID avec _ et .
	$idTmp =~ s/%U/_/g;
	$idTmp =~ s/%P/./g;         
	my $trouve = 0;
	my $seq;
	# open(IN,$fasta."$species{$p}.fa") or die ("Impossible d'ouvrir $fasta.$species{$p}.fa");
	# On cherche la s�quence correspondant au pid en question
	if (-e $fasta."new/".$species{$p}.".fa"){
		open(IN,$fasta."new/".$species{$p}.".fa");
	}
	elsif(-e $fasta."old/".$species{$p}.".fa"){
		 open(IN,$fasta."old/".$species{$p}.".fa");
	}
	else{
		if(-e $fasta.$species{$p}.".fa"){
			open(IN,$fasta.$species{$p}.".fa");
		}
		else{
			die ("Impossible d'ouvrir $species{$p}.fa");
		}
	}
	

	while($trouve == 0)
	{
   	   	$_ = <IN> || die ("PB with $fasta.$species{$p}.fa $p $trouve ($fichierTemp\n");
# if ($_ = <IN>)
# {

#print $_ if ("$fasta/$species{$p}.fa" eq './data/genome_fa/done/Podan.fa');
		if ($_ =~ /^\>(.*)$/)
		{
                        chomp($_);
# 			@Header = split(/\|/,$_);
			$id = $1;
# 			$id =~ s/>//g;
# # 			$id =~ s/ //g;
# #                         $id = $1;
# #  			$id =~ s/_1$/%P1/ if (($species{$p} eq 'Mycgr2') || (($species{$p} eq 'Phybl2'))); #  a suppreimer ensuite car pb dans length .1 alors que _1 normalement
# 			$id =~ s/_/%U/g;
# 			$id =~ s/\./%P/g;

# print "$id et $p\n";

#print "id ".$id if ($species{$p} eq 'Podan');
		}
		if ($id  eq $idTmp)
		{
			$seq = '';
			while ((!$trouve) && ($_ = <IN>))
			{
				if ($_ =~ /^>/)
				{ 	
					$pidTemp->{$nbPid} = $idTmp;
					$specie = $species{$p};
					print OUT ">".$specie."_".$nbPid."\n";
					$seq =~ s/\*$//;
					print OUT $seq."\n";
					$nbPid++;
					$trouve = 1;
				}
				else
				{
					$seq = $seq.$_;
				}
			}
			if ( ($trouve == 0) && ($seq ne '') ) #fin du fichier atteint donc notre seq derniere seq du fichier
			{
				$pidTemp->{$nbPid} = $idTmp;
				$specie = $species{$p};
				print OUT ">".$specie."_".$nbPid."\n";
				$seq =~ s/\*$//;
				print OUT $seq."\n";
				$nbPid++;
				$trouve = 1;
			}
		}

# }
# else
# {
# print "PB avec $id et $p ";
#     exit();
# 
# }
	}
	close(IN);
# print "$p : ".$_ if ($p =~ /Myg.*35189/);exit() if ($p =~ /Myg.*35189/);
    }
    close(OUT);
}

# analyse le fichier d'entree ($file)
# sub parseGenomes
# {
#     # $file: sting: fichier � parser
#     # $pid: ref vers @: o� stocker les donn�es
#     my ($file,$pid,$species,$groupe) = @_;
# 
#     open(IN,$file);
#     # On r�cup�re le No du groupe
#     my @f = split(/\_|\..*/,$file);
#     $$groupe = $f[$#f];
#     while(<IN>)
#     {
# 	chomp;
# #	if(/\<protein.*pid=\"(.*?)\".*?species=\"(.*?)\">/){
# 	if(/\(nodes/){
# 	    s/nodes|\(|\)//g;
# 	    s/^\s|\s$//g;
# 	    my @temp = split(/\s/);
# 	    #my $_pid = $1;
# 	    #my $_species = $2;
# 	    foreach my $p (@temp){ #pour chaque prot
# 		my @t2 = split(/\_/,$p); #recup du nom de la prot: espece_ID
# 		push @$pid,$t2[1]; #on ajoute l'id
# 		$species->{$t2[1]} = $t2[0]; #on associe � l'id l'espece
# 	    }
# 	}
#     }
#     close(IN);
# }

sub parseGenomes
{
    # $file: sting: fichier � parser
    # $pid: ref vers @: o� stocker les donn�es
    my ($file,$pid,$species,$groupe) = @_;
    my %Id = ();
print "$file\n";
    open(IN,$file);
    # On r�cup�re le No du groupe
    my @f = split(/\_|\..*/,$file);
    $$groupe = $f[$#f];
    while(<IN>)
    {
	chomp;
#	if(/\<protein.*pid=\"(.*?)\".*?species=\"(.*?)\">/){
	if(/\(nodes/){
	    s/nodes|\(|\)//g;
	    s/^\s|\s$//g;
	    my @temp = split(/\s/);
	    #my $_pid = $1;
	    #my $_species = $2;
	    foreach my $p (@temp){ #pour chaque prot
		my @t2 = split(/\_/,$p); #recup du nom de la prot: espece_ID
		push @$pid,$t2[1] unless (exists($Id{$t2[1]})) ; #on ajoute l'id
                $Id{$t2[1]} = '';
		$species->{$t2[1]} = $t2[0]; #on associe � l'id l'espec'
	    }
	}
    }
    close(IN);
}

# Convertit le fichier en entree (infile) du format fasta au format phylip
sub fasta2phylip
{
    my ($infile)=@_;
    
    my $nom;
    my @tab=();
    my %h1=();
    my $nbseq=0;
    
    open IN,"$infile";
    #boucle de lecture du fichier d'entr�e
    while(<IN>)
    {
	#supprime le saut de ligne afin d'aligner le nom et la s�quence
	chomp;
	#recherche ">" comme le d�but du nom d'une s�quence
	if(/^>/)
	{
	    #enl�ve le signe ">" dans le nom
	    s/>//;
	    $nom=$_;
	    #le scalaire $nom2 devient la cl� du hash
	    $h1{$nom}="";
	    #cr�e un tableau correspondant � la liste des noms
	    push(@tab,$nom);
	}
	else
	{
	    s/\s//g;
	    #concat�ne plusieurs lignes de s�quence si besoin est
	    $h1{$nom}.=$_;
	}
    }
    close IN;
    $nbseq=@tab;

    my $outfile=$infile;
    #$outfile=~s/\./_/g;
    #$outfile=~s/-/_/g;
    #$outfile.='.phy';
    
    open OUT ,"> $outfile";	
    my $longueur=length($h1{$tab[0]});
    #indique au d�but du fichier le nombre de s�quences align�es et la longueur de ces s�quences
    print OUT "\t$nbseq\t$longueur\n";
    
    #initialisation d'un curseur de position sur la s�quence
    my $pos = 0;
    #execution d'une boucle tant que le curseur n'a pas atteint le dernier �l�ment de la s�quence
    while($pos<=$longueur)
    {
	#parcourt les cl�s du hash
	foreach(@tab)
	{
	    if ($pos < 60)
	    {
		#affiche l'un en face de l'autre le nom de la s�quence et la partie de s�quence 
		#extraite et retourne � la ligne pour afficher les s�quences l'une en dessous de l'autre
		print OUT substr($_,0,10)."\t\t".substr($h1{$_},$pos,60)."\n";
	    }
	    else
	    {
		#affiche l'un en face de l'autre le nom de la s�quence et la partie de s�quence
		#extraite et retourne � la ligne pour afficher les s�quences l'une en dessous de l'autre
		print OUT "\t\t\t".substr($h1{$_},$pos,60)."\n";
	    }
	}
	#ajoute 60 au curseur pour que le foreach puisse afficher le reste de la s�quence
	$pos+=60;
	#saute une ligne entre chaque groupe de s�quence
	print OUT "\n";
    }
    close OUT;
    return($outfile);
}

# Cette fonction ouvre le fichier d'arbre ($in), et remet les v�ritables pids des prot�ines
# Donc de taille > 10 caract�res
sub retablirPID
{
    my ($in,$pidTemp) = @_;
    my $total="";
    open(IN,$in);
    while(<IN>){
	s/([^_]*_)([0-9]*)/$1$pidTemp->{$2}/g; #pour genome non modifie, abrev = 5 lettres
	#s/([a-zA-Z]{5}2_)([0-9]*)/$1$pidTemp->{$2}/g; #pour genome modifie, abrev = 5 lettres + "2"
	$total.=$_;;
    }
    close(IN);

    open(OUT,">".$in);
    print OUT $total;
    close(OUT);
}
