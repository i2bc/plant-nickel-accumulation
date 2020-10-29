#!/usr/bin/perl -w
use strict;
#use DBI;

# programme qui permet de parser les fichiers blast generes pour reecrire les resultats dans le format attendu par les programmes d'analyse d'arbres ecrits par Frederic Lemoine
# parser pour phylogeny
#modifié le 19mai pour ne plus avoir d'erreur liee a une taille de seq abs
#############
# Declaration
#############

# recuperation des arguments
die ("Erreur : six arguments doivent etre donnes: le dossier de sortie, l'evalue seuil, le pourcentage d'alignement seuil, le nom de la premiere espece, le nom de la seconde espece, le dossier qui contient les blasts\n") unless ($#ARGV == 5); 
my $dirOut = $ARGV[0]; # OUTPUT : Repertoire ou seront enregistres les fichiers generes par le programme
my $evalue = $ARGV[1]; # INPUT : seuil sur la evalue pour former familles d'homologues
my $palign = $ARGV[2]; # INPUT : seuil sur le pourcentage d'alignement pour former familles d'homologues
my $esp1 = $ARGV[3] ; # INPUT : nom de la premiere espece
my $esp2 = $ARGV[4] ; # INPUT : nom de la seconde espece
my $dirBlast = $ARGV[5] ; # INPUT : repertoire contenant les fichiers BLAST de la comparaison des 2 genomes
#my $dirBh =  $ARGV[6]; # INPUT : repertoire contenant les fichiers BH
#my $dossierTMP = $ARGV[7];

my @Genome = ($esp1,$esp2);
my @Line = ();
my $File = "";
my %LongSeq = (); # a chaque id correspond la longueur de la sequence
my %ScoreMax = (); #a chaque id correspond le score Blast avec sa sequence,
my %Id = ();
my $id = '';
my $longAlign1 = 0;
my $longAlign2 = 0;
my $dist = 0;
my $id1 = '';
my $id2 = '';

###########
# Programme
###########

#Pour chacun des 2 génomes, lecture du score de reference (pour calcul du ratio) et de la longueur de la sequence
chdir($dirBlast);
foreach my $geno (@Genome) {
	$File = $geno."_".$geno."\.blast"; 
	open(IN,$dirBlast."/".$File) or die ("impossible d'ouvrir le fichier blast $File\n");
	while (<IN>) {
		unless(/^#/){
			chomp($_);
			@Line = split(/\t/,$_);
			if($Line[0] eq $Line[1]){
				$ScoreMax{$Line[0]} = $Line[11];
				$LongSeq{$Line[0]} = $Line[12];
			}
		}	
	}
	close(IN);
}

#ouverture en ecriture du fichier resultat
open OUT,">$dirOut/".$esp1.'_'.$esp2.".tri";

$File = "$esp1\_$esp2.blast"; 
open IN,$dirBlast.'/'.$File;
while(<IN>) {
	unless(/^#/){
		chomp ($_);
		($id1,$id2,$longAlign1,$longAlign2,$dist) = extractRes($_,$esp1,$esp2,$evalue,$palign,\%Id,\%ScoreMax,\%LongSeq);
		if (($longAlign1 != 0) && ($longAlign2 != 0)) {
			print OUT $esp1."_".$id1."\t".$esp2."_".$id2."\t \t \t \t \t".$LongSeq{$id1}."\t".$LongSeq{$id2}."\t \t \t \t$longAlign1\t$longAlign2\t".$dist."\n";
		}
	}
}
close (IN);

#Si les 2 especes ne sont pas les memes et que fichier .tri n'existe pas, on recupere les resultats inverses
if ($esp1 ne $esp2) {
    if (-e "$dirOut/".$esp2."_".$esp1.".tri") {
	print "File $dirOut/".$esp2."_".$esp1.".tri already exist. No new calcul.\n";
    }
    else {
        $File ="$esp2\_$esp1.blast";  
	open IN,$dirBlast.'/'.$File;
	while(<IN>) {
		unless(/^#/){
        		chomp ($_);
			($id1,$id2,$longAlign1,$longAlign2,$dist) = extractRes($_,$esp2,$esp1,$evalue,$palign,\%Id,\%ScoreMax,\%LongSeq);
			 if (($longAlign1 != 0) && ($longAlign2 != 0)) {
				print OUT $esp1."_".$id2."\t".$esp2."_".$id1."\t \t \t \t \t".$LongSeq{$id2}."\t".$LongSeq{$id1}."\t \t \t \t$longAlign2\t$longAlign1\t".$dist."\n";
			}
		}
	}
	close(IN);
    }
}
close(OUT);


###########
# Fonction
###########


sub extractRes {
    my ($line,$esp1,$esp2,$evalue,$palign,$refId,$refScoreMax,$refLongSeq) = @_;
    my @Line = split(/\t/,$line);
 
    my $id1 = $Line[0];
    my $id2 = $Line[1];
    my $percAlign1 = 0;
    my $percAlign2 = 0;
    my $dist = 0;
    #a moins que couple deja enregistre, si 2 hits pour meme prot, on ne prend que le 1er
    #2eme condition dans le cas ou esp1 = esp2
    unless (exists($$refId{$id1.$id2}) || exists($$refId{$id2.$id1})) {
	if ( ($Line[10] < $evalue) && ($id1 ne $id2) ) {
        	if (exists($$refScoreMax{$id1}) && ($$refScoreMax{$id1} != 0)) {
			#Calcul du ratio
			$dist =  sprintf("%.2f", $Line[11] / $$refScoreMax{$id1});
			#die "$id1\t$esp1\n$id2\t$esp2\nNo length for $id1 ($esp1)\n" unless (exists($$refLongSeq{$id1}));
			#die "$id1\t$esp1\n$id2\t$esp2\nNo length for $id2 ($esp2)\n" unless (exists($$refLongSeq{$id2}));
			unless(exists($$refLongSeq{$id1})){
				$$refLongSeq{$id1}=$Line[12];
			}
			unless(exists($$refLongSeq{$id2})){
				$$refLongSeq{$id2}=$Line[13];
			}
			$percAlign1 = ($Line[7]-$Line[6] + 1)*100/$$refLongSeq{$id1};
			$percAlign2 = ($Line[9]-$Line[8] + 1)*100/$$refLongSeq{$id2};

			#si une des 2 seq a % d'alignement superieur au seuil fixe
          		if (($percAlign1 > $palign) || ($percAlign2 > $palign)) {
           	     		$$refId{$id1.$id2}='ok' ;
                 		return($id1,$id2,($Line[7]-$Line[6] + 1),($Line[9]-$Line[8] + 1),$dist);
          		}
        	}
        	else {
            		print "Warning No max score for $id1 ($esp1 et $esp2) $id1 non considéré\n";
            		#exit();
			return ($id1,$id2,0,0,0);
	    	}
    	}
    }
    return ($id1,$id2,0,0,0);
}
