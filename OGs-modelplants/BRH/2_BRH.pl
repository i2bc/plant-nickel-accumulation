#!/usr/bin/perl -w
use strict;
use Parallel::ForkManager;
#use threads;

# Extraction des paires BRH entre chaque paire de génomes
# Les paires doivent repondre aux seuils fixés
# perl 2_BRH.pl [repertoire de travail] [repertoire contenant les genomes]
if($#ARGV <2)
{
	print "usage: 2_BRH.pl <repertoire de travail> <repertoire contenant les genomes> <nombre de process> ";
	exit();
}
my $dirWork = $ARGV[0]; # INPUT : Répertoire de travail
my $dirGenome = $ARGV[1]; # INPUT : Repertoire contenant les genomes
#my $threadsAllowed=$ARGV[3];
my $pm=Parallel::ForkManager->new($ARGV[2]);

# Seuil a respecter pour les BH
my $cutRatio = 0.2; #cutoff sur le ratio
my $cutPerc = 60; #cutoff sur le % d'alignement

# Creation des repertoires OUTPUT
my $dirOutBRH = $dirWork."/BRH/STEP2";
mkdir($dirOutBRH) unless (-e $dirOutBRH);

# Liste des genomes
chdir($dirGenome);
my @TmpGeno = glob("*.fa");
my @AllGeno; 
foreach my $geno (@TmpGeno) 
{
	$geno =~ s/\.fa$//;
	print "$geno\n";
	push(@AllGeno, $geno);
}
# Tri des genomes par ordre alphabetique
@AllGeno = sort(@AllGeno);

#Pour chaque paire de génomes, extraction des BRHs
for (my $geno1 = 0; $geno1 < $#AllGeno; $geno1 ++) 
{
	for (my $geno2 = $geno1+1; $geno2 <= $#AllGeno; $geno2++) 
	{
		print "BRH between $AllGeno[$geno1] and $AllGeno[$geno2]\n";
		my $pid = $pm->start and next;	
		system("perl $dirWork/BRH/my2_BRH_job.pl $dirWork $dirGenome $AllGeno[$geno1] $AllGeno[$geno2]");
		$pm->finish;
  	}
}
$pm->wait_all_children;

