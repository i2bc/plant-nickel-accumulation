#!/usr/bin/perl -w

# Script qui lance un job par famille proteique dont la taille depasse une maximale
# Le job permet d'executer MCL afin de reduire la taille du groupe

use strict;
#use threads;
use Parallel::ForkManager;

if ( $#ARGV != 2) 
{ 
	print "lancement recursif de MCL sur les groupes d'effectif superieur au cutoff\n";
	print "3_phylogeny_launch_jobs.pl <nombre de genome> <nombre de coeur> <dossier de travail>\n" ;
	exit();
}

# definition de la taille limite des groupes (5*nombre de genomes)
my $sizeMax = 5*$ARGV[0]; 
my $pm=Parallel::ForkManager->new($ARGV[1]);

my $dirWork= $ARGV[2].'/Phylogeny/';
my $dirGroup = $dirWork.'STEP2/'; # INPUT : repertoire ou se situent les familles proteiques
my $dirTmp = $dirWork.'STEP3/'; # OUTPUT : repertoire temporaire
my $dirIn = $dirTmp.'mcl_input/'; # OUTPUT : repertoire ou seront stockés les fichiers d'entrée mcl générés à partir des fichiers groupes
my $dirOut = $dirTmp.'mcl_output/'; # OUTPUT : repertoire de sortie mcl

# Autres variables
my @Line = ();
my @File = ();
my $file = '';
my $nbProt = 0;
my @jobs=();

#############
# Programme
#############

#creation du repertoire d'entree et de sortie si besoin
system("mkdir ".$dirTmp) unless (-e $dirTmp);
system("mkdir ".$dirIn) unless (-e $dirIn);
system("mkdir ".$dirOut) unless (-e $dirOut);

#recup de la liste de tous les groupes
chdir("$dirGroup");
@File = glob("groupe*");
print @File." groupes\n"; 

#pour chaque groupe
foreach $file (@File) {
	open GROUP,$dirGroup.$file; 
	$_ = <GROUP>; #on lit uniquement la 1ere ligne qui contient les id du groupe (nodes id )
	@Line = split (/ /,$_);
	close (GROUP);
	pop(@Line); # suppresion de '(nodes'
	shift(@Line); # supprresion de ' )'
	$nbProt = $#Line + 1; #nb de prot dans le groupe
	#si nombre de prot superieur au seuil, on realise MCL
	if ($nbProt > $sizeMax) 
	{
	   #push(@jobs,threads->create(sub{		
	   my $pi=$pm->start and next;
	   print "Analyse mcl de ".$dirGroup.$file." ($nbProt proteins)\n";
		system("cp ".$dirGroup.$file." ".$dirIn.$file);
		system('perl '.$dirWork."3_phylogeny_job.pl $file $dirIn $dirOut $dirGroup $sizeMax 3 $dirWork $dirTmp") ;    
	   $pm->finish;
	   #}));
	   #sleep(1) while(threads->list(threads::running)>=$threadsAllowed);
	}
}
$pm->wait_all_children;
#$_->join for @jobs;
