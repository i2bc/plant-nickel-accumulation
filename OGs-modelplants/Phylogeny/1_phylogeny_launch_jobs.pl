#!/usr/bin/perl -w
use strict;
#use threads;
#use lib "/mnt/data/team/cecile/meth_init_multithreads/Phylogeny_multithreads_v2/lib/";
use Parallel::ForkManager;

# chaque job permet de parser un fichiers blast generes pour l'ensemble des methodes.
# la procedure extrait les resultats repondant aux criteres fixes:
# Cutoff values: sandrine avait donne 0.001 pour la evalue  et 0 pour le coverage [ ou 0.001 et 70?]
# Sortie: fichiers *.tri dans le repertoire trifiles contenant le blast trié et reformatés pour input dans Phylogeny 


#############
# Declaration
#############

# Repertoires
if ($#ARGV != 3) 
{ 
	print "1_phylogeny_lauch_jobs.pl <repertoire de travail> <genome directory> <blast directory> <cpu number>\n";
	exit();
}
my $work_dir =  $ARGV[0];
my $phylodir = "$work_dir/Phylogeny/";
chdir ($phylodir);

my $dirGenome= $ARGV[1];
my $dirBlast= $ARGV[2];
my $cpu = $ARGV[3];
my $dirOut = "STEP1";
unless ( -e $dirOut) { mkdir($dirOut);}

my $evalue = 0.001; #seuil sur la evalue pour former familles d'homologues
my $palign = 0; #seuil sur le pourcentage d'alignement pour former familles d'homologues

my $pm=Parallel::ForkManager->new($cpu);

# affichage seuil
print "Resultat enregistre dans :".$dirOut;
print "\nEvalue seuil : ".$evalue;
print "\npourcentage alignement seuil : ".$palign."\n";

my @AllGenome;
my %lance=();
# Recup du nom des genomes a analyser
if (-e $dirGenome)
{
	chdir($dirGenome);     # Devient le repertoire de travail
	my @File = glob("*.fa"); # Etabli liste du nom des genomes
	foreach my $sp (@File)
	{
		  
		  $sp =~ s/\.fa$//;
		  push(@AllGenome,$sp);
	}
}

# execution des tris des blast
chdir($phylodir);
for (my $i = 0; $i <= $#AllGenome; $i ++) 
{
	my $esp1 = $AllGenome[$i];
	# Comparaison de $esp1 avec tous les autres genomes
	for (my $j = $i; $j <= $#AllGenome; $j ++) 
	{
	  
		my $esp2 =  $AllGenome[$j];
		unless(exists($lance{"$esp2\_$esp1.tri"})) # 
		{
			
				$lance{"$esp2\_$esp1.tri"} = '';
				my $pid = $pm->start and next;
				#Lancement du tri des blast dans le fichier query = esp1 et database = esp2 : extraction des resultats repondant aux criters entre ces deux genomes
				system ("perl 1_phylogeny_job.pl $esp1 $esp2 $dirBlast $evalue $palign");
				$pm->finish;
				print "$esp2\_$esp1.tri lance\n";
			
		}
	   
	}
}
$pm->wait_all_children;

print "FIN\n";

