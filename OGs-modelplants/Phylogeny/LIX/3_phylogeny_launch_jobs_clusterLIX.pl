#!/usr/bin/perl -w

# Script qui lance un job par famille proteique dont la taille depasse celle maximale
# Le job permet d'executer MCL afin de reduire la taille du groupe
# perl 3_phylogeny_launch_jobs.pl /home/cecile/Bureau/miseajourBD/local/Phylogeny/ 35
use strict;

#############
# Declaration
#############
print $#ARGV."\n";
die ("Erreur : le dossier de depart et la taille maximale d'un groupe doit etre donne en argument.\n") unless ($#ARGV==1);
my $sizeMax = $ARGV[1]; # INPUT : nombre maximal de proteines dans une famille ==> si supérieur a ce seuil, on execute MCL pour dimunuer la taille de la famille
#5*nombre de genomes: 35
# Parametres
#my $dirWork = '/projet/extern/sgrossetete/'; # INPUT : repertoire de travail
my $dirWork=$ARGV[0];

my $dirGroup = $dirWork.'result/groupDist_1e3_70test7esp/'; # INPUT : repertoire ou se situent les familles proteiques
my $dirTmp = $dirWork.'tmp/'; # OUTPUT : repertoire temporaire
my $dirIn = $dirTmp.'mcl_input/'; # OUTPUT : repertoire ou seront stockés les fichiers d'entrée mcl générés à partir des fichiers groupes
my $dirOut = $dirTmp.'mcl_output/'; # OUTPUT : repertoire de sortie mcl
my $nbJobMax = 1000; # INPUT : nombre maximal de jobs a executer simultanement (limite fixée sur le cluster)

# Autres variables
my @Line = ();
my @File = ();
my $file = '';
my $nbProt = 0;

#############
# Programme
#############

#creation du repertoire d'entree et de sortie si besoin
system("mkdir ".$dirTmp) unless (-e $dirTmp);
system("mkdir ".$dirIn) unless (-e $dirIn);
system("mkdir ".$dirOut) unless (-e $dirOut);

#recup de la liste de tous les groupes
chdir($dirGroup);
@File = glob("groupe*");
print @File." groupes\n";

#pour chaque groupe
#attention ici, je n'ai pas le droit à plus de 1000 jobs en meme temps => voir combiens je dois mettre de groupes par jobs
my $nb=@File;#$nb contient le nombre de jobs à executer au max
$nb=int($nb/$nbJobMax)+1;#$nb contient le nombre de jobs a faire par script
my $n=$nb+1;;
my $job=0;#numero du jobs dans lequel va etre place l'une des requettes

foreach $file (@File) {
	$n++;
	if($n>$nb){
		$n=0;
		unless($job==0){
			close(OUT);
			system('qsub '.$dirTmp.'phylo3_'.$job.'.sh') ;
			unlink($dirTmp.'phylo3_'.$job.'.sh');
		}
		$job++;
		open (OUT,'>'.$dirTmp.'phylo3_'.$job.'.sh');    
		print OUT '#$ -S /bin/bash'."\n";
   	 	print OUT '#$ -V'."\n";
    		print OUT '#$ -m a'."\n";
    		print OUT '#$ -o '.$dirTmp.'3phylo'.$job.'.out'."\n";
    		print OUT '#$ -e '.$dirTmp.'3phylo'.$job.'.err'."\n";
    		print OUT '#$ -cwd'."\n";
	}
	open GROUP,$dirGroup.$file; 
	$_ = <GROUP>; #on lit uniquement la 1ere ligne qui contient les id du groupe (nodes id )
	@Line = split (/ /,$_);
	close (GROUP);
	pop(@Line); # suppresion de '(nodes'
	shift(@Line); # supprresion de ' )'
	$nbProt = $#Line + 1; #nb de prot dans le groupe
	#si nombre de prot superieur au seuil, on realise MCL
	if ($nbProt > $sizeMax) {
		print OUT "echo 'Analyse mcl de ".$dirGroup.$file." ($nbProt proteins)'\n";
		system("cp ".$dirGroup.$file." ".$dirIn.$file);
		print OUT 'perl '.$dirWork."3_phylogeny_job_clusterLIX.pl $file $dirIn $dirOut $dirGroup $sizeMax 3 $dirWork $dirTmp\n";   
	}
}
close(OUT);
system('qsub '.$dirTmp.'phylo3_'.$job.'.sh') ;
unlink($dirTmp.'3phylo'.$job.'.sh');

