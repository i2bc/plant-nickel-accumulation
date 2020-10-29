#!/usr/bin/perl -w

# programme qui va executer autant de jos qu'il n'y a de familles proteiques construites
# chaque job va generer l'arbre enracine de la famille

#perl 4_phylogeny_launch_jobs_clusterLIX.pl /home/cpereira/

#############
# Declaration
#############

# Parametres
#my $dirWork = '/projet/extern/sgrossetete/'; #INPUT : Repertoire de travail
my $dirGENERAL=$ARGV[0];
my $dirWork = $dirGENERAL.'Phylogeny/';
my $dirData = $dirGENERAL.'datatest7esp/'; #INPUT : Repertoire contenant genomes au format FASTA
my $dirExec = $dirWork.'exec/'; #INPUT : Repertoire contenant les programmes a utiliser (phyml, muscle...)
my $dirTmp = $dirWork.'tmp/'; #INPUT : Repertoire temporaire
my $dirOut = $dirWork.'result/'; # INPUT : repertoire contenant les fichiers generes par la methode phylogenie
my $dirGroup = $dirOut.'groupDist_1e3_70test7esp/'; #INPUT : Repertoire contenant famille proteique
my $dirTree = $dirOut.'arbre_1e3_70test7esp/'; #OUTPUT : Repertoire contenant les arbres enracines
my $dirAnalyse = $dirOut.'analyse_arbretest7esp/'; #OUTPUT : Repertoire contenant les analyses des arbres

# Autres variables
my @File = ();
my $nbJobMax = 1000;
my $nb=0;
my $n=0;
my $job=0;#numero du jobs dans lequel va etre place l'une des requettes

###########
# Programme
###########
mkdir($dirTmp);
mkdir($dirTree);
mkdir($dirAnalyse);
chdir($dirGroup);
@File = glob("group*");
$nb=@File;#$nb contient le nombre de jobs à executer au max
$nb=int($nb/$nbJobMax)+1;#$nb contient le nombre de jobs a faire par script
$n=0;

chdir($dirWork);
# pour chaque groupe, on va construire l'alignement et l'arbre
# on le fera via des jobs envoyés sur le cluster.
# attention, on ne doit pas avoir plus de 1000 jobs envoyés
open(OUT,">".$dirTmp."/job$job.qsub");
print OUT "#!/bin/bash\n";
print OUT '#$ -N '."job$job\n";
print OUT '#$ -cwd'."\n";
print OUT '#$ -o '.$dirTmp."/job$job.out\n";
print OUT '#$ -e '.$dirTmp."/job$job.err\n\n";
foreach my $file (@File) {
    my $tmpFile = $file;
    $tmpFile =~ s/groupe/arbre/;
    #print 'perl '.$dirWork."4_phylogeny_job1.pl $dirExec $dirTmp $dirTree $dirGroup"."$file ".$dirData."\n";
    print OUT 'perl '.$dirWork."4_phylogeny_job1_clusterLIX.pl $dirExec $dirTmp $dirTree $dirGroup"."$file ".$dirData."\n";
    print OUT 'perl '.$dirWork."4_phylogeny_job2_clusterLIX.pl $dirTree $dirGroup $tmpFile".".txt > ".$dirAnalyse.$tmpFile."\n";
    $n++;
    if($n>$nb){
    	close(OUT);
    	system("qsub $dirTmp/job$job.qsub");
    	unlink("$dirTmp/job$job.qsub");
    	$job++;
    	open(OUT,">".$dirTmp."/job$job.qsub");
			print OUT "#!/bin/bash\n";
			print OUT '#$ -N '."job$job\n";
			print OUT '#$ -cwd'."\n";
			print OUT '#$ -o '.$dirTmp."/job$job.out\n";
			print OUT '#$ -e '.$dirTmp."/job$job.err\n\n";
			$n=0;
    }
}
#pour le dernier job
system("qsub $dirTmp/job$job.qsub");
unlink("$dirTmp/job$job.qsub");
