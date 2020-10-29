#!/usr/bin/perl -w
#use threads;
use Parallel::ForkManager;
if ($#ARGV != 2)
{
	print "The script create fasta files for each group, run MUSCLE, PhyML and retree\n";
	print "my4_phylogeny_launch_jobs.pl <repertoire de travail (git) > <repertoire des genomes> <nombre de coeur>\n";
exit();
}
	
my $dirWork = $ARGV[0];
my $dirGENERAL = $dirWork.'/Phylogeny';
my $dirData = $ARGV[1]; #INPUT : Repertoire contenant genomes au format FASTA
my $dirExec = $dirGENERAL.'/exec/'; #INPUT : Repertoire contenant les programmes a utiliser (phyml, muscle...)
my $pm=Parallel::ForkManager->new($ARGV[2]);

my $dirTmp=$dirGENERAL.'/tmp/';#repertoire contenant les fichiers fasta des groupes, alignements, arbres phyML
my $dirGroup = $dirGENERAL.'/STEP2/'; #INPUT : Repertoire contenant famille proteique
my $dirTree = $dirGENERAL.'/STEP4-1/'; #OUTPUT : Repertoire contenant les arbres enracines
my $dirAnalyse = $dirGENERAL.'/STEP4-2/'; #OUTPUT : Repertoire contenant les analyses des arbres

# Autres variables
my @File = ();
my $i = 0;
my $nbJob = 0;
my @jobs=();

###########
# Programme
###########
mkdir($dirTmp);
mkdir($dirTree);
mkdir($dirAnalyse);
chdir($dirGroup);
# lecture des fichiers a traiter (le traitement d'un fichier est une tache parallelisable)
@File = glob("group*"); 

chdir($dirGENERAL); # dossier Phylogeny du git
# pour chaque groupe, on va construire l'alignement et l'arbre
foreach my $file (@File) 
{
    my $pid=$pm->start and next;
   
    #print 'perl '.$dirGENERAL."/my4_phylogeny_job1.pl $dirExec $dirTmp $dirTree $dirGroup"."$file ".$dirData."\n";
    system('perl '.$dirGENERAL."/my4_phylogeny_job1.pl $dirExec $dirTmp $dirTree $dirGroup"."$file ".$dirData);
    my $tmpFile = $file; # fichier de sortie de l'etape 4-1
    $tmpFile =~ s/groupe/arbre/;
    system('perl '.$dirGENERAL."/my4_phylogeny_job2.pl $dirTree $dirGroup $tmpFile".".txt > ".$dirAnalyse.$tmpFile);
    $pm->finish; # création d'un fichier plat décrivant les relations entre les proteines du groupe courrant
}
$pm->wait_all_children;

