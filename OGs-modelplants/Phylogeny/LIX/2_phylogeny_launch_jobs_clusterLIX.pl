#!/usr/bin/perl -w

# programme qui permet de construire les familles protéiques necessite 2_phylogeny_job.pl

#11/05/2012

#perl 2_phylogeny_launch_jobs_clusterLIX.pl /home/cpereira/

#############
# Declaration
#############
# Repertoire
#my $dirWork = '/projet/extern/sgrossetete/'; # INPUT : repertoire de travail 
my $dirWork = $ARGV[0].'Phylogeny/';
my $dirLib = $dirWork.'lib/'; # INPUT : repertoire avec les librairies utilisees
my $dirOut = $dirWork.'result/'; # INPUT : repertoire contenant les fichiers generes par la methode phylogenie
my $dirGenomeNew = $ARGV[0].'data/'; # INPUT : Repertoire contenant les nouveaux genomes au format fasta
my $dirInComparison = $dirOut.'dataPhylo_1e3_0/'; # INPUT : repertoire avec les fichier .tri genere dans les versions precedentes
my $dirOutFamily = $dirOut.'groupDist_1e3_70/'; #OUTPUT : repertoire avec les librairies utilisees
my $dirTmp = $dirWork.'tmp/'; # OUTPUT : repertoire temporaire
my $dirGenomeOld=$ARGV[0].'old';
my @Genome = ();
my @File = ();
my $file = '';
my $stop = 0;

#############
# Programme
#############

die ("Error: No $dirLib directory\n") unless (-e $dirLib);
die ("Error: No $dirGenomeNew directory\n") unless (-e $dirGenomeNew);
die ("Error: No $dirInComparison directory\n") unless (-e $dirInComparison);
mkdir ($dirOutFamily) unless (-e $dirOutFamily);
mkdir($dirTmp)unless(-e $dirTmp);
# Verification que l'on a bien tous les fichiers
if (-e  $dirGenomeNew){
	chdir($dirGenomeNew);
	@File = glob("*fa");
	foreach $file (@File){
		$file =~ s/\.fa$//;
		push(@Genome,$file);
	}
}
@Genome = sort(@Genome);

# Test que tous les fichiers INPUT existent
for (my $i = 0; $i <= $#Genome; $i ++){
	for (my $j = $i; $j <= $#Genome; $j ++){
		unless ((-e $dirInComparison.$Genome[$i].'_'.$Genome[$j].'.tri') || (-e $dirInComparison.$Genome[$j].'_'.$Genome[$i].'.tri')) {
			$stop = 1;
			print "Error : No input files for $Genome[$i] and $Genome[$j]\n";
		}
	}
}
die ("Error : input files missing.\n") if ($stop == 1);
chdir($dirWork);
open OUT,'>'.$dirTmp.'BuildGroup.sh';
print OUT '#$ -S /bin/bash'."\n";
print OUT '#$ -V'."\n";
print OUT '#$ -m a'."\n";
print OUT '#$ -o '.$dirTmp.'BuildGroup.out'."\n";
print OUT '#$ -e '.$dirTmp.'BuildGroup.err'."\n";
print OUT '#$ -cwd'."\n";
print OUT "perl -I $dirLib 2_phylogeny_job_clusterLIX.pl $dirInComparison 70  $dirOutFamily $dirGenomeNew $dirGenomeOld\n" ;
close(OUT);

system('qsub '.$dirTmp.'BuildGroup.sh') ;
unlink($dirTmp.'BuildGroup.sh');
