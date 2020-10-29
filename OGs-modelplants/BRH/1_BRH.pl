#!/usr/bin/perl -w

#use threads;
use Parallel::ForkManager;
#load my1_BRH_job.pl for each file of the blast directory
#cecile pereira (21/04/2012)
#Usage: trier les resultats de blastp, sorties du type: ID du depart ainsi que de son meilleur hit, des longueurs alignees et du score
# usage in Lancer_BRH.pl: system("perl my1_BRH.pl $blastp $work_dir $cpu");
if ($#ARGV < 2) 
{
	print "perl 1_BRH.pl <repertoire de travail> <repertoire des blasts> <nombre de process> ";
	exit();
}
$dirSCRIPT=$ARGV[0]."/BRH/";
$blast_dir = $ARGV[1];
$dirOUT = "$dirSCRIPT/STEP1/";
mkdir($dirOUT) unless (-e $dirOUT);
opendir(BL,$blast_dir) or die("No blast files!");

#$threadsAllowed=$ARGV[3];

@fb=readdir(BL);
closedir(BL);
#my @jobs=();
$pm=Parallel::ForkManager->new($ARGV[2]);

foreach $b(@fb)
{
	unless(($b eq ".") || ($b eq ".."))
	{
	   #push(@jobs,threads->create(sub{
		my $pid = $pm->start and next;
		system("perl $dirSCRIPT/my1_BRH_job.pl $blast_dir/$b $dirOUT");
		$pm->finish;
	   #}));
	   #sleep(1) while(threads->list(threads::running)>=$threadsAllowed);
	}
}
$pm->wait_all_children;
#$_->join for @jobs;
