#!/usr/bin/perl -w

use strict;
use threads;
use multiblastp;

# Repertoire contenant les fichiers in (genomes au format fasta) et out (blasp au format table)
my $in_dir = $ARGV[0]; # input (where are proteomes)
my $out_dir = $ARGV[1]; # output (to store blast results)
my $nb_process = $ARGV[2];
# Verification que le repertoire contenant les genomes existe
if ($#ARGV <2) { print ("allblast_with_threads.pl <fasta directory> <blast directory> <nb of threads>\n"); exit();}
die ("Error : $in_dir) not found.\n") unless (-e $in_dir);

# Creation des repertoires resultats
mkdir($out_dir) or die("writing permission denied") unless (-e $out_dir);
# lecture des fichiers fasta
my @fastafile = glob($in_dir."/*.fa");
print scalar(@fastafile) * scalar(@fastafile);

my $nb_process = $ARGV[2];

my @running = ();
my @threads;
my $i = 0;
# creation d'une database pour blastp pour chaque fichier fasta
while (scalar @threads < scalar(@fastafile))
{
	@running = threads->list(threads::running);
	if (scalar @running < $nb_process) 
	{
 		
		my @arg = ($fastafile[$i]);
		my $job1 = threads->create(\&multiblastp::run_formatdb, @arg);
		push (@threads, $job1);
		my $tid = $job1->tid;
		print "  - starting thread $tid :\t";
		print "format ".$fastafile[$i]."\n";
		$i++;
	}
	@running = threads->list(threads::running);
	foreach my $thr (@threads) 
	{
		if ($thr->is_joinable()) 
		{
			my $tid = $thr->tid;
			$thr->join;
			print "  - Results for thread $tid:\t";
			print "  - Thread $tid has been joined\n";
		}
	}

	@running = threads->list(threads::running);
}
# formatage de tous les genomes avant le lancement des blasts
print "\nJOINING pending threads\n";
while (scalar @running != 0) 
{
	foreach my $thr (@threads) 
	{
		$thr->join if ($thr->is_joinable());
	}
	@running = threads->list(threads::running);
}

@threads = ();	
# blast entre les genomes
$i = 0;
my $j = 0;
while (scalar @threads < scalar(@fastafile)*scalar(@fastafile)) #controle du nombre de threads
{
 	@running = threads->list(threads::running);
	if (scalar @running < $nb_process) 
	{
 		if(-e $fastafile[$j].".phr" and -e $fastafile[$j].".pin" and -e $fastafile[$j].".psq")
		{
			my @arg = ($fastafile[$i], $fastafile[$j], $out_dir);
			my $job2 = threads->create (\&multiblastp::run_blastp,@arg); 
			push (@threads, $job2);
			my $tid = $job2->tid;
			print "  - starting thread $tid :\t";
			print "starting task ".$fastafile[$i]."_versus_".$fastafile[$j]."\n";
			if ($j == $#fastafile) { $i++; $j=0;}
			else
			{
				$j++;
			}
		}
	}	
	@running = threads->list(threads::running);
	foreach my $thr (@threads) 
	{
		if ($thr->is_joinable()) 
		{
			my $tid = $thr->tid;
			$thr->join;
			print "  - Results for thread $tid :\t";
			print "  - Thread $tid has been joined\n";
		}
	}

	@running = threads->list(threads::running);
}

print "\nJOINING pending threads\n";
while (scalar @running != 0) 
{
	foreach my $thr (@threads) 
	{
		$thr->join if ($thr->is_joinable());
	}
	@running = threads->list(threads::running);
}

print "NB started threads = ".(scalar @threads)."\n";
print "End of main program\n";


