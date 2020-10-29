#!/usr/bin/perl -w
# Code by Christine Drevet (2016)
# Replace queries and hits names by numbers (<genome number>_<protein number>) in blast output (table formatted) for input in BRH, Inparanoid Phylogeny and  OrthoFinder. The blasts ouput are generated at the run_blast step.
use strict;
use File::Basename;
my $work_dir = $ARGV[0];
my $blast_dir = $ARGV[1];

if ($#ARGV<1)
{
	print "numerisation of queries and hits names of blastp files\n";
	print "perl numerotation_blast.pl <work directory> <dossier des Blasts directory>\n"; 
	exit();
}
open NG,"$work_dir/OrthoFinder/input/SpeciesIDs.txt"; # see orthofinfinder documentation for details
my %num_genome;
while (my $lg=<NG>)
{
	chomp($lg);
	$lg =~ s/\.fa//;
	my @lg = split(/: /, $lg);
	$num_genome{$lg[1]} = $lg[0];
}
close NG;

mkdir $work_dir."/numblast" unless (-e work_dir."/numblast");
foreach my $f ( glob("$blast_dir/*.blast") ) 
{
	my @suffix = (".blast");
	my $blastfile = basename($f, @suffix);
	my @blastfile = split(/_vs_/,$blastfile);
	my $query = $blastfile[0]; print $query."\t";
	my $num_query = $num_genome{$query}; print $num_query."\n";
	my $hit = $blastfile[1];print $hit."\t";
	my $num_hit = $num_genome{$hit};	print $num_hit."\n";
	
	open NQ, $work_dir."/numfiles/$num_query.name";	
	my %nb_query;
	while (my $lq=<NQ>)
	{
		chomp($lq);
		my @lq = split(/: /, $lq);
		$nb_query{$lq[1]} = $lq[0];
	}
	close NQ;
	
	open NH, $work_dir."/numfiles/$num_hit.name";
	my %nb_hit;
	while (my $lh=<NH>)
	{
		chomp($lh);
		my @lh = split(/: /, $lh);
		$nb_hit{$lh[1]} = $lh[0];
		print $lh[0].'->'.$lh[1].'->'.$nb_hit{$lh[1]} ."\n" ;
	}
	close NH;
	
	open BLAST , $f;
	
	open NUMBLAST, ">$work_dir/numblast/$num_query\_vs_$num_hit.blast";

	while (my $l=<BLAST>) 
	{ 
		unless( $l =~ /^#/)
		{
			print $l;
			chomp $l;
			my @l = split (/\t/, $l);
			print  NUMBLAST $nb_query{$l[0]}."\t".$nb_hit{$l[1]} ;		
			for (my $i=2 ; $i <= $#l ; $i++ ) { print NUMBLAST "\t".$l[$i];}
			print  NUMBLAST "\n";
		}
	}
	system "ln -s $work_dir/numblast/$num_query\_vs_$num_hit.blast $work_dir/OrthoFinder/input/Blast$num_query\_$num_hit.txt";
}
