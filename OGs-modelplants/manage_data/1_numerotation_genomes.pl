#!/usr/bin/perl -w
# This script build the numbered genomes and proteins to input in FUNGIpath. 
# it prepares the input directory for OrthoFinder (no blast option) writing SpeciesIDs.txt, SequenceIDs.txt files and for each genome the Species<N>.fa symbolic link to the numbered genome file.
use File::Basename;

if ($#ARGV <1) 
{	 
	 print "1_numerotation_genomes.pl <git directory> <fasta files directory>\n";	 
	 exit();
} 
my $work_dir = $ARGV[0];
my $fasta_dir = $ARGV[1];
chdir($work_dir);

mkdir ("numfasta");
mkdir ("numfiles");
mkdir ("OrthoFinder/input");
my $i = -1; # rank of the genome (first rank is 0)
# open the genome list (input in OrthoFinder)  
open NF,">OrthoFinder/input/SpeciesIDs.txt";

foreach my $f ( glob("$fasta_dir/*.fa") ) 
{
	$i++;
	my @suffix = (".fa");
	my $genome = basename($f, @suffix);
	#print $genome."\n";
	open IN, $f;
	#lecture du fichier fasta du genome 	
	open NP, ">numfiles/".$i.".name"; # numbered genomes
	
	my $out;
	$j = -1; # rank of the protein (first rank is 0)
	while (my $l=<IN>) 
	{ 
		if($l =~ /^>/) 
		{
			$j++;
			@l = split(/\s+/, $l);
			$l[0] =~ s/>//;
			$out .=  '>'.$i."_".$j."\n";
			print NP $i."_".$j.": ".$l[0]."\n";	# write the .name file (pull files for Orthofinder)		
		} 
		else 
		{ 
			$out .=  $l;
		}
	} 
	print NF "$i: $genome\.fa\n";
	# write the numbered fasta file	
	open OUT, '>numfasta/'.$i.'.fa'; 
	print OUT $out;
	close OUT;
	close NP;
	# symbolic links for OrthoFinder
	system "ln -s $work_dir/numfasta/$i.fa $work_dir/OrthoFinder/input/Species$i.fa";
	# SequenceIDs file for OrthoFinder
	system "cat numfiles/* >  OrthoFinder/input/SequenceIDs.txt" ; # (pull files for Orthofinder)	
}
close NF;
