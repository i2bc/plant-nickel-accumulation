#!/usr/bin/perl
use strict;
use warnings;
use metacyc2fasta;
# Extraction des sequences des proteines de MetaCyc et ecriture d'un fichier au format fasta en sortie
# Certains liens de MetaCyc vers des sequences de Uniprot ne fonctionnent pas (deprecated IDs), exceptions gerees par un test prealable

my $dirWork = $ARGV[0]; # Repertoire de travail (git clone)
my $tag = '_'.$ARGV[1]; # tag du run (facultatif)
if ($ARGV[1]) { my $tag = '_'.$ARGV[1];}
my $file = $dirWork.'/annot'.$tag.'/parsing_proteins'; 
open IN,$file or die ("pb avec l'ouverture du fichier $file\n");
open OUT, ">annot_files/metacyc.fa";	
while (my $l = <IN>) 
{
	chomp($l);				
	my @l = split (/\t/, $l);
	my $protein = $l[0];	print "metacyc_protein=$protein\n";		
	my $metacyc_function = $l[1];	
	$metacyc_function =~ s/<.*>//g; # elimination des balises?
	print "metacyc_function=$metacyc_function\n";							
	#get data from EBI using MetaCyc IDs (XML format)
	if ( $l[3] =~ /UNIPROT=/) 
	{ 
		my $ID= $l[3];
		$ID =~ s/UNIPROT=//;
		print "Uniprot $ID\n";	
		my ($dataset, $ref, $org, $fullname, $ec_list, $sequence) = metacyc2fasta::uniprot2fasta ($ID);
		print OUT ">$ID|$ref metacyc:$metacyc_function;$dataset:$fullname;$ec_list\n$sequence\n"; 
	}
	#get data from NCBI using MetaCyc IDs (XML format)
	if ($l[3] =~ /REFSEQ=/ or $l[3] =~ /PID=/)
	{
		my $ID= $l[3];
		my $dataset;
		if ($ID =~ /REFSEQ=/ ) { $dataset = "RefSeq";} else {$dataset = "GenBank";}
		$ID =~ s/REFSEQ=|PID=//;
		print "NCBI $ID\n";	
		my ($ref, $org, $ncbi_function, $product, $ec_list, $sequence) = metacyc2fasta::ncbi2fasta ($ID);
		print OUT ">$ID|$ref metacyc:$metacyc_function;$dataset:$ncbi_function;$dataset product:$product;$ec_list\n$sequence\n"; 
	}
}
close IN;
