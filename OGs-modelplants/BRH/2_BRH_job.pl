#!/usr/bin/perl -w

# Usage: Extraction des paires BRH entre chaque paire de génomes selon les seuils définis
# jobs elementaire de my2_BRH.pl system("perl my2_BRH.pl $work_dir $genomes $cpu ");
# system("perl $dirWork/BRH/my2_BRH_job.pl $dirWork $dirGenome $AllGeno[$geno1] $AllGeno[$geno2]"

my $dirWork = $ARGV[0];
my $genomeA = $ARGV[1];
my $genomeB = $ARGV[2];

my $input_dir = "$dirWork/BRH/STEP1";
my $output_dir = "$dirWork/BRH/STEP2";

# Seuil a respecter par les BHs
my $cutratio = 0.2; #cutoff sur le ratio
my $cutperc = 60; #cutoff sur le % d'alignement

my %coverage;
my %score;  
# Recuperation des scores de reference pour le calcul du ratio des scores
open (BRH, ">".$output_dir.'/'.$genomeA.'_vs_'.$genomeB.".brh");
my %linkAB;	


open (AA, $input_dir.'/'.$genomeA.'_vs_'.$genomeA.".blast.besthit") or die($genomeA.'_vs_'.$genomeA.'.blast.besthit not found');

my %maxscoreA;
while (my $line= (<AA>))
{
	chomp($line);
	my @line = split ("\t",$line);
	$maxscoreA{$line[0]} = $line[4];
}
close AA;
open (AB, $input_dir.'/'.$genomeA.'_vs_'.$genomeB.".blast.besthit") or die($genomeA.'_vs_'.$genomeB.'.blast.besthit not found');
while (my $line = (<AB>)) 
{
	chomp $line;
    my @line = split ("\t",$line); 
    #[0] id query, [1] id BH, [2] longueur alignée query, [3] longueur alignée Hit, [4] score
	my $coverage = $line[2]/$line[5]*100;
	my $score_ratio;
	if( defined($maxscoreA{$line[0]}))
	{ 
		$score_ratio = $line[4]/$maxscoreA{$line[0]};
	}
	else { $score_ratio = 0;} 

	if ($score_ratio> 0.2 and $coverage> 60)
	{                                                              
		$linkAB{$line[0]} = $line[1]; 
		$coverage{$line[0].'_vs_'.$line[1]} = $coverage;
		$score{$line[0].'_vs_'.$line[1]} = $score_ratio;
    }
}
close (AB);

my %linkBA;          
open (BA, $input_dir.'/'.$genomeB.'_vs_'.$genomeA.".blast.besthit") or die($genomeB.'_vs_'.$genomeA.'.blast.besthit not found');
open (BB, $input_dir.'/'.$genomeB.'_vs_'.$genomeB.".blast.besthit") or die($genomeB.'_vs_'.$genomeB.'.blast.besthit not found');

my %maxscoreB;
while (my $line= (<BB>))
{
	chomp($line);
	my @line = split ("\t",$line);
	$maxscoreB{$line[0]} = $line[4];
}
close BB;

while (my $line = (<BA>)) 
{
	chomp $line;
    my @line = split ("\t",$line);
	$coverage = $line[2]/$line[5]*100; 
	if(defined($maxscoreB{$line[0]})) 
	{
		$score_ratio = $line[4]/$maxscoreB{$line[0]};
		if ($score_ratio> 0.2 and $coverage> 60)
		{                                                              
			$linkBA{$line[0]} = $line[1];
			$coverage{$line[0].'_vs_'.$line[1]} = $coverage;
			$score{$line[0].'_vs_'.$line[1]} = $score_ratio;
		}	
	}
	#else # test bug division par 0
	#{
		#print $line."\n";
		#print "$. and $line[0]\n";
	#}
}
close (BA);
# write the output in the.brh file
foreach my $a (keys %linkAB)
{
	if ($linkBA{$linkAB{$a}} and $linkBA{$linkAB{$a}} eq $a)
	{
		my $b = $linkAB{$a};
		my $coverage = sprintf("%.0f",min($coverage{ $a.'_vs_'.$b},$coverage{ $b.'_vs_'.$a}));
		my $score = sprintf("%.2f",min($score{ $a.'_vs_'.$b},$score{ $b.'_vs_'.$a}));
		print  BRH $a."\t".$b."\t". $linkBA{$linkAB{$a}}."\t".$coverage."\t".$score."\n";
	}	
}

#recherche valeur minimale entre deux valeurs
sub min 
{ 
	my ($val1, $val2) = @_;
	if ($val1 < $val2) 
	{
		return $val1;
	}
	else 
	{
		return  $val2;
	}
}
