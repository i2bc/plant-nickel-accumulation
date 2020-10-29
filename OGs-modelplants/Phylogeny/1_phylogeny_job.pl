#!/usr/bin/perl -w
use strict;

# programme qui permet de parser les fichiers blast generes pour reecrire les resultats dans le format attendu par les programmes d'analyse d'arbres ecrits par Frederic Lemoine
# parser pour phylogeny
# Blast++ files format: blastp -query <genome fasta file> -db <genome blast DB> -outfmt \"7 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen\" -seg yes -out  <blast file>
#############
# Declaration
#############

# recuperation des arguments
die ('perl 1_phylogeny_job.pl esp1 esp2 dirBlast evalue palign:variables non renseignées \n') unless ($#ARGV == 4); 

my $spA = $ARGV[0] ; # INPUT : nom de la premiere espece
my $spB = $ARGV[1] ; # INPUT : nom de la seconde espece
my $blast_dir = $ARGV[2] ; # INPUT : repertoire contenant les fichiers BLAST de la comparaison des 2 genomes
my $result_dir = "STEP1/";
print $result_dir."\n";
my $evalue_cutoff = $ARGV[4]; # INPUT : seuil sur la evalue pour former familles d'homologues
my $alignment_cutoff = $ARGV[5]; # INPUT : seuil sur le pourcentage d'alignement pour former familles d'homologues

#Pour chacun des 2 génomes, lecture du score de reference (pour calcul du ratio) et de la longueur de la sequence
my $blast_spA_vs_spA = $blast_dir."$spA\_vs_$spA.blast";
my %ScoreMax; # bitscore of each sequence against itself (for spA and spB) 
my %LongSeq; # length of each sequence (for spA et spB) ???? data are in field 12 and 13 of tabular blast output
my %Id;

# get maximum bit scores for proteins of specie A from the blast file of specieA genome against itself
open (IN,"$blast_spA_vs_spA"); 
print $blast_spA_vs_spA."\n";
while (my $l= (<IN>)) 
{
	unless($l =~ /^#/)
	{
		chomp($l);
		my @line = split(/\t/,$l);
		if($line[0] eq $line[1]) # search maximum bit score for each sequence ( blast against itself)
		{
			$ScoreMax{$line[0]} = $line[11]; # field 0 is query_id; field 11 is bit score 
			$LongSeq{$line[0]} = $line[12];	 # field 12 is query length
			
		}	
	}	
}
close(IN);

# get maximum bit scores for proteins of specie B from the blast file of specie B genome against itself
my $blast_spB_vs_spB = $blast_dir."$spB\_vs_$spB.blast";
print $blast_spB_vs_spB."\n";
open (IN, "$blast_spB_vs_spB"); 
while (my $l = (<IN>)) 
{
	unless($l =~ /^#/)
	{
		chomp($l);
		my @line = split(/\t/,$l);
		if($line[0] eq $line[1]) # search maximum bit score for each sequence ( blast against itself)
		{
			$ScoreMax{$line[0]} = $line[11]; # field 0 is query_id; field 11 is bit score 
			$LongSeq{$line[0]} = $line[12];	 # field 12 is query length
			
		}	
	}
}
close(IN);

#ouverture en ecriture du fichier resultat
open (OUT,">$result_dir$spA\_vs_$spB.tri");
print "$spA\_vs_$spB.tri\n";
my $blast_spA_vs_spB = "$blast_dir/$spA\_vs_$spB.blast"; 

# parse the blast file (whole genome A against whole genome B; tabular output) 
open IN, $blast_spA_vs_spB;
while( my $l = <IN>) 
{
	unless($l =~ /^#/)
	{
		chomp ($l);
		my ($query,$hit,$longAlign1,$longAlign2,$dist) = parse_blast_line($l,$spA,$spB,$evalue_cutoff,$alignment_cutoff,\%Id, \%ScoreMax,\%LongSeq);
		
		if (($longAlign1 != 0) && ($longAlign2 != 0)) # validated blast according to cutoff valueshead 
		{
			print OUT $query."\t".$hit."\t \t \t \t \t".$LongSeq{$query}."\t".$LongSeq{$hit}."\t \t \t \t$longAlign1\t$longAlign2\t".$dist."\n";
		}
	}
}
close (IN);

#Si les 2 especes ne sont pas les memes et que fichier .tri n'existe pas, on recupere les resultats inverses
if ($spA ne $spB) 
{
	if (-e $result_dir.'/'.$spB."_vs_".$spA.".tri") 
	{
		print "File ".$result_dir.'/'.$spB."_vs_".$spA.".tri already exist. No new calcul.\n";
	}
	else 
	{
		my $File="$spB\_vs_$spA.blast";
		open (IN,$blast_dir.'/'.$File);
		while(my $l = <IN>) 
		{
			unless($l =~ /^#/)
			{
				chomp ($l);
				my ($query,$hit,$longAlign1,$longAlign2,$dist) = parse_blast_line($l,$spB,$spA,$evalue_cutoff,$alignment_cutoff,\%Id,\%ScoreMax,\%LongSeq);
				if (($longAlign1 != 0) && ($longAlign2 != 0)) 
				{
			      	print OUT $hit."\t".$query."\t \t \t \t \t".$LongSeq{$hit}."\t".$LongSeq{$query}."\t \t \t \t$longAlign2\t$longAlign1\t".$dist."\n";
				}
			}
		}
		close(IN);
	}
}
close(OUT);

###########
# Fonction
###########


sub parse_blast_line 
{
    my ($line,$spA,$spB,$evalue_cutoff,$alignment_cutoff, $refId, $refScoreMax,$refLongSeq) = @_;
    my @Line = split(/\t/,$line);
 
    my $query = $Line[0]; 
    my $hit = $Line[1];
    my $percAlign1 = 0;
    my $percAlign2 = 0;
    my $dist = 0;
    
	#a moins que couple deja enregistre, si 2 hits pour meme prot, on ne prend que le 1er
    #2eme condition dans le cas ou esp1 = esp2
	unless (exists($$refId{$query.$hit}) || exists($$refId{$hit.$query})) 
	{
		if ( ($Line[10] < $evalue_cutoff) && ($query ne $hit) )  # field 10 = tabular blast evalue
		{
			if (exists($$refScoreMax{$query}) && ($$refScoreMax{$query} != 0)) 
			{
				#Calcul du ratio
				$dist =  sprintf("%.2f", $Line[11] / $$refScoreMax{$query}); # field 11 = tabular blast bit score
				die "$query\t$spA\n$hit\t$spB\nNo length for $query ($spA)\n" unless (exists($$refLongSeq{$query}));
				die "$query\t$spA\n$hit\t$spB\nNo length for $hit ($spB)\n" unless (exists($$refLongSeq{$hit}));
				my $val_align_query = ($Line[7]-$Line[6] + 1)*100/$$refLongSeq{$query}; # field 6 = tabular blast q.start; field 7 = tabular blast q. end 
				my $val_align_hit = ($Line[9]-$Line[8] + 1)*100/$$refLongSeq{$hit};  # field 8 = tabular blast s.start; field 9 = tabular blast s. end 

				#si une des 2 seq a % d'alignement superieur au seuil fixe
				if (($val_align_query > $alignment_cutoff) || ($val_align_hit > $alignment_cutoff)) 
				{
					$$refId{$query.$hit}='ok' ;
					return($query,$hit,($Line[7]-$Line[6] + 1),($Line[9]-$Line[8] + 1),$dist); # query, hit, alignment length of query, alignment length of hit, bit score ratio= bit score / maximum bit score of blast (query against itself)
				}
			}        		
			else 
			{
				print "Error No max score for $query ($spA et $spB)\n";
				exit();
			}
		}
	}
    return ($query,$hit,0,0,0);
}
