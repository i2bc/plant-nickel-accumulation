#!/usr/bin/perl -w
use threads;
use File::Basename;
# Description:
# parser les sorties de blast de maniere a avoir des fichiers au bon format pour inparanoid
# format tabulé avec les informations suivantes:
# * Query id
# * Hit id
# * Bit score
# * Query length
# * Hit length
# * Length of longest segment on query. This is defined as the length of the segment from the
#  	first position od the first hsp to the last position of the last hsp. I.e. if the
# 	hsps are 1-10 and 90-100, the length of the longest segment will be 100.
# * Length of longest segment on hit, see explanation above.
# * Total match length on query. This is defined as the total length of all hsps. If the hsps
# 	are 1-10 and 90-100, the length will be 20.
# * Total match length on hit. see explanation above
# * Positions for all segments on the query and on the hit. Positions on the query is written
# 	as p:1-100 and positions on the hit is writteh as h:1-100. Positions for all hsps are
#	  specified separated nt tab. Positions for query and hut for a segment is separated by 
# 	one whitespace. 

#score_cutoff utilisé par default par inparanoid: 40
#format blast d'entree : tabulé
#qseqid0 sseqid1 pident2 length3 mismatch4 gapopen5 qstart6 qend7 sstart8 send9 evalue10 bitscore11 qlen12 slen13

#17/06/2014
#perl parser_pour_Inparanoid.pl [score_cutoff] [dossier contenant les blasts] [dossier de sortie]
#perl parser_pour_Inparanoid.pl 40 /home/cecile/Bureau/miseajourBD/local/result/blastp/ /home/cecile/Bureau/miseajourBD/local/Inparanoid/inparanoid_cecile/result/blast/
#perl parser_pour_Inparanoid.pl 40 /home/cecile/Bureau/miseajourBD/local/Inparanoid/inparanoid_cecile/test/ /home/cecile/Bureau/miseajourBD/local/Inparanoid/inparanoid_cecile/test/

if ($#ARGV <3) 
{ 
	print "<score limite = 40 par default> <repertoire des blasts> <tag> <cpu>\n";
	exit();
}
$score_cutoff=$ARGV[0];
$dirBlast=$ARGV[1];
$dirBlast .= '/';
$tag=$ARGV[2];
$dirOut = 'STEP1_'.$tag.'/';
$threadsAllowed=$ARGV[3];
@jobs=();
@running = 
mkdir($dirOut)unless(-e $dirOut);

# on recupere tous les noms des fichiers blasts
opendir(BL,$dirBlast);
@fb=readdir(BL);
#pour chaque fichier blast

foreach $f(@fb)
{
	my $name = basename($f); 
	my @name= split (/_vs_|\./, $name);
	my $genomeA = $name[0]; # infile genomeA_vs_genomeB.blast
	my $genomeB = $name[1];
	my $processed_file = "$dirOut/$genomeA\.fa-$genomeB\.fa"; 
	if (-e $processed_file) 
	{
		print $f, $processed_file."\n";
	} 
	unless(($f eq "." )||($f eq "..") or (-e $processed_file))
	{   		
		push(@jobs,threads->create(sub{system("perl my_parser_pour_Inparanoid_job.pl $score_cutoff $dirBlast $dirOut $f");}));
		sleep(1) while(threads->list(threads::running)>=$threadsAllowed);		
	}
	@running = threads->list(threads::running);
}
closedir(BL);
#$_->join for @jobs;

while (scalar @running != 0) 
{
	foreach my $job (@jobs) 
	{
		$job->join if ($job->is_joinable());
	}
	@running = threads->list(threads::running);
}
