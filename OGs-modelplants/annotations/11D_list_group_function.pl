#!/usr/bin/perl -w
use warnings;

#extraction des infos sur la fonction de chaque groupe à partir des résultats de hmmsearch (hits des profils des groupes dans metacyc.fa et swissprot.fa)
# output dans un fichier par groupe
if ($#ARGV <0) { print "usage: perl 11_list_group_function.pl <git directory> [run tag]\n"; exit();}

$workdir = $ARGV[0];
$tag = $ARGV[1]; 
if ($tag) { $tag= '_'.$tag;} else { $tag = '';}

$hmm_dir = "$workdir/annotations/group_function$tag/";
$annot_dir = "$workdir/annotations/annotated_group_reload$tag/";
if (-e "$workdir/annotations/annotated_group_reload$tag/") {} else { system "mkdir $workdir/annotations/annotated_group_reload$tag/";}
# open(OUT,">$annot_dir/group_function.sql") or die("can not open $annot_dir/group_function.sql");
# my $cutoff = 1e-20;


#my @cutoff = (1e-10, 1e-15, 1e-20 , 1e-25, 1e-30, 1e-35, 1e-40, 1e-45, 1e-50, 1e-55,1e-60 , 1e-65, 1e-70, 1e-75, 1e-80, 1e-85, 1e-90, 1e-95);
open (BESTMC, ">besthits_metacyc.tab");
open (BESTSP, ">besthits_swissprot.tab");
my @cutoff = (1e-45);

foreach my $cutoff (@cutoff)
{
	unlink("$annot_dir/*_functions\_$cutoff");
	#parsing of hmmsearch results against (metacyc case)
	chdir($hmm_dir);
	@metacyc_files=glob("MC*");
	foreach $mcf(@metacyc_files)
	{
		my $group;
		if($mcf=~/MC(\d+)\.hmmsearch/) {$group=$1;} else { die ("missing MetaCyc infiles");}
	
		my @hits = ();
		my @desc = ();
		my @evalues = ();
		my $nb_hits = -1;
		open(IN,$mcf)or die;
		while($line = <IN>)
		{
			#0:target name - 1:accession - 2:query name - 3:accession - 4:E-value - 5:score - 6:bias - 7:E-value - 8:score - 9:bias - 10:exp - 11:reg - 12:clu - 13:ov - 14:env - 15:dom - 16:rep - 17:inc - 18:description of target
			unless($line =~ /^#/)
			{
				chomp ($line);
				@line = split(/metacyc:/, $line);
				@hmm = split(/\s+/,$line[0]);		
				if ($cutoff >= $hmm[4]) # test on the "full sequence E-value" field
				{
					$hit_ID = $hmm[0];
					$hit_desc = $line[1];
					push (@hits, $hit_ID); 
					push (@desc, $hit_desc);
					push (@evalues, $hmm[4]);
				}
			}
		}
		$nb_hits = $#hits; 
		if ($nb_hits != -1) 
		{
			print BESTMC "$group;$hits[0];$desc[0];$evalues[0]\n"; 
			#print ">$annot_dir/$group\_functions_$cutoff\n";
			open(OUT,">$annot_dir/$group\_functions_$cutoff") or die("error creating outfile");
			my %ec_count;
			my %func_count;
			foreach $d (@desc)
			{
				@d = split (';', $d);
				$metacyc_function = $d[0]; 
				$otherDB = $d[1]; 
				@otherDB = split (':', $otherDB);
				$source = $otherDB[0];
				$source_function = $otherDB[1];
				if (lc($metacyc_function) eq lc($source_function)) 
				{ 
					print OUT $metacyc_function."[allmetacyc] \n";
				}
				elsif ( $metacyc_function eq 'none')
				{
					print OUT $source_function."[metacyc2$source]\n";
				}				
				else 
				{ 
					print OUT $metacyc_function."[metacyc]\n".$source_function."[metacyc2$source]\n";
				}						
			}
			close OUT;
		}
		else { print BESTMC "$group;-;-;-;-;100\n";}
	}
	close(BESTMC);		
	@swissprot_files=glob("SP*");
	foreach $swf(@swissprot_files)
	{
		my $group;
		if($swf=~/SP(\d+)\.hmmsearch/) {$group=$1;} else { die ("missing Swiss-Prot infiles");}
		my @hits = ();
		my @desc = ();
		my @evalues = ();
		my $nb_hits = -1;
		open(IN,$swf)or die;
		while($line = <IN>)
		{
			#0:target name - 1:accession - 2:query name - 3:accession - 4:E-value - 5:score - 6:bias - 7:E-value - 8:score - 9:bias - 10:exp - 11:reg - 12:clu - 13:ov - 14:env - 15:dom - 16:rep - 17:inc - 18:description of target
			unless($line =~ /^#/)
			{
				chomp ($line);
				@line = split(/swissprot:/, $line);
				@hmm = split(/\s+/,$line[0]);		
				if ($cutoff >= $hmm[4]) # test sur le champ "full sequence E-value"
				{
					$hit_ID = $hmm[0];
					$hit_desc = $line[1];
					push (@hits, $hit_ID); 
					push (@desc, $hit_desc);
					push (@evalues, $hmm[4]);
				}
			}
		}
		$nb_hits = $#hits; 
		if ($nb_hits != -1)
		{ 
			print BESTSP "$group;$hits[0];$desc[0];$evalues[0]\n"; 
			open(OUT,">>$annot_dir/$group\_functions_$cutoff") or die("error creating outfile");
		
			foreach $swissprot_function (@desc)
			{
				print OUT "$swissprot_function\[swissprot]\n";			
			}		
			close OUT;		
		}
		else { print BESTSP "$group;-;-;100\n";}                
	}
	close BESTSP;
	close IN;
}
