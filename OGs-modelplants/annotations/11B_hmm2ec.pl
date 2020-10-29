#!/usr/bin/perl -w
use strict;
#use DBI;
use Array::Utils qw(:all);
# Description
# Annotation des groupes d'orthologues
# diff avec la version 11_intersection.pl : permet de conserver la sequence de SPMC dont l'annotation a été transferee


# Parametres
if ($#ARGV <1) { print "perl 11_hmm2ec.pl <dossierGIT> <dossier des groupes>\n";exit();} 
my $work_dir = $ARGV[0]; # dossier GIT
my $mario_dir = $ARGV[1];# dossier contenant les groupes non annotés (output de MaRiO)


# input
my $references_dir = $work_dir.'/annotations/annot_files/'; # dossier contenant les sequences de reference avec EC de  SwissProt/MetaCyc
my $annotated_sequences_fungipath_tab = $references_dir.'annotated_sequences_fungipath.tab'; # INPUT : Fichier contenant les sequences similaires entre FUNGIpath et SwissProt/MetaCyc
my $hmmsearch_dir = $work_dir.'annotations/group_EC/'; # INPUT : repertoire contenant les resultats de hmmsearch pour chaque groupe
my $outfile_dir = $work_dir.'/annotations/EC_annotated_group-v1/';

system("mkdir -p $outfile_dir");

# output
my $SwissProt_MetaCyc_ec_sql =  $references_dir.'SwissProt_MetaCyc_ec.sql';  # INPUT : fichier des sequences de reference dans SwissProt/MetaCyc (tabular)
my $SwissProt_MetaCyc_ec_fa = $references_dir.'SwissProt_MetaCyc_ec.fa'; # INPUT : fichier des sequences de reference dans SwissProt/MetaCyc (fasta)
my $NoGroupEC_sql =  $outfile_dir.'NoGroupEC.sql'; # OUTPUT : sequences depourvues d'orthologues mais annotees dans SwissProt/MetaCyc 
my $GroupEC_sql=  $outfile_dir.'GroupEC.sql';  # OUTPUT : annotation des groupes d'orthologues 
my $ID_EC_Predit_sql =  $outfile_dir.'ID_EC_Predit.sql';  # OUTPUT : annotation de chaque sequence presente dans un groupe
my $group_ref_file =  $outfile_dir.'group_ref.tab';  # OUTPUT : annotation de chaque sequence presente dans un groupe


my $seuilBestEval = 1e-45; # 1er seuil sur la Evalue
my $seuilMinEval = 1e-20; # 2nd seuil sur la Evalue

# Construction d'un hash entre les display_ids des references et leurs séquences protéiques
open(IN,$SwissProt_MetaCyc_ec_fa) or die("$SwissProt_MetaCyc_ec_fa not found");
my %ref_seq = ();
my $id = "none";
my $seq = "";
while( my $line = <IN>)
{
  chomp($line);
  if($line =~ /^>/)
  { 
	 if($id ne "none")
	 { 
		$ref_seq{$id} = $seq; #store last sequence
	 }
	$line =~ s/>//;
	my @line = split(/\|/,$line);
	$id = $line[0];
	$seq = "";
	}
  else
  {
    $seq .= $line;
  }
}
$ref_seq{$id} = $seq;
close(IN);

# Recuperation des EC numbers associés a des sequences de reference de SwissProt et/ou MetaCyc
my %ref_ecs;
if (-e $SwissProt_MetaCyc_ec_fa)
{
	open IN,$SwissProt_MetaCyc_ec_fa or die;
	while(my $line = <IN>) #AMPA_CHRSD	3.4.11.1,3.4.11.10	 GO:0004177, IEA:UniProtKB-HAMAP.;...
	{
		if ($line =~ /^>/)
		{
			$line =~ s/>//;
			my @line = split(/\|/,$line);
			$ref_ecs{$line[0]} = $line[1];
		}
	}
	close(IN);
}
else
{
	die ("Error : No file $SwissProt_MetaCyc_ec_fa\n");
}

# Recuperation des sequences de FUNGIpath annotatees dans SwissProt et/ou MetaCyc
my %AnnotatedSequences_FUNGIpath;
my %aliasID;

if (-e $annotated_sequences_fungipath_tab) 
{
    open IN,$annotated_sequences_fungipath_tab or die;
    while (my $line = <IN>) 
    {       
        chomp($line);
        my @line = split(/\t/,$line);
        
        if ($line[2] ne '-') 
        { 
			$ref_ecs{$line[1]} = $line[2];
			$aliasID{$line[1]} = $line[0];			
		}
    }
    close(IN);
}
else 
{
    die ("Error : file $annotated_sequences_fungipath_tab not found\n");
}

#test
#foreach my $k (keys(%{$ref_ec{'fungipath'}}))
#{
	#foreach my $l (keys(%{$ref_ec{'fungipath'}{$k}}))
	#{
		#print "proteine:EC number".$k.':'.$l."\n";
	#}
#}
#exit();

# Creation des fichiers OUTPUT
open GROUPEC,'>'.$GroupEC_sql;
open IDEC,'>'.$ID_EC_Predit_sql;
open GROUPREF,'>'.$group_ref_file;
# analyse des resultats de hmmsearch
chdir($hmmsearch_dir);
my @File = glob("EC*.hmmsearch");

# Pour chaque fichier résultat
foreach my $file (@File) 
{
	#if ($file =~ /^EC([0-9]+)\.hmmsearch$/ and $1 == 3039) 
	if ($file =~ /^EC([0-9]+)\.hmmsearch$/ ) 
	{ 
		my $group_id = $1 ; # recuperation du numero du groupe courrant		
		#recuperation des ID des sequences du groupe dans FUNGIpath 
		my @group_orthologs = ();
		open(IN,"$mario_dir/$group_id.fa") or die("impossible d'ouvrir le fichier de sequence du groupe");
		while(my $l= <IN>)
		{
			if($l =~ />(.+)\t.+/)
			{
				push(@group_orthologs,$1);
			}
		}
		close(IN);		
		#foreach my $p (@group_orthologs) { print $p."\n";}
		#print "//\n";
		
		#recuperation des ID des sequences du group	
							
		# Ouverture et lecture du fichier résultat de hmmsearch
		open IN,$hmmsearch_dir.$file or die("impossible d'ouvrir $hmmsearch_dir/$file");
		open HMM, $hmmsearch_dir.$file or die("impossible d'ouvrir $hmmsearch_dir/$file");
			
		my @hits = ();
		my @best_hits = ();			
		my @group_ecs;
		my @refs_ecs;
		my $evalue = 1e-45;
		while (my $l = <HMM>) 
		{
			#print "current$l";
			chomp $l;				
			unless ($l =~ /^#/)
			{
				my @results = split(/ +/,$l);										
				if ($results[4] <= $evalue)
				{
					#print "res: ".$results[4]."\n";
					push (@hits, $l);
					#foreach my $case (@hits) { print $case."\n";} print "\n";
					$evalue = $results[4];
				}
				else {last;}	
			}
			#print $#hits."\n";
			if($#hits == 0)
			{
				@best_hits = ($hits[0]);
				
			}
			elsif ($#hits >0)
			{
				my $evalue_dom = 100;
				my $score_dom = 0;
				foreach my $h (@hits)
				{
					my @h = split(/ +/,$h);
					#print $h[7]."\n";
					if ($h[7] < $evalue_dom)
					{
						@best_hits = ($h);
						$evalue_dom = $h[7];
						$score_dom = $h[8];
					}
					elsif ($h[7] == $evalue_dom)
					{						
						if( $h[8] > $score_dom)
						{
							@best_hits = ($h);
							$score_dom = $h[8];
						}
						elsif ($h[8] == $score_dom)
						{
							push( @best_hits, $h);
						}
					}
				}
			}
			if (scalar(@best_hits) >= 1) 
			{
				foreach my $c (@best_hits) 
				{ 
					my @res = split (/ +/,$c);
					my @hitname = split (/\|/, $res[0]);
					my $gec = $ref_ecs{$hitname[0]};	
					my @gec = split (/,/,$gec);
					foreach my $ec (@gec) 
					{ 
						if (is_not_in ($ec,@group_ecs) )
						{ 	
								push (@group_ecs, $ec);
						}
					}						
				}
				
				
			}
		}
		
		foreach my $o (@group_orthologs)
		{
			if (exists($aliasID{$o})) # la proteine est une reference
			{
				print $o.'->'.$aliasID{$o}."\n";
				if (exists($ref_ecs{$aliasID{$o}}))
				{								
					print "Reference EC numbers:".$ref_ecs{$aliasID{$o}}."\n";
					my @rec = split(/,| /, $ref_ecs{$aliasID{$o}});
					foreach my $ec (@rec)
					{
						print GROUPREF  $group_id.';'.$o.';'.$ec.';'; 
						if(is_not_in ($ec,@group_ecs) )	{ print GROUPREF "conflict\n";} else { print GROUPREF "consensus\n";}
						if (is_not_in ($ec,@refs_ecs) )
						{	
							push (@refs_ecs, $ec);
						}	
					}					
				}								
			}
		}
		if (scalar(@hits) >0)
		{	
			# group
			print "group: $group_id\n";
			print "nb hits: ". scalar(@best_hits)."\n";
			# EC predictions
			print "group EC: \n";
			foreach my $ec (@group_ecs) 
			{ 
				print $ec.'|';
				print GROUPEC "$group_id;$ec;";
				if (is_not_in($ec, @refs_ecs)) 
				{ 
					print "predicted\n";
					print GROUPEC  "predicted\n";
				} 
				else 
				{ 
					print "referenced\n";
					print GROUPEC "referenced\n";
				}
			}
			# EC in group reference 
			print "Ref EC: \n";
			my $ref_ec_field;			
			foreach my $ec (@refs_ecs) 
			{ 
				$ref_ec_field .= $ec.'|';
			}
			if (defined($ref_ec_field)) {chop $ref_ec_field;}
			if(defined ($ref_ec_field))
			{
				print $ref_ec_field."\n";
			}
			# comment about the EC number prediction
			if (scalar(@refs_ecs == 0))
			{
				print "predit\n";
			}	
			elsif (array_diff(@refs_ecs, @group_ecs) ) 
			{
				print "conflict\n";			
			}
			else 
			{
				print "consensus\n";
			}
			print "//\n";
		}
	}
}						
			
sub is_not_in
{
	my ($item, @list) = @_;
	foreach my $case (@list)
	{
		if ($case eq $item) {return 0;}
	}
	return 1;
}
