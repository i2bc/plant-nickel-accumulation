#!/usr/bin/perl -w

#Programme qui prend en entrée un graphe dans le repertoire STEP2, le decoupe avec le logiciel MCL pour et écrit plusieurs fichiers de sous-graphe (au même format) en sortie
use strict;

my $group = $ARGV[0]; # nom du graphe/ groupe_{wildcard} dans snakemake
my $dirWork = $ARGV[1]; #repertoire de travail (git clone)
my $dirOut = $dirWork."/Phylogeny/STEP3/"; # output de l'étape (fichiers originaux non traités (copie Snakefile) + sortie de MCL)
my $dirMclIn = $dirWork."/Phylogeny/STEP3/mcl_input/"; # input MCL
my $dirMclOut = $dirWork."/Phylogeny/STEP3/mcl_output/";  # output MCL
my $nbProtMax =  $ARGV[2]; # nombre maximum de proteines dans un graphe
my $I_para =  $ARGV[3]; # Parametre d'inflation a utilise pour le programme MCL (premier run)

my @Line = ();
my $i = 0;
my @FileFn = ();
my %ToReDo = ();
my @Id = ();
my %Edge = ();

#si le fichier d'entree existe bien
if (-e $dirMclIn.$group) 
{
	
	open GROUP, $dirMclIn.$group; #lecture du fichier copier dans mcl_input
	open DATAMCL, ">".$dirMclIn.$group.".mcl"; #reformatage du fichier pour MCL (dans le même repertoire mcl_input)
	while (<GROUP>) {
		if ($_ =~ /\(edge ([^ ]*) ([^ ]*) ([^ ]*)\)/) {
			print DATAMCL "$1 $2 $3\n";;
                        $Edge{$1}{$2} = $_;
		}
	}
	close (GROUP);
	close (DATAMCL);

	# premier lancement de MCL avec l'option -I (input dans mcl_input=groupe_{wildcard}.mcl ; output dans mcl_output=mcl_groupe_{wildcard}.out
	system("mcl ".$dirMclIn.$group.".mcl -I ".$I_para." --abc -o ".$dirMclOut."mcl_".$group.".out");

	#lecture du resultat genere par MCL
	chdir($dirMclOut);
	#system ("ls");
	@FileFn = glob("mcl*$group*.out"); #lecture des fichiers de sortie de MCL (dans le repertoire mcl_output <!>)
	
	# le fichier mcl_groupe_{name}.out comprend un groupe par ligne
	foreach my $fileFn (@FileFn) 
	{ 
		$i = 0; # reinitilisation de l'indice du groupe 
		
		open IN,$dirMclOut.$fileFn;
		$fileFn =~ s/.out//;
		$fileFn =~ s/mcl_//; # recuperation du nom original du groupe groupe_{wildcard}
		while (<IN>) 
		{
			chomp;
			@Line = split(/\s/,$_); # premier groupe
			if ($#Line > 0) 
			{
			    # si le nombre de proteines est toujours superieur au seuil (nombre de genomes x 5), on realise MCL avec le seuil max si non deja fait
			    if ($#Line >= $nbProtMax) 
			    {
					print scalar(@Line)." prot ".$fileFn."-mcl-".$i."\n";
                    $ToReDo{$fileFn."-mcl-".$i} = '';   # stockage du nom du fichier dans un hash pour reload de MCL
                }
			    print $dirMclOut.$fileFn."-mcl-".$i.": ".($#Line+1)." proteines\n";
				
			    # sinon on reformate le fichier genere par MCL au format de graphe de depart groupe_{wildcard}-mcl-n°ligne
			    open OUT, ">".$dirMclOut.$fileFn."-mcl-".$i;
                $_ =~ s/\t/ /g;
			    print OUT "(nodes ".$_.")\n"; # premiere ligne = nodes
                @Id = split(/ /,$_);
			    foreach my $id1 (@Id) 
			    {
					foreach my $id2 (@Id) 
					{
						print OUT $Edge{$id1}{$id2} if (exists($Edge{$id1}) && exists($Edge{$id1}{$id2}));
                    }
                }
			    close(OUT);
			    # on copie le nouveau fichier dans le repertoire contenant tous les groupes apres MCL (Phylogeny/STEP3), input de l'etape suivante (STEP4-1) [modification pour Snakemake]
			    system("cp ".$dirMclOut.$fileFn."-mcl-".$i." ".$dirOut );
			    $i ++;
			}
		}
		close (IN);
	}
}
else 
{
	print "No file ".$dirMclIn.$group."\n";
}

#suppression du fichier initial qui a ete remplace par les fichiers generes avec MCL
unlink($dirOut.$group) if ((-e $dirOut.$group."-mcl-0") && !(-z $dirOut.$group."-mcl-0")); 

#si certaines familles ont encore une taille trop importante, on utilise le parametre 'inflation max
foreach my $file (keys(%ToReDo)) # Second run de MCL
{
    system("cp ".$dirMclOut.$file." ".$dirMclIn);
	system('perl '.$dirWork."/Phylogeny/3_phylogeny_snakemake.pl $file $dirWork $nbProtMax 5") if ($I_para != 5);	
}

