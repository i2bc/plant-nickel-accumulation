#!/usr/bin/perl -w

#Programme qui permet d'ecrire le fichier d'entree au format attendu par MCL et d'executer ensuite le programme MCL
#les resultats generes par MCL (i.e. les nouvelles familles) sont reformates au format de depart
use strict;

#############
# Declaration
#############
die ("8 arguments doivent etre donnes en entree (".($#ARGV+1)." ont ete donnes).\n") unless ($#ARGV == 7);

my $group = $ARGV[0];
my $dirMclIn = $ARGV[1];
my $dirMclOut = $ARGV[2]; 
my $dirGroup = $ARGV[3]; 
my $nbProtMax =  $ARGV[4]; 
my $I_para =  $ARGV[5]; # Parametre dinflation a utilise pour le programme MCL
my $dirWork = $ARGV[6]; #repertoire de travail
my $dirTmp = $ARGV[7]; #repertoire de travail

my @Line = ();
my $i = 0;
my @FileFn = ();
my %ToReDo = ();
my @Id = ();
my %Edge = ();

#si le fichier d'entree existe bien
if (-e $dirMclIn.$group) {
	#lecture du fichier pour le formater au format attendu par MCL
	open GROUP, $dirMclIn.$group; 
	open DATAMCL, ">".$dirMclIn.$group.".mcl"; #ecriture du fichier data pour mcl
	while (<GROUP>) {
		if ($_ =~ /\(edge ([^ ]*) ([^ ]*) ([^ ]*)\)/) {
			print DATAMCL "$1 $2 $3\n";;
                        $Edge{$1}{$2} = $_;
		}
	}
	close (GROUP);
	close (DATAMCL);

	#lancement de MCL
	system("/home/cpereira/Phylogeny/bin/mcl ".$dirMclIn.$group.".mcl -I ".$I_para." --abc -o ".$dirMclOut."mcl_".$group.".out");

	#lecture du resultat genere par MCL
	chdir($dirMclOut);
	@FileFn = glob("mcl*$group*.out");
	#pour chaque groupe genere par MCL
	foreach my $fileFn (@FileFn) { 
		$i = 0; #reinitilisation

		#lecture du groupe
		open IN,$dirMclOut.$fileFn;
		$fileFn =~ s/.out//;
		$fileFn =~ s/mcl_//;
		while (<IN>) {
			chomp;
			@Line = split(/\s/,$_);
			if ($#Line > 0) {
			    #si tjs superieur au seuil, on realise mcl avec seuil max si non deja fait
			    if ($#Line >= $nbProtMax) {
					print scalar(@Line)." prot ".$fileFn."-mcl-".$i."\n";
                    $ToReDo{$fileFn."-mcl-".$i} = '';                               
                }
			    print $dirMclOut.$fileFn."-mcl-".$i.": ".($#Line+1)." proteines\n";
				
			    #on reformatte le fichier genere par MCL au meme format que celui de depart
			    open OUT, ">".$dirMclOut.$fileFn."-mcl-".$i;
                $_ =~ s/\t/ /g;
			    print OUT "(nodes ".$_.")\n";
                @Id = split(/ /,$_);
			    foreach my $id1 (@Id) {
					foreach my $id2 (@Id) {
						print OUT $Edge{$id1}{$id2} if (exists($Edge{$id1}) && exists($Edge{$id1}{$id2}));
                    }
                }
			    close(OUT);
			    #on deplace le fichier dans le repertoire contenant tous les groupes car c'est ce nouveau fichier qui sera utilise par la suite
			    system("cp ".$dirMclOut.$fileFn."-mcl-".$i." ".$dirGroup );
			    $i ++;
			}
		}
		close (IN);
	}
}
else {
	print "No file ".$dirMclIn.$group."\n";
}

#suppression du fichier initial qui a ete remplace par les fichiers generes avec MCL
unlink($dirGroup.$group) if ((-e $dirGroup.$group."-mcl-0") && !(-z $dirGroup.$group."-mcl-0")); 

#si certaines familles ont encore une taille trop importante, on utilise le parametre 'inflation max
foreach my $file (keys(%ToReDo)) {
    system("cp ".$dirMclOut.$file." ".$dirMclIn);
		system('perl '.$dirWork."3_phylogeny_job_clusterLIX.pl $file $dirMclIn $dirMclOut $dirGroup $nbProtMax 5 $dirWork $dirTmp\n") if ($I_para != 5);	
}

