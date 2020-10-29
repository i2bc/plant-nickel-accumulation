#!/usr/bin/perl -w
# changement de format des fichiers d'output de BRH pour input dans MARIO
# CD petites modifications du script original pour un output identique
# <!> Dans le cas de groupes avec des proteines en commun, un des groupes est totalement ecarte , d'ou la suppression de groupes d'orthologues potentiels
if ($#ARGV<0) { print "brh2mario.pl <repertoire de travail>\n"; exit();}
my $work_dir = $ARGV [0];
$brh_file = "$work_dir/BRH/STEP3/brh_orthologous_groups.tab";
mkdir ("$work_dir/mario_converted") unless (-e "$work_dir/mario_converted");
$brh_mario_converted_file = "$work_dir/mario_converted/brh_groups";
parserBRH($brh_file , $brh_mario_converted_file);

sub parserBRH
{
	($infile, $outfile) = @_;
	open IN, $infile;
	# lecture des groupes au format de BRH
	while ($l = <IN>)
	{	 
		#ex de ligne du fichier :canbrvCABR.00029-j161	cancavCACA.00015-j229	canglgnlGLVCAGL0G06666g	cannivCANI.00004-j76	klubavKLBA.00010-j91	kludevKLDE.00012-j221	sacceSCRT_03239	7
		chomp($l);
		# Recuperation des differentes colonnes dans un tableau
		@group = split(/\s+/,$l); 
		$size = pop(@group); # discard the last field (= group size)
		# lecture des proteines du groupe courrant
		foreach $prot (@group) 
		{
			# Discard empty field ( filled with -)
			if ($prot ne '-') 
			{
				# On enregistre le contenu des groupes dans un hash (numero de ligne/ID du groupe => proteines)
				# On enregistre les hash  dans %group_content (ID group => (groupe => proteines))
				$group_content{$.}{$prot} = '1'; 
			}
		}
		# on enregistre dans un hash (effectif groupe => numero de ligne/ID des groupes de proteines) les lignes de même effectif
		# on enregistre  les  effectifs comme clefs du hash % groupe_size (effectif => hash (effectifs=> ID groupes) )
		$group_size{$size}{$.} = '1'; 
	}
	close(IN);
  
	# Avec la methode BRH, une sequence peut appartenir à  plusieurs groupes, on ne va retenir que le groupe dont la taille est la plus importante
	print "Nettoyage BRH\n";

	# Tri des effectifs des groupe par taille decroissante (max= nombre d'especes)
	@sizes = keys(%group_size);  # liste des effectifs des groupes
	@sizes = sort { $b <=> $a } @sizes; # liste des effectifs classés par ordre decroissant

	open OUT,  ">".$outfile;
	# On parcourt la liste des effectifs classés par ordre decroissant
	for ($i = 0; $i <= $#sizes; $i++) 
	{
		# Pour chaque $n (= groupe ID) d'effectif $i 
		foreach $n (keys(%{$group_size{$sizes[$i]}}))
		{
			$test = 0; # Booleen qui vaut 1 si une sequence du groupe a deja ete vue 
			@group_proteins = keys(%{$group_content{$n}}); # Liste des proteines du groupe $n
			print "protein_list: \n";
			foreach $c (@group_proteins) { print $c." as c\t";} ; print "\n";# affiche les proteines du groupes courrant
			
			# parcourt de la liste des proteines d'un groupe d'effectif $i
			$j = 0; 
			# si aucune proteine du groupe n'a déjà été vue dans un groupe de taille plus grande (ou egale) ($test == 0) on conserve le groupe, 
			# sinon ($test == 1) le groupe sera ecarte (<!>aleatoirement si l'effectif des deux groupes est le même)
			while (($test == 0) && ($j <= $#group_proteins))
			{
				# pour chaque proteine du groupe courrant
				print "while $j protein_list[j] is ".$group_proteins[$j]."\n";
				
				# hash %protein_list (protein =>  ID group?				
				if (exists($protein_list{$group_proteins[$j]})) { $test = 1 ; }
				else
				{
					$protein_list{$j} = '';
	 				$j++;
				}
			}
			if ($test == 0)
			{
				# ajout des proteines du groupe courrant dans le hash des proteines déjà enregistrée lors de l'analyse d'un groupe precedent
				# % protein_list ( protein des groupes finaux => ligne du fichier input)
				foreach $k (keys(%{$group_content{$n}})) 
				{
					$protein_list{$k} = $n; # Enregistrement des IDs conservés pour ne pas les enregistrer une deuxieme fois
				}
				# ecriture de la ligne au format MARIO (separation par des ;)
				# les proteines d'un groupe sont les clefs du hash  $group_content{numero de ligne}
				$res=join(";",keys(%{$group_content{$n}}));
				print OUT $res."\n";
			}
		}			
	}
	close(OUT);
}



