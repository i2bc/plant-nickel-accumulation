#!/usr/bin/perl -w
#Script transmis par Cecile Pereira avec quelques modifications (CD) ne concernant que les noms des dossiers d'input/output)
#le repertoire de travail est un repertoire cree en début de projet par git clone (il contient les scripts des 4 methodes initiales et de MaRiO)
#le tag du run est une etiquette des fichiers de sortie (permet de lancer plusieurs fois les scripts dans le même repertoire de travail
use strict;
#Construction des groupes d'orthologues pas la méthode du lien simple à partir des resultats generes par inparanoid a l'etape 2

if ($#ARGV <0)
{
	print "my3_inparanoid.pl <repertoire de travail> > > <outfile>\n";
	exit();
}
# Parametres
my $dirWork=$ARGV[0];
my $dirData = $dirWork.'/Inparanoid/STEP2/sqltable/'; # INPUT : Repertoire contenant les tables SQL generes par inparanoid entre chaque paire de Genomemes dans 2_inparanoid.pl
my $dirRes = $dirWork.'/Inparanoid/STEP3/'; # OUTPUT : Repertoire contenant les groupes d'orthologues
my $fileOutOrtho = $dirRes.'inparanoid_orthologous_groups.sql'; # OUTPUT : table sql des groupes d'ortholoques
my $fileOutInpara = $dirRes.'inparanoid_inparalogs.sql'; # OUTPUT : table sql des inparalogues

# Autres variables
my $score = 0;
my $id1 = '';
my $id2 = '';
my $Genome1 = '';
my $Genome2 = '';
my $nb = 0;
my $i = '';
my $ortho = '';
my $nbGroup = '';
my @Genome = ();
my @File = ();
my @Line = ();
my @Inpara1 = ();
my @Inpara2 = ();
my @Ortho = ();
my %PairOrtho = ();
my %IdGenome = ();
my %Inpara = ();
my %Genome = ();
my %TmpGroup = ();
my %NewId = ();
my %OldId = ();
my %IdSave = ();



############
# Programme
############

# Creation du repertoire resultat
mkdir($dirRes) unless (-e $dirRes);
# Verification que le repertoire contenant les tables SQL existe bien
die("Error in data repertory ($dirData)\n") unless(-e $dirData);

# Suppression des fichiers resultats si deje crees
unlink($fileOutInpara) if (-e $fileOutInpara);
unlink($fileOutOrtho) if (-e $fileOutOrtho);


# Recupération des résultats entre paire de génomes
###################################################
print "Lecture des comparaisons entre paire de Genomes...\n";
chdir($dirData);
@File = glob("sqltable*"); 
#Lecture des résultats Inparanoid
foreach my $file (@File) {
    if (-z $dirData.$file) {
	print "Warning : empty file ($file)\n";
    }
    else {
		# Recuperation des noms de fichiers pour verifier qu'il n'en manque aucun
		if ($file =~ /sqltable\.(.*)-(.*)$/) {#modif cecile, au depart /sqltable\.([a-z]*)-([a-z]*)$/
			my $v1=lc($1);#modif cecile
			my $v2=lc($2);#modif cecile
			$v1=~s/\.fa//;#modif cecile
			$v2=~s/\.fa//;#modif cecile
			$Genome{$v1}{$v2} = '';
			$Genome{$v2}{$v1} = '';
		}
		# Initialisation
		$nb = 0; # Numéro du groupe : au départ nul
		$id1 = '';  # ID appartenant a l'espece $Genome1
		$id2 = ''; # ID appartenant a l'espece $Genome2
		$Genome1 = ''; # espece $Genome1
		$Genome2 = ''; # espece $Genome2
		@Inpara1 = (); # Liste des Inparalogues de $Genome1 dans le groupe $nb
		@Inpara2 = (); # Liste des Inparalogues de $Genome2 dans le groupe $nb
		$score = 0; # Score du groupe $nb
		
		# Lecture du fichier
		open (IN, $dirData.'/'.$file) or die("impossible d'ouvrir le fichier $dirData.'/'.$file");
		while (<IN>) {
			#ligne type du fichier :1	4959	canbr.fa	1.000	canbrvCABR.00021-j221
			chomp($_);
			$_ =~s/\.fa//g;
			#Récupération de chaque colonne dans un tableau
			#Col 0 : Numéro du groupe
			#Col 1 : Score Blast
			#Col 2 : Génome de la séquence
			#Col 3 : Confidence value 
			#Col 4 : ID		
			@Line = split(/\t/,$_);

			#A chaque ID on lui associe son génome si non déjà fait
			$IdGenome{$Line[4]} = lc($Line[2]) unless (exists($IdGenome{$Line[4]}));
		
			#Si nouveau groupe, on réinitialise les valeurs du groupe
			if ($nb  != $Line[0]) {
				#s'il ne s'agit pas du premier groupe,on enregistre les resultats du groupe precedant
				if ($nb != 0) {
					#pour chaque paire du groupe, on enregistre le score
					#orthologues
					$PairOrtho{$id1}{$id2} = $score;  
					$PairOrtho{$id2}{$id1} = $score;
					#inparalogues a $id1
					foreach $i (@Inpara1) {
						$Inpara{$id1}{$i}{lc($Genome2)} = '';
					}
					#inparalogues a $id2
					foreach $i (@Inpara2) {
						$Inpara{$id2}{$i}{lc($Genome1)} = '';
					}
				}
				
				# Reinitialisation des valeurs pour le nouveau groupe
				$id2 = '';
				$Genome2 = '';
				@Inpara1 = ();
				@Inpara2 = ();
				$Genome1 = $Line[2]; #Genome
				$score = $Line[1]; #score
				$id1 = $Line[4]; #id
			}
			else { #sinon
				# si meme Genomeme que la premiere seq, on enregistre la seq comme inparalogue
				if ($Line[2] eq  $Genome1) {
					push(@Inpara1,$Line[4]) if ($Line[3] == 1); # On ne garde que les inparaloques ave une "confidence value" maximale i.e. 1
				}
				else {
					#sinon si deuxieme Genomeme
					$Genome2 = $Line[2];
					#si c'est la premiere seq de $Genome2, on l'enregistre comme orthologue
					if ($id2 eq '') {
						$id2 = $Line[4];
					}
					else { #  on l'enregistre comme inparalogue
						push(@Inpara2,$Line[4]) if ($Line[3] == 1);
					}
				}
			}
			#ID groupe
			$nb = $Line[0];
		}
		close(IN);
		# Enregistrement des resultats du dernier groupe
		$PairOrtho{$id1}{$id2} = $score;
		$PairOrtho{$id2}{$id1} = $score;
		foreach $i (@Inpara1) {
			$Inpara{$id1}{$i}{lc($Genome2)} = '';
		}
		foreach $i (@Inpara2) {
			$Inpara{$id2}{$i}{lc($Genome1)} = '';
		}
    }
}

# Verification que l'on a bien tous les fichiers avant de construire les groupes
#################################################################################
@Genome = keys(%Genome);
@Genome = sort(@Genome);
$i = $#Genome + 1;
print "$i Genomes\n";
foreach $Genome1 (@Genome) {
	die("Lack ".($i - scalar(keys(%{$Genome{$Genome1}})) -1 )." file(s) with $Genome1\n") if (scalar(keys(%{$Genome{$Genome1}})) != ($i-1));
}


# Les sequences inparalogues peuvent varier selon les comparaisons
# Seul les inparalogues trouves dans toutes les comparaisons sont conservees
############################################################################

# Pour chaque orthologue $id1 posssedant au moins une sequence inparalogue
@Ortho = keys(%Inpara);
foreach $id1 (@Ortho) {
	# Si la sequence n'a pas deja ete traitee
	if (exists($Inpara{$id1})){
		#pour chaque sequence inparalogue a $id1
		foreach $id2 (keys(%{$Inpara{$id1}})) {
			$i = 0; # vaudra 1 si id1 et id2 ne sont pas inparalogues
			# Si id2 n'a pas ete trouve comme inparalogue a chaque fois que $id1 a ete trouve comme orthologue, 
			# les deux sequences ne seront pas considérés comme inparalogues ($i = 1)
			$i = 1 unless (keys(%{$Inpara{$id1}{$id2}}) == keys(%{$PairOrtho{$id1}}));
			# De meme, si $id2 a ete trouve comme orthologue (exists($PairOrtho{$id2}))
			if (exists($PairOrtho{$id2})){
				# Si $id1 a bien ete trouve au moins une fois comme inparalogue a $id2
				if (exists($Inpara{$id1}) && exists($Inpara{$id2}{$id1})){
					# On verifie que $id2 a ete trouve comme etant inparalogue a $id1 autant de fois que $id2 a ete trouve comme orthologues
					$i = 1 unless ( keys(%{$PairOrtho{$id2}}) == keys(%{$Inpara{$id1}{$id2}}) );
				}
				else{ # si $id2 n'a jamais ete trouve comme inparalogues
					$i = 1; # $id1 et $id2 ne sont pas inparalogues
				}
			}
			# if (exists($PairOrtho{$id2}) && exists($Inpara{$id1}) && exists($Inpara{$id2}{$id1})){
			  # $i = 1 unless ( keys(%{$PairOrtho{$id2}}) == keys(%{$Inpara{$id1}{$id2}}) );
			# }

			#si les inparalogues n'ont pas ete trouvés a chaque fois, on supprime la relation d'inparalogues entre les deux sequences
			if ($i == 1) {
			   delete($Inpara{$id1}{$id2}) if (exists($Inpara{$id1}) && exists($Inpara{$id1}{$id2}));
			   delete($Inpara{$id2}{$id1}) if (exists($Inpara{$id2}) && exists($Inpara{$id2}{$id1}));
			   delete ($Inpara{$id1}) if (keys(%{$Inpara{$id1}}) == 0);
			   delete ($Inpara{$id2}) if (keys(%{$Inpara{$id2}}) == 0);
			}
		}
	}
}

# On enregistre toutes les sequences inparalogues comme orthologues car n'a plus d'importance pour la construction des groupes
foreach $id1 (keys(%Inpara)) {
    foreach $id2 (keys(%{$Inpara{$id1}})) {
        $PairOrtho{$id2}{$id1} = '';
        $PairOrtho{$id1}{$id2} = '';
        foreach $ortho (keys(%{$PairOrtho{$id1}})) {
            $PairOrtho{$id2}{$ortho} = $PairOrtho{$id1}{$ortho};
            $PairOrtho{$ortho}{$id2} = $PairOrtho{$id1}{$ortho};
        }
    }
}


print "Construction des groupes...\n";
# Construction des groupes par liens simples 
############################################
# pour chaque sequence ortho
foreach $id1 (keys(%PairOrtho)) {
    %TmpGroup = ();
    %NewId = ();
    # A moins que la sequence ait deja ete traite (car une seq de son groupe a deja ete analysee)
    unless (exists($IdSave{$id1})) {
	$IdSave{$id1} = '';

        # Initialisation du groupe avec $id1
        $TmpGroup{$id1} = '';

        # On recupere tous les ortholoques a $id1
        foreach $id2 (keys(%{$PairOrtho{$id1}})) {
            $NewId{$id2} = '' unless (exists($TmpGroup{$id2}));
        }

        # On recupere tous les ortholoques des ortholoques tant qu'on peut
        while (keys(%NewId) != 0) {
            %OldId = %NewId; # orthologues precedemment ajouter
            %NewId = (); # nouveaux orthologues 
			
            # Enregistrement des ID à ajouter dans le groupe
            foreach $id2 (keys(%OldId)) {
                $TmpGroup{$id2} = '';
                $IdSave{$id2} = '';
            } 
			
            # Enregistrement des orthologues aux sequences precedemment ajoutees si la sequence n'appartient pas encore au groupe
            foreach $id2 (keys(%OldId)) {
                foreach $ortho (keys(%{$PairOrtho{$id2}})) {
                    $NewId{$ortho} = '' unless (exists($TmpGroup{$ortho}));
                }
            }
        }

        print "GROUPE $nbGroup\t$id1\t".keys(%TmpGroup)." sequences (initially)\n";
        $nbGroup = searchGroup(\%TmpGroup,\%PairOrtho,\%IdGenome,\@Genome,$nbGroup,$fileOutOrtho,$fileOutInpara);
    }
}

print "Instruction pour la base de données : \n";
print "CREATE TABLE inparanoidortho (";
#a chaque Genomeme on associe un numero qui lui est specifique
for ($i = 0; $i <= $#Genome; $i ++) {
    print "$Genome[$i] character varying(80), ";
}
print "size integer PRIMARY key)\n";
print "CREATE TABLE inparanoidinpara (inpara character varying(80), nb integer  REFERENCES inparanoidortho (nb), CONSTRAINT inparanoidinpara_pkey PRIMARY KEY (inpara,nb));\n";
print "\\COPY inparanoidortho FROM '".$fileOutOrtho."' WITH DELIMITER AS '\\t'\n";
print "\\COPY inparanoidinpara FROM '".$fileOutInpara."' WITH DELIMITER AS '\\t'\n";

############
# Fonction
############

sub searchGroup {
    my ($refTmpGroup,$refPairOrtho,$refIdGenome,$refGenome,$nbGroup,$fileOutOrtho,$fileOutInpara) = @_;
		my $nbId = keys(%$refTmpGroup); # Nombre de sequences dans le groupe
    my $min = 100;
    my $id1 = '';
    my $id2 = '';
    my $perc = 0;
		my $i = 0;
    my $j = 0;
		my $Genome = '';
		my %Link = ();
    my %TmpGroupNew = ();
    my %IdSave = ();
    my %NewId = ();
    my %OldId = ();
    my %Ortho = ();
    my @Ortho = ();
    my @Inpara = ();
    my @Min = ();
    my $seuilpourcentageliengroupe=20;
	
    # Pour chaque seq, on calcule le nombre de liens de la seq / au nombre total de seq dans le groupe
    foreach $id1 (keys(%$refTmpGroup)) {
        $perc = int(keys(%{$$refPairOrtho{$id1}})*100/($nbId-1));
        $min = $perc if ($min > $perc); # valeur minimale
    }
    #si au minimum les proteines ont au moins $seuilpourcentageliengroupe% de lien entre elles, on enregistre le groupe tel quel
    if ($min >= $seuilpourcentageliengroupe) {
		#on distingue orthologues et inparalogues
        foreach $id1 (keys(%$refTmpGroup)) {
            if (exists($$refIdGenome{$id1})) {
                $Genome = $$refIdGenome{$id1}; # Genomeme de $id1
		#si une sequence a deja ete enregistree pour le Genomeme $Genome, on enregsitre $id1 comme inparalogue
                if (exists($Ortho{$Genome})) {
                    push(@Inpara,$id1);
                }
		else { #sinon orthologues
                    $Ortho{$Genome} = $id1;
                }
            }
            else {
                die( "PB with ID $id1.\n");
            }
        } 
		#on enregistre resultat dans le fichier si plus d'une espece appartient au groupe
        if (keys(%Ortho) > 1) {
            $nbGroup ++;
			# Enregistrement des inparalogues du groupe
            open INPARA,'>>'.$fileOutInpara;
            foreach $id1 (@Inpara) {
                print INPARA "$id1\t$nbGroup\n";
            }
            close(INPARA);
			# Enregsitrement des orthologues du groupe
            open ORTHO,'>>'.$fileOutOrtho;
            for ($i = 0 ; $i <= $#$refGenome; $i ++) {
                if (exists($Ortho{$$refGenome[$i]})) {
                    print ORTHO "$Ortho{$$refGenome[$i]}\t";
                }
                else {
                    print ORTHO "-\t";
                }
            }
            print ORTHO "$nbGroup\n";
            close(ORTHO);
        }
    }
	#si valeur minimale inferieure à $seuilpourcentageliengroupe, 
    else {
		# Si une sequence a tres peu de lien avec les sequences du groupe, 
		# on va supprimer les liens les moins fiables qui fusionnent des groupes proches 
		# et qui provroquent des groupes de taille tres importante
        %Link = ();
        $min = 100;
		# Pour chaque sequence du groupe
        foreach $id1 (keys(%$refTmpGroup)) {
			#recuperation de tous ces orthologues et inparalogues
            @Ortho = keys(%{$$refPairOrtho{$id1}});
			#pour chacune de ces sequences
            for ($i = 0; $i <= $#Ortho; $i++) {
				# on enregistre le lien dans un tableau 
                $Link{$id1}{$Ortho[$i]} = 0;
				# et on va compter le nombre d'orthologues communs entre $Ortho[$i] et chaque autre orthologue $Ortho[$j] de $id1
                for ($j = 0; $j <= $#Ortho; $j ++) {
                    unless ($i == $j) {
                        $Link{$id1}{$Ortho[$i]} ++ if (exists($PairOrtho{$Ortho[$i]}{$Ortho[$j]}));
                    }
                }
				#calcul du pourcentage de sequences communes 
				$perc = int(100 * ($Link{$id1}{$Ortho[$i]}) / ($#Ortho+1 + keys(%{$PairOrtho{$Ortho[$i]}}) - $Link{$id1}{$Ortho[$i]}) );
                $Link{$id1}{$Ortho[$i]} = $perc;
				# on enregistre sequences partageant le score minimal
				if ($min > $perc){
					$min = $perc;
					@Min = ();
					push(@Min,$id1);
					push(@Min,$Ortho[$i]);
				}
				elsif($min == $perc){
					push(@Min,$id1);
					push(@Min,$Ortho[$i]);	
				}
            }
        }
		#suppression du lien qui relie les paires de sequences partageant le moins d'orthologues communs (deux seq paralogues)
		for ($i = 0; $i <= $#Min; $i = $i + 2){
			$id1 = $Min[$i]; $id2 = $Min[$i+1];
			if (exists($$refPairOrtho{$id1}{$id2})){
				#print "Suppression du lien entre $id1 et $id2 (score : $$refPairOrtho{$id1}{$id2} et pourcentage $min)\n";
				delete($$refPairOrtho{$id1}{$id2});
			}
            delete($$refPairOrtho{$id2}{$id1}) if (exists($$refPairOrtho{$id2}{$id1}));

            # recontruction du groupe sans ce lien et on refait le test
            foreach $id1 (keys(%$refTmpGroup)) {
                unless (exists($IdSave{$id1})) {
                    #initialisation
                    %TmpGroupNew = ();
                    %NewId = ();
                    $TmpGroupNew{$id1} = '';
                    $IdSave{$id1} = '';
                    #on recupére tous ces ortho
                    foreach $id2 (keys(%{$$refPairOrtho{$id1}})) {
                        $NewId{$id2} = '' unless (exists($TmpGroupNew{$id2}));
                    }

                    #on récupère tous les orthos des orthos
                    while (keys(%NewId) != 0) {
                        %OldId = %NewId;
                        %NewId = ();
                        #enregistrement des id à ajouter
                        foreach $id2 (keys(%OldId)) {
                            $TmpGroupNew{$id2} = '';
                            $IdSave{$id2} = '';
                        }
                        #on verifie que les orthos des id rajoutés sont dans le groupe sinon on les rajoute
                        foreach $id2 (keys(%OldId)) {
                            foreach $i (keys(%{$$refPairOrtho{$id2}})) {
                                $NewId{$i} = '' unless (exists($TmpGroupNew{$i}));
                            }
                        }
                    }
                    $nbGroup = searchGroup(\%TmpGroupNew,\%$refPairOrtho,\%$refIdGenome,\@$refGenome,$nbGroup,$fileOutOrtho,$fileOutInpara) if (keys(%TmpGroupNew) != 1);
                }
            }
		}
    }
    return $nbGroup; #ID groupe
}
