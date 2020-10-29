#!/usr/bin/perl -w
use strict;
use Data::Dumper;
# Description
# Construction des groupes d'orthologues par la méthode des liens multiples à partir des BRH entre paire de genomes
#perl 3_BRH.pl [repertoire de travail] [repertoire contenant genomes_abrege.fa ]

if ($#ARGV <1) 
{ 
	print "Usage my3_BRH.pl <repertoire de travail> <repertoire des genomes> ";
	exit();
}
my $dirWork = $ARGV[0]; #repertoire de travail
my $dirInGenome = $ARGV[1]; # INPUT : Repertoire contenant les genomes
my $dirInBrh= $dirWork."/BRH/STEP2/";
my $dirOut = $dirWork."/BRH/STEP3/";
my $fileOut = $dirOut.'brh_orthologous_groups.tab'; # OUTPUT : fichier contenant les groupes d'orthologues par liens multiples


# Autres varables
my @File = ();
my @Line = ();
my @Geno = ();
my @TmpGroup = ();
my @TmpId = ();
my @Brh_prot = ();

my %Brh2 = ();
my %Geno = ();
my %Save = ();
my %NewGroup = ();
my %AllGroup = ();


my $size = 0;
my $nb_group_with_prot = 0;
my $test = 0;
my $test1 = 0;
my $test2 = 0;
my $diff = 0;
my $ref = '';
my $tmpId = '';
my $num_group = 0;
my $id1 = '';
my $id2 = '';

#######################################################################
# Programme
###########

# Verificiation des repertoires INPUT
die("No directory $dirInGenome\n") unless (-e $dirInGenome);
die("No BRH directory $dirInBrh\n") unless (-e $dirInBrh);

# Creation du repertoire OUTPUT
mkdir($dirOut) unless (-e $dirOut);

# Recuperation des noms de genomes
##################################
if (-e $dirInGenome){
	chdir($dirInGenome);
	@File = glob("*fa");
	foreach my $file (@File) 
	{
		$file =~ s/\.fa$//;
		push(@Geno,$file);
	}
}

#tri des genomes par ordre alphabetique
@Geno = sort(@Geno);

#pour verifier qu'il ne manque pas de fichier
for (my $i = 0; $i < $#Geno; $i ++) 
{
	$Geno{$Geno[$i]} = $i; # A chaque genome on associe son ID (= position par ordre alphabetique)
    for (my $j = $i + 1; $j <= $#Geno; $j ++) 
    {
		$Save{$Geno[$i].' and '.$Geno[$j]} = '';
	}
}

# A chaque genome on associe son ID (= position par ordre alphabetique)
for (my $i = 0; $i <= $#Geno; $i ++) 
{
	$Geno{$Geno[$i]} = $i;
}

#Lecture de toutes les paires de BRH 
####################################*

# Les nouvelles
chdir($dirInBrh);
print "dirinbrh $dirInBrh\n";
@File = glob("*.brh");
foreach my $file (@File) 
{
    print "$file\n";
    if ($file =~ /(.*)_vs_(.*).brh/) {
       my $spA = $1;
       my $spB = $2;
	   delete ($Save{$spA.' and '.$spB});
       open IN, $dirInBrh.$file ;
       while (<IN>) 
       {
            chomp($_);
            @Line = split(/\t/,$_);
			# on enregistre dans une tabke de hash toutes les paires de BRH associées à une sequence
			# la valeur correspond a l'identifiant du genome de la deuxieme sequence
			# paire ID1 ID1 = genome 1
            $Brh2{$Line[0]}{$Line[0]} = $Geno{$spA};
            $Brh2{$Line[1]}{$Line[1]} = $Geno{$spB};
			# paire ID1 ID2 = genome 2
            $Brh2{$Line[0]}{$Line[1]} = $Geno{$spB};
            $Brh2{$Line[1]}{$Line[0]} = $Geno{$spA};
        }
        close(IN);
    }
}
#si il manque certains fichiers de paires de BRH, on arrete
if (keys(%Save) != 0){
	my $tmp=join(' ',keys(%Save));
	print "$tmp\n";
	print "Not all BRH files :\n".join("\n",keys(%Save))."\n";
	die ("STOP : Not all BRH files\n");
}

%Save = ();

# Construction des groupes par liens multiples
##############################################

open OUT,'>'.$fileOut;

# Pour chaque sequence $prot formant au moins une paire de BRH avec une autre sequence
foreach my $prot (keys(%Brh2)) 
{ 
    #Array des proteines ayant un BRH avec $prot trie par ordre alphabetique
    @Brh_prot = sort(keys(%{$Brh2{$prot}})); #cle2 = liste de tous les BRH de $prot
    $ref = join('',@Brh_prot); # construction d'une variable contenant une liste triées des protéines ($prot inclus)
    # enregistrement de la suite  dans %save (a moins qu'une autre protéine du groupe ai déjà été vu dans la boucle avec les memes BRHs)
    unless (exists($Save{$ref})) 
    {	
		# enregistrement de la liste $ref dans %save (a moins qu'une autre protéine du groupe ai déjà été vu dans la boucle et possede les meme BRHs)
        $Save{$ref} = '';
		# Initialisation des variables
        $test1 = 0;
        my $j = 0;
        # On parcourt l array des proteines du groupe courrant
        # $test = 0 si tous les BRHs de $prot possèdent les meme BRHs que $prot
        while (($test1 == 0) && ($j <= $#Brh_prot)) 
        {
			# On extrait un array des BRHs de chaque proteine, on construit la liste triée correspondante
            $tmpId = join('',sort(keys(%{$Brh2{$Brh_prot[$j]}}))); 
            $test1 = 1 if ($ref ne $tmpId); # $test1 vaut un si les listes de BRHs sont differentes
            $j ++; 
        }		
        # si toutes les listes de BRHs sont identiques, c'est un groupe
        if ($test1 == 0) 
        {
            # on enregistre le groupe dans un array
            # format table (1 cellule par genome, par default -)
			# intialisation: on construit un groupe vide (que des '-' pour tous les genomes)
            @TmpGroup = ();
            # On parcourt la liste des genomes
            for (my $j = 0; $j <= $#Geno; $j ++) 
            {
                $TmpGroup[$j] = '-';
            }
			# Pour chaque genome present dans le groupe, on remplace le tiret par l'identifiant de la sequence 
            foreach my $k (keys(%{$Brh2{$prot}})) 
            {    
                $TmpGroup[ $Brh2{$prot}{$k} ] = $k;	
                # $k est une proteine du groupe ($prot incluse)  Brh2{$prot}{$k} est le numero (la colonne) du genome de $k 
            }
			#on ecrit la ligne du fichier de sortie correspondant au groupe 
			# le nom des proteines du goupe ou - dans les cellules correspondant à leur genome (classement alphabetique des champs abrev) 
            print OUT join("\t",@TmpGroup)."\t".($#Brh_prot + 1)."\n";
        }		
		#Si il il y a des différences dans les paires de BRH (i.e. certaines seq ont des liens en plus ou en moins), on construit tous les groupes contenant la seq de depat $i
        else { #si brh differe
			@TmpGroup = (); # liste de tous les groupes contenant $i
			$nb_group_with_prot = 0; # nb groupe differents contenant $i
			
			#pour chaque sequence formant une paire de BRH avec $i
			for (my $l = 0; $l <= $#Brh_prot; $l ++) {
				$id1 = $Brh_prot[$l]; #sequence BRH à $i
				if ($id1 ne $prot) {
					# on recherche si la sequence $id1 possède des liens BRH avec toutes les sequences d'un groupe deja construit
					$test1 = 0; #vaut 1 si $id1 peut etre associe à au moins un groupe existant
					#pour tous les groupes deja construits avec $i
					for (my $k = 0; $k < $nb_group_with_prot; $k ++) {
						$test2 = 0; #valeur test : vaut 1 si manque un lien entre $id1 et le groupe $k
						$id2 = 0; #numero des identifiants du groupe $k
						@TmpId = keys(%{$TmpGroup[$k]}); #on recupere les id du groupe $k
						#tant que $id1 possede un lien de BRH avec les sequences du groupe $k; on continue
						while (($test2 == 0) && ($id2 <= $#TmpId )) {
							$test2 = 1 unless (exists($Brh2{$id1}{$TmpId[$id2]})); # arret si $id1 et $TmpId[$id2] ne forment pas une paire de BRH
							$id2 ++; #on prend la sequence suivante du groupe $k
						}
						# si la sequence $id1 possede des liens avec chaque sequence du groupe $k, on l'ajoute au groupe
						if ($test2 == 0) {
							$TmpGroup[$k]{$id1} = '';
							$test1 = 1; # pour signaler que la sequence a pu etre attribuee a au moins un groupe
						}
					}
					#si $id1 n'a pu etre attribuee a aucun groupe deja defini, on cree un nouveau groupe
					if ($test1 == 0) {
						%NewGroup = ();
						my $k = 0;
						#au depat, juste la sequence $id1
						$NewGroup{'1'}{$id1} = '';
						#on regarde dans les ID deja traites si certains ne pourraient pas etre ajoutes au groupe
						while ($k < $l) {
							#si la seq de depart ne correspond pas à la seq de depart $i (que l'on integrera dans chaque groupe a la fin
							if ($Brh_prot[$k] ne $prot) {
								#on recommence les test avec les nouveaux groupes
								$test1 = 0;
								#POur chaque nouveau groupe, on regarde si on ne peut pas ajouter des ID
								foreach $num_group (keys(%NewGroup)) {
									$test2 = 0; # vaut 1 si $Brh_prot[$k] a un lien manquant avec au moins une seq du nouveau groupe
									#pour toutes les sequences $id2 du groupe $num_group
									foreach $id2 (keys(%{$NewGroup{$num_group}})) {
										$test2 = 1 unless(exists($Brh2{$Brh_prot[$k]}{$id2})); # un si un lien BRH manquant avec $Brh_prot[$k]
									}
									#si tous les liens ont été trouvés, on peut ajouter l'ID $Brh_prot[$k] au groupe $num_group
									if ($test2 == 0) {
										#ajout de l'id
										$NewGroup{$num_group}{$Brh_prot[$k]} = '';
										$test1 = 1;
									}
								}
								# SI $Brh_prot[$k] associé a aucun groupe ($test1 == 0)
								# MAIS qu'il possede un lien BRH avec $id1 (exists($Brh2{$id1}{$Brh_prot[$k]})), 
								# ALORS on cree un nouveau groupe avec ces deux sequences
								if (($test1 == 0) && exists($Brh2{$id1}{$Brh_prot[$k]})) {
									$num_group = keys(%NewGroup) + 1;
									$NewGroup{$num_group}{$Brh_prot[$k]} = '';
									$NewGroup{$num_group}{$id1} = '';
									$k = 0; #on recommence à 0 pour voir si des sequences ne pourraient pas appartenir a ce groupe		
								}
								else { #si aucun nouveau groupe cree on continue avec la sequence suivante
									$k ++;
								}
							}
							# si on traite de la sequence de depart $i, on passe directement a la sequence suivante car elle appartient a tous les groupes
							else {
							  $k ++;
							}
						}
						# une fois que tous les groupes contenant $i et $id1 ont ete crees, on les enregistre dans la table globale
						foreach $num_group (keys(%NewGroup)) 
						{
							foreach $id2 (keys(%{$NewGroup{$num_group}})) {
								$TmpGroup[$nb_group_with_prot]{$id2} = '';
							}
							$nb_group_with_prot ++;
						}
					}		
				}
			}
			# une fois que tous les groupes d'orthologues contenant $i ont ete crees, 
			# on formate chacun des groupes
			for (my $k = 0; $k < $nb_group_with_prot; $k ++) 
			{
				$TmpGroup[$k]{$prot} = ''; #on ajoute notre ID de départ a chacun des groupes
				@TmpId = sort(keys(%{$TmpGroup[$k]}));
				unless (exists($AllGroup{join(',',@TmpId)}))
				{
					$AllGroup{join(',',@TmpId)} = '';
					@TmpId = (); #intialisation
					$size = 0;
					for (my $j = 0; $j <= $#Geno; $j ++) 
					{
						$TmpId[$j] = '-';
					}
					#chaque id est place dans le tableau en fonction de son genome
					foreach $id1 (keys(%{$TmpGroup[$k]})) 
					{
						$TmpId[$Brh2{$id1}{$id1}] = $id1; # INITIALISED VALUE HERE ($Brh2{$id1}{$id1})
					}
					for (my $j = 0; $j <= $#Geno; $j ++) 
					{
						#$size ++ if ($TmpId[$j] eq '-');
						$size ++ if ($TmpId[$j] ne '-');#modification cecile, la taille de certains groupe était fausse
						print OUT "$TmpId[$j]\t";
					}
					print OUT "$size\n";					
				}
			}		
        }
    }
}
close(OUT);
 
# Instruction sql
print "CREATE TABLE brh (";
for (my $j = 0; $j <= $#Geno; $j ++) {
	print "$Geno[$j] varchar(80), "; 
}
print "size integer)\n";
print "\\COPY brh FROM '$fileOut' WITH DELIMITER AS '\\t'\n";

