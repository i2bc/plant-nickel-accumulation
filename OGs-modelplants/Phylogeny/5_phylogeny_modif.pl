#!/usr/bin/perl -w
use strict;

#script qui cree les groupes d'orthologues a partir des resultats issus de l'analyse des arbres
#modif par cecile le 18/09/2013
#perl 5_phylogeny.pl /home/cecile/Bureau/miseajourBD/local/ > 5_phylogeny.sql
if($#ARGV < 1) 
{ 
	print "perl my5_phylogeny_modif.pl <repertoire de travail> <repertoire des fasta> \n";
	exit(); 
}
my $workDir = $ARGV[0];
my $dirGeno=$ARGV[1];
my $outd = $workDir.'/Phylogeny';
my $dirTreeAnalysis = $outd.'/STEP4-2/'; # INPUT : repertoire contenant les paires d'ortho
print $dirTreeAnalysis."\n";

my $dirOut = $outd.'/STEP5/'; # OUTPUT : repertoire ou sera enregistre les groupes d'orthologues au format SQL
print $dirOut."\n";
# Fichiers 
my $fileOutOrtho = $dirOut.'phylogeny_orthologous_groups.sql'; # OUTPUT : fichiers contenant les groupes d'orthologues
my $fileOutInpara = $dirOut.'phylogeny_inparalogs.sql'; # OUTPUT : fichiers contenant les inparalogues

# Autres variables
my $file = '';
my $id1 = '';
my $id2 = '';
my $idRef = '';
my $i = 0;
my $j = 0;
my $k = 0;
my $geno1 = '';
my $geno2 = '';
my $geno = '';
my $test = 0;
my $nbGroup = 0;
my $maxId = 0;
my $maxGeno = 0;

my %PairOrtho = (); 
my %AlreadySave = ();
my %IdGenome = ();
my %IdOld = ();
my %IdNew = ();
my %PairOk = ();
my %TmpGroup = ();

my @Specie = (); 
my @Line = ();
my @Ortho1 = ();
my @Ortho2 = ();
my @Tmp = ();
my @Id = ();

###########
# Programme
###########

# Creation du repertoire OUPUT si il n'existe pas
mkdir ($dirOut) unless (-e $dirOut);

# Verification que les fichiers d'entree existent bien
die ("No tree directory $dirTreeAnalysis\n") unless (-e $dirTreeAnalysis);
die ("No input directory $dirGeno\n") unless (-e $dirGeno);

# Recuperation du nom des espèces et tri par ordre alphabetique
if (-e $dirGeno){
	chdir($dirGeno);
	@Tmp = glob("*.fa");
	foreach $file (@Tmp) {
		$file =~ s/.fa$//;
		push(@Specie,$file);
	}
}
else{
	die ($dirGeno." not found.\n");
}
@Specie = sort(@Specie);

#Lecture des paires d'ortho par la methode de phylogenie et construction des groupes d'orthologues
##################################################################################################

# Ouverture des fichiers ou serons enregistres les groupes
open OUTORTHO, ">".$fileOutOrtho;
open OUTINPARA, ">".$fileOutInpara;


# Recuperation des donnees et construction des groupes
chdir ($dirTreeAnalysis);
@Tmp = glob("arbre*");
# lecture des fichiers plats issus de l'analyse des arbres ( contenu du dossier STEP4-2_<tag>)
foreach $file (@Tmp)
{
  #Initialisation des variables pour la nouvelle famille proteique
  %PairOrtho = ();
  %AlreadySave = ();
  %IdGenome = ();

  #Lecture des paires d'orthologues pour la famille proteique $file
  print "$file...\n";
  open IN, $dirTreeAnalysis.$file;
  while(<IN>) 
  {
    @Line = split(/\t/,$_); #Recuperation du couple orthologues
    if ($#Line >= 2 ) 
    {
        $id1 = ''; $id2 = ''; $geno1 = ''; $geno2 = '';
		if($Line[1]=~/^((.+?)_.+)$/)
		{
			$id1 = $1;
			$geno1 = $2;		
			$IdGenome{$id1} = $geno1;
		}
		if($Line[2]=~/((.+?)_.+)$/)
		{            
            $id2 = $1;
            $geno2 = $2;
            $IdGenome{$id2} = $geno2;
		}
		if (($id1 ne '') && ($id2 ne '')) 
		{
	    	$PairOrtho{$id1}{$id2} = '';
	    	$PairOrtho{$id2}{$id1} = '';
		}
		else 
		{
			die("Error 1 in $i with $_");
		}
    } 	
    else 
    {
		print("Error 2 in $i with $_");
    }
  }
  close(IN);

  # Construction des groupes d'orthologues associes a cette famille
  #################################################################
  # Pour chaque sequence de la famille
  foreach $idRef (keys(%PairOrtho))
  {
      # A moins que la sequence est deja ete analysee car une autre sequence de son groupe a ete analysee
      unless (exists($AlreadySave{$idRef}))
      {
		%PairOk = ();
		%IdOld = ();
		%IdNew = ();

		#On recupere tous les orthos à $idRef,et tous les orthos à ces orthos ... (car peuvent etres differents)
		$IdNew{$idRef} = ''; # au depart seulement $idRef
		# Tant qu'un nouvel ortho a ete trouve, on recupere tous ces orthos
		while (keys(%IdNew) != 0) 
		{
			#construction d'un array de proteines orthologues a idRef
            @Id = keys(%IdNew); # Nouveaux orthologues
			# Chaque nouvel orthologue,devient un ancien
            foreach $id1 (@Id) 
            {
                $IdOld{$id1} = ''; #  IdOld stocke les proteines deja vues
				$AlreadySave{$id1} = 'ok'; 
            }
            %IdNew = ();
			# pour chaque orthologue precedemment ajoute, on regarde si il ne possede pas un orthologue non deja sauve
            foreach $id1 (@Id) 
            {
				# pour chaque orthologue à $id1 ($id2)
                foreach $id2 (keys(%{$PairOrtho{$id1}})) 
                {
					#on l'enregistre si non deja present
                    $IdNew{$id2} = '' unless (exists($IdOld{$id2}));
                }
            }
		}
		@Id = keys (%IdOld);

	  #pour chaque paire de sequence, 
	  for($i = 0; $i < $#Id; $i ++)
	  {
	    $id1 = $Id[$i]; #ID1
	    $geno1 = $IdGenome{$id1}; #genome de $id1
            
	    for($j = $i + 1; $j <= $#Id; $j ++)
	    {
	      $id2 = $Id[$j]; #ID2
	      $geno2 = $IdGenome{$id2}; #genome de $id2
	      #initialisation
	      @Ortho1 = ();
	      @Ortho2 = ();
	      $test = 0;
	      # Si les genomes sont identiques, on verifie que les seq sont des inparaloques et donc partagent exactement les memes ortholoques
	      if ($geno1 eq $geno2)
	      {
			@Ortho1 = sort(keys(%{$PairOrtho{$id1}}));	
			@Ortho2 = sort(keys(%{$PairOrtho{$id2}}));	
			# si meme nombre, on regarde si c'est bien les memes sequences
			if ($#Ortho1 == $#Ortho2)
			{
		      $k = 0;
		      while (($test == 0) && ($k <= $#Ortho1))
		      {
				$test = 1 if ($Ortho1[$k] ne $Ortho2[$k]);
				$k ++;
		      }
			}
	      }
	      #si les genomes sont differents,on verifie que les deux sequences sont orthologues et qu'elles partagent les memes orthologues avec les autres genomes
	      #ex: si $id2 a un inparalogue, $id1 aura deux orthologues avec le genome $id1
	      else{
		  @Ortho1 = keys(%{$PairOrtho{$id1}});
		  $k = 0;
		  while (($test == 0) && ($k <= $#Ortho1)){
		    #si l'ortho à $id1 n'appartient pas au meme genome que $id2 et qu'il n'a pas ete trouve comme orthoa $id2, on arrete
		    $test = 1 if (($IdGenome{$Ortho1[$k]} ne $geno2) && !(exists($PairOrtho{$id2}{$Ortho1[$k]})));
		    $k ++;
		  }
		  #si tous les orthos à $id1 ont ete trouves, on verifie que l'inverse est vrai
		  if ($test == 0){
		    @Ortho2 = keys(%{$PairOrtho{$id2}});
		    $k = 0;
		    while (($test == 0) && ($k <= $#Ortho2)){
		      #si l'ortho à $id2 n'appartient pas au meme genome que $id1 et qu'il n'a pas ete trouve comme orthoa $id1, on arrete
		      $test = 1 if (($IdGenome{$Ortho2[$k]} ne $geno1) && !(exists($PairOrtho{$id1}{$Ortho2[$k]})));
		      $k ++;
		    }
		  }
	      }
	      # SI les sequences sont des inparalogues ou appartiennent au meme groupe d'orthologues, on enregistre la paire
	      if ($test == 0){
			$PairOk{$id1}{$geno2}{$id2} = '';
			$PairOk{$id2}{$geno1}{$id1} = '';
	      }
	    }
	  }
	  # Construction de tous les groupes possibles, 
	  # On ne selectionnera ensuite que celui contenant le plus de genomes et le plus de sequences 
	  # car sinon on ferait dans l'exemple ci dessous 3 groupes d'orthologues ce qui pourrait fausser les groupes lors de l'intersection
	  # donc on fait le choic de prendre un groupe par les 3, celui contenant le plus de genomes et le plus de genes par mi A, B et C
	  #  ______________ Basidiomycota genes C
          # |    
          # |	 _________  Ascomycota genes B
          # |___|
          #     |_________  Ascomycota genes A
          #	
	  $maxId = 0; # Nombre de seq max dans un groupe
	  $maxGeno = 0; # Nombre de genomes max dans un groupe
	  %TmpGroup = (); # Groupe correspondant aux conditions optimales
	  foreach $id1 (keys(%PairOk)){
	      @Ortho2 = keys(%{$PairOk{$id1}}); #liste des genomes possedant un orthologue avec $id1
	      $geno1 = $IdGenome{$id1} ; #espèce de $id1
	      push(@Ortho2,$geno1) unless (exists($PairOk{$id1}{$geno1})); #on ajoute $geno1 a la liste d'espece si non deja ajoute - ce qui est la cas si $id1 possede un imparalogue
	      #Si groupe possedent le plus de genomes, on l'enregistre comme ref
	      if (($#Ortho2 + 1) > $maxGeno){
		  $maxGeno = $#Ortho2 + 1;
		  %TmpGroup = ();
		  $TmpGroup{$geno1}{$id1} = '';
		  $maxId = 1;
		  foreach $geno2 (@Ortho2){
		      foreach $id2 (keys(%{$PairOk{$id1}{$geno2}})){
			  $TmpGroup{$geno2}{$id2} ='';
			  $maxId ++;
		      }
		  }
	      }
	      #Si autant de genomes, on regare le nombre de sequences
	      elsif (($#Ortho2 + 1) == $maxGeno){
		  $i = 1;
		  foreach $geno2 (@Ortho2){
		      foreach $id2 (keys(%{$PairOk{$id1}{$geno2}})){
			  $i ++;
		      }
		  }
		  if ($i > $maxId){
		      $maxId = $i;
		      %TmpGroup = ();
		      $TmpGroup{$geno1}{$id1} = '';
		      foreach $geno2 (@Ortho2){
			foreach $id2 (keys(%{$PairOk{$id1}{$geno2}})){
			  $TmpGroup{$geno2}{$id2} ='';
			}
		      }
		  }
	      }
	      
	  }
	  #Enregistrement du groupe
	  $nbGroup ++;
        for ($i = 0; $i <= $#Specie; $i ++) {
            if (exists($TmpGroup{$Specie[$i]})) {
				@Ortho1 = keys(%{$TmpGroup{$Specie[$i]}});
                print OUTORTHO $Ortho1[0]."\t";
				#si plusieurs sequences pour meme genome, les autres sont enregistrees comme inpara
				for ($j = 1; $j <= $#Ortho1; $j ++){
					print OUTINPARA $Ortho1[$j]."\t".$nbGroup."\n";
				}
            } 
            else {
                print OUTORTHO "-\t";
            }
          }
          print OUTORTHO $nbGroup."\n";
      }
  }
}
close (OUTORTHO);
close (OUTINPARA);

print $nbGroup." orthologous groups.\n";
print "Instruction pour la base de données : \n";
print "CREATE TABLE phyloortho (";
#a chaque genome on associe un numero qui lui est specifique
for ($i = 0; $i <= $#Specie; $i ++) {
    print "$Specie[$i] character varying(80), ";
}
print "size integer PRIMARY key)\n";
print "CREATE TABLE phyloinpara (inpara character varying(80), nb integer  REFERENCES phyloortho (nb), CONSTRAINT phyloinpara_pkey PRIMARY KEY (inpara,nb));\n";
print "\\COPY phyloortho FROM '".$fileOutOrtho."' WITH DELIMITER AS '\\t'\n";
print "\\COPY phyloinpara FROM '".$fileOutInpara."' WITH DELIMITER AS '\\t'\n";
