#!/usr/bin/perl -w
use strict;
#use DBI;

# Description
# Annotation des groupes d'orthologues

# perl 11_intersection.pl [dossier de travail] [type de BD, Pg ou mysql] [nom de la base de donnée] [utilisateur] > 11_intersection.sql
# perl 11_intersection.pl /home/cecile/git/fungipath/miseajourBD/local/intersection/ Pg cecile cecile > 11_intersection.sql

#######################################################################
# Declaration
###########
if ($#ARGV < 2) { print "perl 11_intersection_GO.pl <dossier git> <dossier des genomes> <dossier des groupes (dans mario)> [tag]"; exit();}
# Parametres
my $dirWork=$ARGV[0]; # repertoire GIT
my $dirFA=$ARGV[1]; # repertoire des genomes fasta
my $dirGroupsSeqs=$ARGV[2]; # sortie de mario (groupes au format fasta)
my $tag = $ARGV[3]; # tag du run
if ( $tag) { $tag = '_'.$tag;} else {$tag = "";}
my $dirAnnotation = $dirWork.'/annotations/annot_files/'; # INPUT : repertoire contenant les informations necessaires a l'annotation
my $dirComparison= $dirWork.'/annotations/group_GO'.$tag.'/';
my $file_annotated_seq_FUNGIpath = $dirAnnotation.'annotated_sequences_fungipath.tab'; # INPUT : Fichier contenant les sequences similaires entre FUNGIpath et SwissProt/MetaCyc
my $fileInAnnotation =  $dirAnnotation.'SwissProt_MetaCyc_go.sql';  # INPUT : fichier des sequences annotees dans SwissProt/MetaCyc
my $fileInAnnotationfa = $dirAnnotation.'SwissProt_MetaCyc_go.fa'; # INPUT : fichier des sequences annotees dans SP/MC au format fasta
my $outannotation = $dirWork.'/annotations/GO_annotated_group'.$tag.'/';

if (-e $dirWork.'/annotations/GO_annotated_group'.$tag.'/') {} else { system("mkdir $dirWork/annotations/GO_annotated_group$tag/");}

my $fileOutNoGroupEc =  $outannotation.'NoGroupGO.sql'; # OUTPUT : sequences depourvues d'orthologues mais annotees dans SwissProt/MetaCyc 
my $fileOutGroupEc =  $outannotation.'GroupGO.sql';  # OUTPUT : annotation des groupes d'orthologues 
my $fileOutId =  $outannotation.'ID_GO_Predit.sql';  # OUTPUT : annotation de chaque sequence presente dans un groupe

# cutoff
my $seuilBestEval = 1e-45; # 1er seuil sur la Evalue
my $seuilMinEval = 1e-20; # 2nd seuil sur la Evalue

# Autres variables
my @go=();
my $file = '';
my $nbSeq = 0;
my $ec = '';
my $idDb = '';
my $idSp = '';
my $idGroup = '';
my $req = '';
my $Res = ();
my $bestEval = 100;
my $i = 0;
my $j = 0;
my $bestScore = 0;
my %IdEc = ();
my %EcGroupTmp = ();
my %AnnotatedSequences_FUNGIpath = ();
my %ScoreTmp = ();
my %EvalTmp = ();
my %EcTmp = ();
my @Geno = ();
my @Id = ();
my @Res = ();
my @EcTmp = ();	
my @IdSp = ();
my @File = ();
my $vu=0;
my $idSppos='';
my %IdSP_seq=();
my $idtmp="";
my $seqtmp='';
my @tmp=();
#######################################################################
# Programme
###########
#recuperation des seq associees aux id de sp
open(IN,$fileInAnnotationfa)or die("no $fileInAnnotationfa");
while(<IN>){
  chomp;
  if(/^>/){
    if($vu==1){
      $IdSP_seq{$idtmp}=$seqtmp;
      $idtmp="";
      $seqtmp="";
    }
    else{
      $vu=1;
    }
    $_=~s/>//;
    @tmp=split(/\|/,$_);
    $idtmp=$tmp[0];
  }
  else{
    $seqtmp.=$_;
  }
}
$IdSP_seq{$idtmp}=$seqtmp;
close(IN);

# Recuperation des genomes presents dans FUNGIpath
#$req = $db->prepare('SELECT abrev FROM genome');
#$req->execute();
#while($Res = $req->fetchrow_hashref) {
#    push(@Geno,$$Res {'abrev'});
#}

chdir($dirFA);
my @tmpfa=glob("*");
foreach my $tf(@tmpfa){
        my $nomtmp=$tf;
        $nomtmp=~s/\.fa$//;
        push(@Geno,$nomtmp);
}

# Recuperation des annotations de chaque sequence dans SwissProt et/ou MetaCyc
if (-e $fileInAnnotation){
	open IN,$fileInAnnotation;
	while(<IN>) {
		#FABH_MYCTA      2.3.1.180	GOterm
		chomp($_);
		@Res = split(/\t/,$_);
		#print "\n $_ \n".$Res[2]."\n";
		#die;
		@EcTmp = split(/;/,$Res[2]);	
		foreach $ec (@EcTmp) {
			$ec=~s/ //g;
			$IdEc {'swissprot_metacyc'} {$Res[0]} {$ec} = '';
		}
	}
	close(IN);
}
else{
	die ("Error : No file $fileInAnnotation\n");
}

# Recuperation des sequences de FUNGIpath annotatees dans SwissProt et/ou MetaCyc
if (-e $file_annotated_seq_FUNGIpath) {
    open IN,$file_annotated_seq_FUNGIpath;
    while(<IN>) {
        #RIFK_YARLI      YALI0B01826g ec
        chomp($_);
        @Res = split(/\t/,$_);
        @EcTmp = split(/;/,$Res[3]);
        foreach $ec (@EcTmp){
            $AnnotatedSequences_FUNGIpath{$Res[1]} {$ec} = '';
            $IdEc {'db'} {$Res[1]} {$ec} = '';
        }
    }
    close(IN);
}
else {
    die ("Error : file $file_annotated_seq_FUNGIpath not found\n");
}

# Creation des fichiers OUTPUT
open OUTGROUPEC,'>'.$fileOutGroupEc;
open OUTIDEC,'>'.$fileOutId;

# Recuperation des fichiers results
chdir($dirComparison);
@File = glob("GO*hmmsearch");

# Pour chaque fichier résultat
foreach $file (@File) {#adecommenter
#$file="3.hmmsearch";#aretirer
	print "$file\n";
	if ($file =~ /^GO(.*).hmmsearch$/) { 
		$idGroup = $1 ;
		#initialisation
		$nbSeq = 0;
		@Id = ();

		# Récupération des IDs des séquences du groupe $idGroup
# 		$req = $db->prepare("SELECT * FROM ortho4meth WHERE nb = $idGroup");
# 		$req->execute();
# 		while($Res = $req->fetchrow_hashref) {
# 			foreach $i (@Geno) {
# 				 push(@Id,$$Res {$i}) if ($$Res {$i} ne '-');
# 			}
# 		}
		open(IN,"$dirGroupsSeqs/$idGroup.fa")or die("impossible d'ouvrir le fichier seq du groupe");
                while(<IN>){
                        if(/>(.+)\t.+/){
                                push(@Id,$1);
                        }
                }
                close(IN);


#		$req = $db->prepare("SELECT id0 FROM maapf_v1_201401_1e5_a20_i3 WHERE nb = '$idGroup'");
#		$req->execute();
#		while($Res = $req->fetchrow_hashref) {
#		  push(@Id,$$Res {'id0'});
#		}
		
# 		$req = $db->prepare("SELECT inpara FROM inpara4meth WHERE nb = $idGroup");
# 		$req->execute();
# 		while($Res = $req->fetchrow_hashref) {
# 			push(@Id,$$Res {'inpara'});
# 		}
		$nbSeq = $#Id + 1;
		
		if ($nbSeq != 0){
			# Recuperation des eventuelles annotations des sequences du groupe deja annotees initialement dans SwissProt/MetaCyc
			%EcGroupTmp = (); # Liste des annotations attribuees a au moins une sequence du groupe
			foreach $idDb (@Id) {
				if (exists($IdEc {'db'} {$idDb})) {
					# On supprimer au fur et a mesure les sequences annotees possedant des orthologues
					# pour conserver ensuite les annotations des sequences depourvues d'orthologues
					delete($AnnotatedSequences_FUNGIpath{$idDb});
					# Recuperation des annotations de la sequence $idDb
					foreach $ec (keys(% {$IdEc {'db'} {$idDb}})) {
						$EcGroupTmp {$ec} = '';
					}
				}
			}
			$idSp = '';
			$bestEval = 10;
			@IdSp = ();		
			%ScoreTmp = ();
			%EvalTmp = ();
			
			# Ouverture et lecture du fichier résultat
			open IN,$dirComparison.$file;
			while(<IN>) {
				chomp($_);
				# Recherche des lignes correspondant au meilleur résultat 
				# A moins que la ligne soit un commentaire
				unless ($_ =~ /^#/) {
					$_ =~ s/ +/ /g;
					# Recuperation des colonnes dans un tableau
					@Res = split(/ /,$_);
					# Si la E-value est inferieure ou egal a la meilleure Eval
					if ($bestEval >= $Res[4]) {
						# Recuperation de l'ID
						$Res[0] = $1 if ($Res[0] =~ /^(.*?)\|.*/);
						# Si la E-value est egale a la meilleure Eval
						# C'est que plusieurs resultats avec la meme eval ont ete trouves
						if ($bestEval == $Res[4]) {
							# Ajout de la sequence dans le tableau des IDs trouves
							push(@IdSp,$Res[0]);
						}
						# Sinon, c'est le premier resultat
						else{
							# Creation du tableau des IDs trouves
							@IdSp = ($Res[0]);
							$bestEval = $Res[4]; # Enregistrement de la meilleure Evalue
							# Reinitialisation des tableaux
							%ScoreTmp = ();
							%EvalTmp = ();
						}
						$ScoreTmp{$Res[0]} = $Res[8]; # Enregistrement du score du meilleur domaine
						$EvalTmp{$Res[7]}{$Res[0]} = ''; # Enregistrement de la Evalue du meilleur domaine
					}
					# Si resultat ne correspond pas au meilleur resultat, on s'arrete
					else {
						last;
					}
				}
			}
			close(IN);
			
			# Si la meilleure E-vale ne depasse pas le seuil le moins stringent, 
			# le groupe n'est pas annote et on ne va pas plus loin
			
			if ($bestEval < $seuilMinEval){ 
				# Si plusieurs sequences ont ete trouvees aves la meme E-value globale,
				# on prend celle avec la meilleure eval sur le meilleur domaine
				if ($#IdSp > 0){#si on a plusieurs sequences avec la meilleure evalue
					@Res = keys(%EvalTmp);
					@Res = sort { $a <=> $b } @Res;
					$i = shift(@Res);#on prend la meilleur evalue
					# Si plus d'une seq possede la meme Evalue,
					# on prend celle avec le meilleur score sur le meilleur domaine
					if (keys(%{$EvalTmp{$i}}) > 1){#si plusieurs seq avec la meilleur evalue
						@Res = ();#les scores des seq avec les meilleurs evalues
						foreach $idSp (keys(%{$EvalTmp{$i}})){
							push(@Res,$ScoreTmp{$idSp});
						}
						@Res = sort { $a <=> $b } @Res;#trie des scores
						$j = pop(@Res);#on garde le plus grand score
						@IdSp = ();
						foreach $idSp (keys(%{$EvalTmp{$i}})){
							push(@IdSp,$idSp) if ($j == $ScoreTmp{$idSp});#on garde meilleur score et meilleure sequence
						}
					}
					# Sinon, recuperation de l'ID dans le tableau
					else{
						@IdSp = keys(%{$EvalTmp{$i}});
					}
				}
				@EcTmp = ();
				%EcTmp = ();

				# Si la Evalue eval est inferieure au seuil le plus stringent
				# Le groupe est annote directement
				if ($bestEval < $seuilBestEval) {
					# Récupération des ECs correspondant au meilleur resultat
					if($#IdSp>0){
						print "attention, plusieurs id dans IdSp\n";
					}
					foreach $idSp (@IdSp){
						#print $idSp."\n";
						foreach $ec (keys(% {$IdEc {'swissprot_metacyc'} {$idSp}})){
							$EcTmp{$ec}{$idSp} = '';
						}
					}
					# Liste des annotations enzymatiques a transferer au groupe $idGroup
					@EcTmp = keys(%EcTmp);
					if ($#EcTmp >= 0) {
						# On enregsitre les annotations attribuees au groupe
						foreach $ec (@EcTmp) {
							#print "l'ec ligne 283 $ec\n";
							@go=split(/,/,$ec);
							my $nbtest=1;
							foreach $idSppos (keys(%{$EcTmp{$ec}})){
								#print "nbtest = $nbtest\n";
								$nbtest++;
								print OUTGROUPEC "$idGroup\t$go[0]\t$go[1]\t$idSppos\t$IdSP_seq{$idSppos}\n";
							}
						}
						# On enregsitre les annotations attribuees a chaque sequence du groupe
						# On regarde en plus si la sequence etait initiallement annotee au non
						# Pour chaque sequence du groupe
						foreach $idDb (@Id) {
							# Pour chaque annotation du groupe
							foreach $ec (@EcTmp){
								# Si la sequence $idDb possedait une annotation initialement dans SwissProt / MetaCyc
								@go=split(/,/,$ec);
								if (exists($IdEc {'db'} {$idDb})) {
									# Si la sequence $idDb possede l'EC $ec egalement dans SwissProt / MetaCyc
									# on l'indique : SP = SwissProt / MetaCyc et PR = predit
									if (exists($IdEc {'db'} {$idDb} {$ec})) {
										print OUTIDEC "$idDb\t$go[0]\t$go[1]\tSP//PR\n";
									}
									# Sinon on indique que cette annotation est specifique a FUNGIpath (PR = predit)
									else {
										print OUTIDEC "$idDb\t$go[0]\t$go[1]\tPR\n";
									}
								}
								# Sinon on indique que cette annotation est specifique a FUNGIpath (PR = predit)
								else {
									print OUTIDEC "$idDb\t$go[0]\t$go[1]\tPR\n";
								}
							}
							# Si la sequence $idDb possede des annotations presentes SwissProt / MetaCyc mais non trouvees dans FUNGIpath
							# On l'enregistre
							if (exists($IdEc {'db'} {$idDb})) {
								foreach $ec (keys(% {$IdEc {'db'} {$idDb}})) {
									unless($ec eq "-"){
										@go=split(/,/,$ec);
										#print STDERR "idDb=$idDb\n";
										#print STDERR "go0 $go[0]\n";
										#print STDERR "go1 $go[1]\n";
										print OUTIDEC "$idDb\t$go[0]\t$go[1]\tSP\n" unless (exists($EcTmp {$ec}));
									}
								}
							}
						}
					}
					else {
						die("Error 1 in $file : no EC for @IdSp, eval= $bestEval\n");
					}
				}
				else {  
					# Si la meilleure Evalue est comprise entre nos deux seuils
					# On regarde si une sequence du groupe ne possederait pas la meme annotation
					# Liste des annotations enzymatiques attribuees aux sequences du groupe
					@Res = keys(%EcGroupTmp);
					if ($#Res >= 0) {
						# On regarde si l'une des annotations ne serait pas attribuee a l'une des sequences du groupe
						foreach $idSp (@IdSp){ # Pour chaque sequence utilisee pour l'annoation
							# Pour chacune de ces annotations enzymatiques
							foreach $ec (keys(% {$IdEc {'swissprot_metacyc'} {$idSp}})){
								# On enregistre l'EC  si il a ete attribue a une sequence du groupe
								$EcTmp{$ec}{$idSp} = '' if (exists($EcGroupTmp {$ec}));
							}
						}
						# Liste des annotations enzymatiques retenues pour l'annotation
						@EcTmp = keys(%EcTmp);
						# Si au moins un resultat a ete trouve, on enregistre le resultat
						if ($#EcTmp >= 0) {
							# On enregsitre les annotations attribuees au groupe
							foreach $ec (@EcTmp) {
								@go=split(/,/,$ec);
								foreach $idSppos(keys(%{$EcTmp{$ec}})){
									print OUTGROUPEC "$idGroup\t$go[0]\t$go[1]\t$idSppos\t$IdSP_seq{$idSppos}\n";
								}
							}
							# On enregsitre les annotations attribuees a chaque sequence du groupe
							foreach $idDb (@Id) {
								foreach $ec (@EcTmp){
									@go=split(/,/,$ec);
									if (exists($IdEc {'db'} {$idDb})) {
											if (exists($IdEc {'db'} {$idDb} {$ec})) {
												print OUTIDEC "$idDb\t$go[0]\t$go[1]\tSP//PR\n";
											}   
											else {
												print OUTIDEC "$idDb\t$go[0]\t$go[1]\tPR\n";
											}
									}
									else {
											print OUTIDEC "$idDb\t$go[0]\t$go[1]\tPR\n";
									}
								}
								#enregistrement des ec non retrouvés par notre annot
								foreach $ec (keys(% {$IdEc {'db'} {$idDb}})) {
									@go=split(/,/,$ec);
									print OUTIDEC "$idDb\t$go[0]\t$go[1]\tSP\n" unless (exists($EcTmp {$ec}));
								}
							}
						}
					}
				}
			}
		}
	}
	else{
		die ("Error : Problem with $file\n");
	}
}#adecommenter
close(OUTGROUPEC);
close(OUTIDEC);

# On enregistre les sequences depourvues d'orthologues mais annotees dans MEtaCYc ou SwissProt
open OUTNOGROUPEC,'>'.$fileOutNoGroupEc;
foreach $idDb (keys(%AnnotatedSequences_FUNGIpath)){
    foreach $ec (keys(%{$AnnotatedSequences_FUNGIpath{$idDb}})){
	unless($ec eq "-"){
		@go=split(/,/,$ec);
        	print OUTNOGROUPEC "$idDb\t$go[0]\t$go[1]\n";
    
	}
    }
}
close(OUTNOGROUPEC);

print "CREATE TABLE idgopredit_maapf_v1_201401_1e5_a20_i3 (id0 varchar(80) , go varchar(20),evidencecode varchar(80) ,predit varchar(12), CONSTRAINT kidgopr_maapf_v1_201401_1e5_a20_i3_pkey PRIMARY KEY (id0,go,evidencecode));\n";
print "CREATE TABLE groupgo_maapf_v1_201401_1e5_a20_i3 (nb integer REFERENCES idgroup_sequentiel_v1 (nb), go varchar(20),evidencecode varchar(80),idannotSPMC varchar(80),seqannotSPMC text,  CONSTRAINT groupgo_maapf_v1_201401_1e5_a20_i3_pkey PRIMARY KEY (nb,go,evidencecode,idannotSPMC));\n";
print "CREATE TABLE nogroupgo_maapf_v1_201401_1e5_a20_i3 (id0 varchar(80), go varchar(20),evidencecode varchar(80), CONSTRAINT nogroupgo_maapf_v1_201401_1e5_a20_i3_pkey PRIMARY KEY (id0,go,evidencecode));\n";
print "\\COPY idgopredit_maapf_v1_201401_1e5_a20_i3 FROM '$fileOutId' WITH DELIMITER AS E'\\t'\n";
print "\\COPY groupgo_maapf_v1_201401_1e5_a20_i3 FROM '$fileOutGroupEc' WITH DELIMITER AS E'\\t'\n";
print "\\COPY nogroupgo_maapf_v1_201401_1e5_a20_i3 FROM '$fileOutNoGroupEc' WITH DELIMITER AS E'\\t'\n";

