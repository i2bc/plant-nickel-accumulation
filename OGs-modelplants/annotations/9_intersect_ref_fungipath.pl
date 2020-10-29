#!/usr/bin/perl -w
use strict;
#use DBI;
use LWP::Simple;
use Net::FTP;
#dernieres modifs Cecile Pereira 14/06/2015, modifications mineures CD

# Description
# Extraction des sequences annotees dans MetaCyc et/ou SwissProt et presentes dans FUNGIpath
# needle doit etre installé
# perl 9_intersection_ECetGO.pl [repertoire de travail] [type de base de donnée mysql ou Pg] [nom de la base de donnée] [nom de l'utilisateur de la base de donnee]> 9_intersection.out
#perl 9_intersection.pl /home/cecile/git/fungipath/miseajourBD/local/intersection/ Pg cecile cecile > 9_intersection.out

#dans le fichier resultats ligne du type:
#id dans swissprot ou metacyc \t 
#id dans Fungipath \t 
#	Identical (si la sequence a été annotée car elle est identique a une sequence de swissprot ou de metacyc) 
#	Similar (si une seul seq a été trouvée pertinente par blast) 
#	Similar_PlsID (si annotation par similarité et pls sequences ont été trouvés pertinente, on garde alors celle avec la plus basse evalue) \t 
#0 si la seq est de metacyc 1 pour swissprot
if($#ARGV < 2) 
{ 
	print "9_intersect_ref_fungipath.pl <git directory> <liste des genomes> <fasta files directory> [tag] "; 
	exit();
}

# Parametres
my $dirWork = $ARGV[0]; # INPUT : Repertoire de travail
my $fiabrev_nom=$ARGV[1];# INPUT : fichier avec abrev \t nom de l'espece \t taxid
my $dosfa=$ARGV[2]; # fasta input files
my $tag = $ARGV[3]; # INPUT : repertoire contenant les fichiers necessaires a l'annotation des groupes
my $dirIn = $dirWork.'/annotations/annot_files';
if ($tag) 
{
	$dirIn  .= '_'.$tag.'/';	
}
else
{
	$dirIn  .= '/';
}	

my @FileInAnnot = ($dirIn.'SwissProt_ec.sql',$dirIn.'MetaCyc_ec.sql',$dirIn.'SwissProt_go.sql'); # INPUT :  fichier contenant les sequences annotees
my $url_ftp_taxonomy = 'ftp.ncbi.nih.gov'; # URL pour trouver le fichier a telecharger contenant la liste des especes
my $url_dir_taxonomy = 'pub/taxonomy/';  # repertoire
my $url_file_taxonomy = 'taxdump.tar.gz';  # Nom de l'archive
my $file_taxonomy = 'names.dmp';  # Nom du fichier qui nous interesse dans l'archive
my $dirTmp = $dirIn.'tmp2/'; # OUTPUT : repertoire temporaire contenant les BLAST realises (car sequences de notre base peuvent differer de celles de SwissProt et MetaCYc si version differente utilisees
my $file_annotated_seq_FUNGIpath = $dirIn.'annotated_sequences_fungipath.tab'; # OUTPUT : Fichier contenant les sequences similaires entre FUNGIpath et SwissProt/MetaCyc

# Autres variables
my $ftp = '';
my $req = '';
my $Res = ();
my $id = '';
my $name = '';
my $geno = '';
my $file = '';
my $i = '';
my $perc_identity = 0;
my $perc_align = 0;
my $test = 0;
my @Line = ();
my @File = ();
my %IDtaxon_FUNGIpath = ();
my %IdName = ();
my %Name = ();
my %Specie = ();
my %SpecieSeq = ();
my %Id_Blast = ();
my %IdEc = ();
my %ID_equivalent_identical = ();
my %ID_equivalent_similar = ();
my $basederef=0;
my $nomaecrire="";

#######################################################################
# Programme
###########

# Suppression des fichiers OUTPUT si deja crees
unlink($file_annotated_seq_FUNGIpath) if (-e $file_annotated_seq_FUNGIpath);

# Suppression du repertoire temporaire si existe deja
system('rm -rf '.$dirTmp) if (-e $dirTmp);
mkdir($dirTmp);

# lecture des codes, noms et Taxon IDs des genomes
open(IN,$fiabrev_nom) or die("impossible d'ouvrir le fichier de correspondance des noms");
while(<IN>)
{
	print $_;
	chomp;
	my @tmpn=split(/\t/,$_);
	$IDtaxon_FUNGIpath{$tmpn[2]}{'abrev'} = $tmpn[0];
	$IDtaxon_FUNGIpath{$tmpn[2]}{'complet'} = $tmpn[1];
}
close(IN);


chdir($dosfa);
my @tmpfa=glob("*.fa");
my %idseq=();
foreach my $f(@tmpfa)
{
	open(IN,$f)or die("impossible d'ouvrir le fichier fasta $f");
	my $genome0tmp=$f;
	$genome0tmp=~s/\.fa$//;
	my $id0tmp='';
	my $seq0tmp='';
	while(<IN>)
	{
		chomp;
		if(/>(.+)/)
		{
			my @tmpid=split(/\s/,$1);
			unless($id0tmp eq '')
			{
				#print "$genome0tmp $seq0tmp $id0tmp\n";
				$SpecieSeq{$genome0tmp}{$seq0tmp}=$id0tmp;
				$idseq{$id0tmp}=$seq0tmp;
			}
			$seq0tmp='';
			$id0tmp=$tmpid[0];
		}
		else
		{
			$seq0tmp.=$_;
		}
	}
	#derniere sequence
	$SpecieSeq{$genome0tmp}{$seq0tmp}=$id0tmp;
	$idseq{$id0tmp}=$seq0tmp;
	close(IN);
}

# Telechargement de la base de donnees TAXONOMY
unless (-e $dirIn.$file_taxonomy)
{
	chdir($dirIn);
	my($filetype) = "binary"; 
	$ftp = Net::FTP->new($url_ftp_taxonomy,  Timeout => 30, debug => 1) or die "Cannot connect to some.host.name: $@";
	$ftp->login("anonymous",'-anonymous@') or die "Cannot login ", $ftp->message;
	$ftp->$filetype;
	$ftp->cwd($url_dir_taxonomy) or die "Cannot change working directory ", $ftp->message;
	$ftp->get($url_file_taxonomy) or die "get failed ", $ftp->message ;
	$ftp->quit;

	# Desarchivage et extraction du fichier qui nous interesse
	if ($url_file_taxonomy =~ /gz$/)
	{
		system("gunzip $url_file_taxonomy");
		$url_file_taxonomy =~ s/\.gz$//;
		system ("tar -xvf $url_file_taxonomy $file_taxonomy");
		if (-e $dirIn.$file_taxonomy)
		{
			unlink($dirIn.$url_file_taxonomy);
		}
		else
		{
			die ("Error : No file ".$dirIn."$file_taxonomy\n");
		}
	}
}
# Lecture du fichier donnant les noms asscoies achque taxonomy ID
# 85982   |       Serpula lacrymans       |               |       scientific name |
# 85982   |       dry rot fungus  |               |       common name     |

open IN,$dirIn.$file_taxonomy;
while (<IN>){
    chomp($_);
    if ($_ =~ /^([0-9]*)[\s]*\|[\s]*([\S].*[\S])[\s]*\|.*\|(.*)\|$/){
        $id = $1;
        $name = $2;
	$i = $3;
	unless ($i =~ /authority/){
	# si meme nom pour plusieurs especes on supprime ce nom
		if (exists($Name{$name}) && ($Name{$name} != $id)){
			if (exists($IdName{$Name{$name}})){
				delete($IdName{$Name{$name}}{$name}) if (exists($IdName{$Name{$name}}{$name}));
				delete($IdName{$Name{$name}}) if (keys(%{$IdName{$Name{$name}}}) == 0);
			}
			delete($Specie{$name}) if (exists($Specie{$name}));
		}
		# Sinon on l'enregistre
		elsif (exists($IDtaxon_FUNGIpath{$id})) {
			$IdName{$id}{$name} = '' ; # Pour chaque ID taxon, liste de tous les noms associes
			$Specie{$name} = $id; # A chaque nom d'espece on associe son toxonomy ID
			$Name{$name} = $id; # Pour connaitre toute la liste des noms d'especes (et savoir si un nom est utilise plusieurs fois
		}
	}
    }
    else{
		print ("Warning : $_\n");
    }

}
close(IN);

# Pour chaque taxonomy ID present dans FUNGIpath
foreach $id (keys(%IDtaxon_FUNGIpath)) {
	# Si son taxonomy ID a ete trouve dans swissprot, on affiche la liste des noms qui lui sont associes
    if (exists($IdName{$id})){
        print "$id\t$IDtaxon_FUNGIpath{$id}{'complet'}\t".join(',',keys(%{$IdName{$id}}))."\n";
    }
	# Sinon on signale que l'on a rien trouve
    else{
        print "Warning $id\t$IDtaxon_FUNGIpath{$id}{'complet'}\tnot found\n";
    }
}

# Lecture des fichiers de donnees (Sequences annotees dans Swiss-Prot et/ou MetaCyc)
foreach $file (@FileInAnnot){#pour les fichiers SwissProt_ec.sql et MetaCyc_ec.sql et SwissProt_go.sql
    open IN,$file or die("impossible de lire le fichier $file");
    if(($file=~/SwissProt_ec.sql/)||($file=~/SwissProt_go.sql/)){
    	$basederef=1;#0 pour metacyc 1 pour swissprot
    }
    else{
    	$basederef=0;
    }
    while (<IN>){
    #ligne du type:AY292526	3.1.1.14	Ginkgo biloba (maidenhair tree)	Eukaryota_Viridiplantae_Streptophyta_Embryophyta_Tracheophyta_Spermatophyta_Ginkgophyta_Ginkgoales_Ginkgoaceae_Ginkgo	3311	MVLVKDVFSEGPLPVQILAIPQANSSPCSKLADKNGTATTPSPCRPPKPLLIALPSQHGDYPLILFFHGYVLLNSFYSQLLRHVASHGYIAIAPQMYSVIGPNTTPEIADAAAITDWLRDGLSDNLPQALNNHVRPNFEKFVLAGHSRGGKVAFALALGRVSQPSLKYSALVGLDPVDGMGKDQQTSHPILSYREHSFDLGMPTLVVGSGLGPCKRNPLFPPCAPQGVNHHDFFYECVAPAYHFVASDYGHLDFLDDDTKGIRGKATYCLCKNGEAREPMRKFSGGIVVAFLQAFLGDNRGALNDIMVYPSHAPVKIEPPESLVTEDVKSPEVELLRRAVCR	MetaCyc
    #001R_FRG3G	-	Frog virus 3 (isolate Goorha) (FV-3)	Viruses; dsDNA viruses, no RNA stage; Iridoviridae; Ranavirus	654924	MAFSAEDVLKEYDRRRRMEALLLSLYYPNDRKLLDYKEWSPPRVQVECPKAPVEWNNPPSEKGLIVGHFSGIKYKGEKAQASEVDVNKMCCWVSKFKDAMRRYQGIQTCKIPGKVLSDLDAKIKAYNLTVEGVEGFVRYSRVTKQHVAAFLKELRHSKQYENVNLIHYILTDKRVDIQHLEKDLVKDFKALVESAHRMRQGHMINVKYILYQLLKKHGHGPDGPDILTVKTGSKGVLYDDSFRKIYTDLGWKFTPL	Swiss-Prot	 GO:0006351, IEA:UniProtKB-KW.; GO:0006355, IEA:UniProtKB-KW.; GO:0046782, IEA:InterPro.
    #id	ec	nom	trucs	nb	sequence	basededonne	GOterms	
	chomp($_);
        @Line = split(/\t/,$_);
        $Line[2] = $1 if ($Line[2] =~ /^(.*) \(/);
        $Line[2] = $1 if ($Line[2] =~ /^([^ ]* [^ ]*) /);
	#$Line[2] contient le nom de l'espece
	# Si le taxonomy ID ($Line[4])est present dans FUNGIpath
	if (exists($IDtaxon_FUNGIpath{$Line[4]})){
		# On recupere le nom que l'on lui a asscie dans FUNGIpath (car pour un taxonomy ID peut avoir plusieurs noms)
		$Line[2] = $IDtaxon_FUNGIpath{$Line[4]}{'complet'};
	}
	#si l'espece est presente dans FUNGIpath (taxonomy ID identique ou nom)
        if (exists($Specie{$Line[2]})){
	# Nom abrege de l'espece dans FUNGIpath
            $geno = $IDtaxon_FUNGIpath{$Specie{$Line[2]}}{'abrev'}; 
	    # Recuperation de la sequence proteique 
            $Line[5] = lc($Line[5]); # seq prot, Tout minuscule
            $Line[5] =~ s/\*//g; # Suppression des exventuelles *
            #on a fait SpecieSeq comme ça:$SpecieSeq{$$Res{'genome0'}}{$$Res{'seqprot'}} = $$Res{'id0'};
            if (exists($SpecieSeq{$geno})){
		# Si une sequence identique pour le genome $geno est presente dans FUNGIpath
		# On connait directement les IDs equivalents
                if (exists($SpecieSeq{$geno}{$Line[5]})){
                    $ID_equivalent_identical{$geno}{$SpecieSeq{$geno}{$Line[5]}} = $Line[0];
		# Enregistrement des deux sequences equivalentes dans FUNGIpath et SwissProt/MetaCyc ainsi que des annotations trouvees
                    open OUT,'>>'.$file_annotated_seq_FUNGIpath;
                    #cecile: ajout d'une colonne pour indiqué d'ou vient l'annotation, metacyc ou swissprot
                    if($#Line<7){
		      $Line[7]="-";
                    }
                    print OUT "$Line[0]\t$SpecieSeq{$geno}{$Line[5]}\t$Line[1]\t$Line[7]\tIdentical\t$basederef\n";#cecile: basederef=0 pour metacyc et 1 pour swissprot
                    close(OUT);
                }
                else {
                    # Si pas de sequence identique
			# On va rechercher la sequence la plus similaire par BLAST
			# Enregistrement des informations relatives a $Line[0] (utilise losque le BLAST sera fait)
			my $nomtemp=$Line[0]."__$basederef";
			$Id_Blast{$nomtemp}{'seq'} = $Line[5];
			$Id_Blast{$nomtemp}{'ec'} = $Line[1];
			if($#Line<7){
			  $Line[7]="-";
			}
			$Id_Blast{$nomtemp}{'go'} = $Line[7];
										# On ajoute la sequence au fichier des sequences a comparer au genome $geno
                    open OUT,'>>'.$dirTmp.$geno.'_input_blast.fa';
                    print OUT ">$Line[0]__$basederef\n$Line[5]\n";#avec l'ajout a l'id de $basederef je garde l'info de la provenance de la sequence
                    close(OUT);
                }
            }
            else {
                die "Problem with $geno ($_)\n";
            }
        }
    }
    close(IN);
}

print "recuperation des fichiers de sequences a comparer\n";
# Recuperation des fichiers de sequences a comparer
chdir($dirTmp);
@File = glob("*_input_blast.fa");
foreach $file (@File) {
	# Recuperation du genome
    $geno = '';
    $geno = $1 if ($file =~ /^(.*)_input_blast.fa$/);

    # Construction de la base de donnees pour le BLAST
    open OUT,'>'.$dirTmp.$geno.'.fa' or die("impossible d'ouvrir ".$dirTmp.$geno.'.fa');
    foreach $i (keys(%{$SpecieSeq{$geno}})) {
			print "specie seq $geno $i : $SpecieSeq{$geno}{$i}\n";
			print OUT ">$SpecieSeq{$geno}{$i}\n$i\n";
    }
    close(OUT);
	
	# Formattage de la base
    #system("formatdb -p T -i $dirTmp$geno.fa -n $dirTmp"."blast_db_$geno");
	
		system("makeblastdb -dbtype prot -in $dirTmp$geno.fa -out $dirTmp"."blast_db_$geno");
		#unlink($dirTmp.$geno.'.fa');
		system("mv ".$dirTmp.$geno.'.fa'." ".$dirTmp.$geno.'vu.fa');
	
	# Realisation du BLAST
    #system("blastall -p blastp -i $dirTmp$file -d $dirTmp"."blast_db_$geno -m8 -e 0.001 > $dirTmp$geno.blast");
		system("blastp -query $dirTmp$file -db $dirTmp"."blast_db_$geno -evalue 0.001 -outfmt 6 -seg yes -out $dirTmp$geno.blast");
		#Suppression de la base 
		@Line = glob("blast_db_$geno*");
		foreach $i (@Line){
#			unlink($dirTmp.$i);
			system("mv ".$dirTmp.$i." ".$dirTmp.$i."vu");
		}

		$id = '';
		# Analyse des resultats BLAST
    open IN,$dirTmp.$geno.'.blast' or die("impossible d'ouvrir ".$dirTmp.$geno.'.blast');
    while (<IN>) {
        chomp($_);
        @Line = split(/\t/,$_);
				# Si on lit le meilleur resultat associe a $id
        unless ($Line[0] eq $id){
					$test = 0; # Vaut un si un des deux criteres remplis
					# Calul du pourcentage aligne
					print STDERR "$geno $Line[0]\n";
       		$perc_align = int((($Line[7] - $Line[6])/length($Id_Blast{$Line[0]}{'seq'}))*100); 
					# Si plus de 90% aligne et identique, on suppose que les sequences sont equivalentes
       	 	if (($Line[2] >= 90) && ($perc_align >= 90)) {
						$test = 1;
					}
					else{
						# Sinon, calcul du pourcentage total d'identite $i
						$perc_identity =  alignNeedle($Line[0],$Id_Blast{$Line[0]}{'seq'},$Line[1],$dirTmp);
						$test = 1 if ($perc_identity >= 85);
					}
					# Si une sequence correpond au critere
					if ($test == 1){
					# Si la sequence a deja ete associee a une autre sequence de FUNGIpath car identique,
					# On ne l'enregistre pas
            unless (exists($ID_equivalent_identical{$geno}) && exists($ID_equivalent_identical{$geno}{$Line[1]})){
						# On enregistre la sequence
               $ID_equivalent_similar{$geno}{$Line[1]}{$Line[0]} = $Line[10];
            }
        	}			
            # else {
                # print "Warning : Best Hit BLAST of $Line[0] is \t$Line[1]\t(which is inferior to threshold - [$Line[2] - $perc_align] and $perc_identity)\n";
            # }
        }
        $id = $Line[0]; # Pour savoir si c'est le premier resultat de Line[0] qu'on lit ou non
    }
    close(IN);

	# Pour chaque sequence trouvee par similarite
	if (exists($ID_equivalent_similar{$geno})){
		open OUT,'>>'.$file_annotated_seq_FUNGIpath;
		foreach $id (keys(%{$ID_equivalent_similar{$geno}})) {
			# Recuperation des IDs trouves comme etant equivalent a $id dans SwissProt/MetaCyc
			@Line = keys(%{$ID_equivalent_similar{$geno}{$id}});
			$nomaecrire=$Line[0];
			$nomaecrire=~s/__(.*)//;
			$basederef=$1;
			# Si un seul ID trouve, on l'enregsitre
			if ($#Line == 0) {
				print OUT "$nomaecrire\t$id\t$Id_Blast{$Line[0]}{'ec'}\t$Id_Blast{$Line[0]}{'go'}\tSimilar\t$basederef\n";
			}
			# Sinon, on le signale
			else {
				$i = '';
				$test = 100;
				# Recherche de la sequence avec la meilleure E-val
				foreach $name (@Line) {
					if ($ID_equivalent_similar{$geno}{$id}{$name} < $test){
						$test = $ID_equivalent_similar{$geno}{$id}{$name};
						$i = $name;
					}
				}
				#print "Warning: $id\t".join(',',@Line)."\t$i save\n";
				$nomaecrire=$i;
				$nomaecrire=~s/__(.*)//;
				$basederef=$1;
				print OUT "$nomaecrire\t$id\t$Id_Blast{$i}{'ec'}\t$Id_Blast{$i}{'go'}\tSimilar_PlsID\t$basederef\n";#cecile modification de la ligne: si plusieurs ID trouvé on garde celui avec la meilleur evalue
			}
		}
		close(OUT);
	}
}

#######################################################################
# Fonction
###########
# Fonction qui aligne deux sequences et extrait le pourcentage d'identite
sub alignNeedle{
    my ($idSp,$seq,$idDb,$repOut) = @_;
    my @Line = ();
    my @Res = ();
    my $fileOut = $repOut.$idSp.'_'.$idDb.'.needle';
    my $i = 0;
    my $req = '';
	# Realisation de l'alignement avec needle a moins qu'il est deja ete realise
    unless (-e $fileOut) {
		# Ecriture du fichier FASTA de la premiere sequence si fichier non cree
        unless (-e $repOut.$idSp.'.fa'){
            open OUT,'>'.$repOut.$idSp.'.fa';
            print OUT '>'.$idSp."\n$seq\n";
            close(OUT);
        }
		# Ecriture du fichier FASTA de la seconde sequence si fichier non cree
        unless (-e $repOut.$idDb.'.fa') {
#            $req = $db->prepare('SELECT seqprot FROM seq WHERE id0 = \''.$idDb.'\';');
#            $req->execute();
#            while (@Res = $req->fetchrow_array()){
#                $seq = $Res[0];
#            }
	    $seq=$idseq{$idDb};
            $seq =~ s/\*//g;
            open OUT,'>'.$repOut.$idDb.'.fa';
            print OUT '>'.$idDb."\n$seq\n";
            close(OUT);
        }
		# Realisation de l'alignement
        system('needle -aseq '.$repOut.$idSp.'.fa -bseq '.$repOut.$idDb.'.fa -gapopen 10 -gapextend 0.5 -outfile '.$fileOut) unless ((-e $fileOut) && !(-z $fileOut));
    }
	# On extrait le pourcentage d'identite
    $i = `grep Identity $fileOut`;
    return int($2) if ($i =~ /Identity: *[0-9]*\/([0-9]*) *\((.*)%\)/);
    return 0;
}
