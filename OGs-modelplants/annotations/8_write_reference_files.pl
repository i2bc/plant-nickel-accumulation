#!/usr/bin/perl -w
use strict;
use LWP::Simple;
use LWP::UserAgent;

#24/04/2013

#perl 8_intersection.pl [repertoire de travail] [version du dossier de metacyc]> 8_intersection.out
#perl 8_intersection.pl /home/cecile/git/fungipath/miseajourBD/local/intersection/ 16.1> 8_intersection.out

# Description

#######################################################################
# LIST OF VARIABLES
###################

my $dirWork = $ARGV[0]; # INPUT : Repertoire de travail
my $dirIn = $dirWork.'/annot_source/';  # INPUT : Repertoire contenant les bases de données Swiss-Prot et MetaCyc
my $fileIn_SwissProt = $dirIn.'uniprot_sprot.dat';  # INPUT : Swiss-Prot au format EMBL
my $dirIn_MetaCyc = $dirIn.'/'.$ARGV[1].'/';  # Metacyc: all downloaded flat files (meta.tar.gz)
my $dirOut = $dirWork.'/annot_files/'; # OUTPUT : repertoire contenant les fichiers necessaires a l'annotation des groupes
my $fileOut_SwissProt = $dirOut.'SwissProt_ec.sql'; # OUTPUT : Base de données des sequences annotees dans Swiss-Prot avec des ECs complets
my $fileOut_SwissProt_go = $dirOut.'SwissProt_go.sql'; # OUTPUT : Base de données des sequences annotees dans Swiss-Prot avec des Go
my $fileOut_MetaCyc = $dirOut.'MetaCyc_ec.sql'; # OUTPUT : Base de données des sequences annotees dans MetaCyc avec des ECs complets
my $fileOut_All_sql = $dirOut.'SwissProt_MetaCyc_ec.sql'; # OUTPUT : Base de données des sequences annotees dans MetaCyc et SwissProt avec des ECs complets - sequences identiques en un exemplaire
my $fileOut_All_sql_go = $dirOut.'SwissProt_MetaCyc_go.sql'; # OUTPUT : Base de données des sequences annotees dans MetaCyc et SwissProt avec des ECs complets - sequences identiques en un exemplaire
my $fileOut_All_fa = $dirOut.'SwissProt_MetaCyc_ec.fa'; # OUTPUT : Base de données des sequences annotees dans MetaCyc et SwissProt avec des ECs complets - sequences identiques en un exemplaire
my $fileOut_All_fa_go = $dirOut.'SwissProt_MetaCyc_go.fa'; # OUTPUT : Base de données des sequences annotees dans MetaCyc et SwissProt avec des ECs complets - sequences identiques en un exemplaire

# Autres variables
my @File = ();
my @Tmp = ();
my $id = '';
my $ec = '';
my $seq = '';
my $max = 0;
my $i = '';
my %Seq = ();
my %Ec = ();
my %Test = ();
my $goterm="";
my @Tmpprgot=();
my $tmpenr="";



#######################################################################
# Programme
###########
# suppression des anciens fichiers de sortie:
if (-e $fileOut_SwissProt){
 	unlink($fileOut_SwissProt);
}
if(-e $fileOut_MetaCyc){
	unlink ($fileOut_MetaCyc);
} 
if(-e $fileOut_All_sql){
	unlink ($fileOut_All_sql);
}
if(-e $fileOut_All_fa){
	unlink($fileOut_All_fa);
}

# Creation des repertoires
mkdir($dirOut) unless (-e $dirOut);



# Extraire les sequences annotees avec un EC number complet de SwissProt
if (-e $fileIn_SwissProt) {
     ExtractSeqSwissProtbanquesdiff($fileIn_SwissProt,$fileOut_SwissProt,$fileOut_SwissProt_go);
}
else {
     die("No file $fileIn_SwissProt\n");
}
#Extraire les sequences annotees avec un EC number complet de MetaCyc
if (-e $dirIn_MetaCyc) {
   ExtractSeqMetaCyc($dirIn_MetaCyc,$fileOut_MetaCyc);
}
else {
   die("No directory $dirIn_MetaCyc\n");
}


# Comparaison des deux bases de donnees pour ne pas avoir deux sequences identiques
# à faire aussi avec les annotations GO?
# par la suite deux fichiers. Un avec les ec et un avec les go.
@File = ($fileOut_SwissProt,$fileOut_MetaCyc,$fileOut_SwissProt_go);#lit les deux fichiers sql fait à l'étape précédente
foreach my $file (@File){
	if (-e $file){
		# Lecture de chaque fichier
		open IN,$file or die;
		while (<IN>)    {
		chomp($_);
			@Tmp = split(/\t/,$_);
			if ($#Tmp > 4) {
				$Tmp[5] =~ s/\*$//;
				# $Seq{sequence_aa}{liste des ECs}{ID}{liste des GO} = '';
				#si il y a des go term et aps d'ec? =>EC remplace par -
				#si il y a des ec et pas de go term?
				#on ne va pas jusqu'au niv 4
				if($#Tmp>6){
				  $Seq{uc($Tmp[5])}{$Tmp[1]}{$Tmp[0]} = $Tmp[7] if ($Tmp[5] ne '');
				}
				else{
				  $Seq{uc($Tmp[5])}{$Tmp[1]}{$Tmp[0]} = '' if ($Tmp[5] ne '');
				}
				
			}
			else {
				print "Pb with $file ($_)\n".join('/',@Tmp)."\n";
			}
		}
		close(IN);
	}
	else{
		die ("Erreur : le fichier $file n'a pas ete trouve.\n");
	}
}

# Ecriture des fichiers finaux en supprimant les sequences identiques pour ne pas avoir de doublons
open FA,'>'.$fileOut_All_fa;
open SQL,'>'.$fileOut_All_sql;
open FAGO,'>'.$fileOut_All_fa_go;
open SQLGO,'>'.$fileOut_All_sql_go;
# Pour chaque sequence proteique
foreach $seq (keys(%Seq)) {#pour chaque sequence
    #print $seq."\n";
    #Recuperation des listes d'ECs
    @Tmp = keys(%{$Seq{$seq}});
    $id = '';
    $max = -1;
    %Ec = ();
    $goterm="";
	# Si plus d'une liste, c'est qu'elle differe entre les deux bases
    if ($#Tmp > 0){
        print "\nSeveral annotation for the same sequence:\n";
        #recuperation de la liste d'ECs trouves par chaque base
        foreach $i (keys(%{$Seq{$seq}})){
            print "$i\t".join(',',keys(%{$Seq{$seq}{$i}}))."\n";
            @Tmp = split(/,/, $i); # Recuperation des ECs
            @Tmpprgot=keys(%{$Seq{$seq}{$i}});
            $id=pop(@Tmpprgot);
            if(exists($Seq{$seq}{$i}{$id})){
	      $goterm.=$Seq{$seq}{$i}{$id};#je conserve les go terms
            }
            #on va faire l'union des annotations
            foreach $ec (@Tmp){
		unless($ec eq '-'){
		  	$Ec{$ec} = '' ;
		}
            }
	# On conservera l'ID possedant le plus d'ECs
            if ($#Tmp > $max)  {
               $max = $#Tmp;
               @Tmp = keys(%{$Seq{$seq}{$i}});
               $id = pop(@Tmp);
            }
        }
	# Liste des ECs trouves pour la sequence $seq
        $ec = join(',',sort(keys(%Ec)));
        print "Save\t$id\t$ec\t$goterm\n";
    }
		# Si une liste d'ECs trouves
    else{
        $ec = pop(@Tmp); # Recuperation de la liste d'ECs
				# Recuperation de l'ID de la sequence
        @Tmp = keys(%{$Seq{$seq}{$ec}});
        $id = pop(@Tmp);
        if(exists($Seq{$seq}{$ec}{$id})){
	  			$goterm=$Seq{$seq}{$ec}{$id};
				}
				else{
	  			$goterm="-";
				}
    }
# On verifie qu'une sequence n'a pas plusieurs annotations
# car si un ID pst dans les deux dbs mais avec des sequences aa differentes
    if (exists($Test{$id})){
	print "Warning : $id gets several sequences\n";
        unless ($Test{$id} eq $ec){
            print "Problem: $id with different annotation and sequence ($ec and $Test{$id})\n";
        }
    }
    else{
        $Test{$id} = $ec;
        chomp($id);
        if(($ec eq '')||($ec eq '-')){#pas d'ec
        	print FAGO '>'.$id."|$ec|$goterm\n$seq\n";
		print SQLGO "$id\t$ec\t$goterm\n";
        }
        else{
        	if(($goterm eq '-')||($goterm eq '')){#pas de go
        		print FA '>'.$id."|$ec|$goterm\n$seq\n";
        		print SQL "$id\t$ec\t$goterm\n";
        		print "presence d'ec sans go terms\n";
        	}
        	else{#ec et go
        		print FAGO '>'.$id."|$ec|$goterm\n$seq\n";
			print SQLGO "$id\t$ec\t$goterm\n";
			print FA '>'.$id."|$ec|$goterm\n$seq\n";
        		print SQL "$id\t$ec\t$goterm\n";
        	}
        }
    }
}
close(FA);
close(SQL);


###########
#FUNCTION
###########
sub ExtractSeqSwissProtbanquesdiff{
   print STDERR "ExtractSeqSwissProt banques diff pour ec et go\n";
    my ($fileIn,$fileOut,$fileOutGO) = @_;
    my $id = '';
    my $taxonId = 0;
    my $de = '';
    my $taxon = '';
    my $organism = '';
    my $i = '';
    my $seq = '';
    my @Tmp = ();
    my @Ec = ();
    my @Organism = ();
    my @Taxon = ();
    my @GO=();#Goterm,qualité
    my @tmpGO=();
#print STDERR "ExtractSeqSwissProt\n";
# Creation du fichier SQL des sequences annotees avec un EC
    open SQL,'>'.$fileOut;
    open SQLGO,'>'.$fileOutGO;
    #ajouter des infos pour les go sur la meme ligne
    #mais attention prendre aussi la qualité de l'annotation GO
    #lignes du type:
    #DR   GO; GO:0006355; P:regulation of transcription, DNA-dependent; IEA:UniProtKB-KW.
    
    # Lecture du fichier 
    open IN,$fileIn;    
    while (<IN>) {
        chomp($_);
	# Recuperation de l'ID de la sequence
        if ($_ =~ /ID +([^_ ]*_[^ ]*) /) {
            $id = $1;
        }
	# Recuperation des annotations enzymatiques completes
        elsif ($_ =~ /^DE *([^ ].*)$/) {
					$de = $1 ;
					if ($de =~ /EC=/){
						#if ($de =~ /EC=([0-9\.]*);/){
						if($de=~/EC=([0-9\.]*)[ ;]+/){
							push(@Ec,$1);
						}
					else{
#Problem with EC DE            EC=2.4.1.333 {ECO:0000269|PubMed:24647662}; (12OLP_LISIN)
						die("Problem with EC $_ ($id)\n") unless ($de =~ /EC=([0-9\.\-n]*)[ ;]+/);
					}
		}

	}
	# Recuperation de l'organisme
	elsif($_ =~ /^OS *([^ ].*)$/){
		push (@Organism, $1) ;
	}
	# Recuperation du taxon ID
	elsif($_ =~ /^OX *NCBI_TaxID=([0-9]*);/){
		$taxonId = $1 ;
	}
	# Recuperation du taxon
	elsif($_ =~ /^OC *([^ ].*)$/){
		push (@Taxon, $1);
	}
	elsif ($_ =~ /^SQ /)  {
            $_ = <IN>;
            while ($_ =~ /^ /) {
		chomp($_);
                $seq .= $_;
                $_ = <IN>;                
            }
            $seq =~ s/ //g; # Suppression des espaces dans la sequence
        }
        elsif ($_ =~/^DR\s*GO;/){
	    @tmpGO=split(/;/,$_);
	    $tmpenr="$tmpGO[1],$tmpGO[3]";
	    push(@GO,$tmpenr);
        }
	if ($_ =~ /^\/\//) { 
		# On conserve la sequence si au moinq un EC complet ou un GO
		if (($#Ec >= 0)||($#GO>=0)){
		  $organism = join(' ', @Organism);
		  $taxon = join(' ', @Taxon);
		  # $organism = $1 if ($organism =~ /^([^\(]*) \(/);
		  $taxon =~ s/\.$//;
		  $organism =~ s/\.$//;
		  chomp($id);
		  if($#Ec==-1){#si il n'y a pas d'ec
		    print SQLGO "$id\t-\t$organism\t$taxon\t$taxonId\t$seq\tSwiss-Prot\t".join(';',sort(@GO))."\n" ;
		  }
		  else{
		  	if($#GO==-1){#si il n'y a pas de go
		  		print "ec sans go trouve\n";
		    	$GO[0]="-";
		    	print SQL "$id\t".join(',',sort(@Ec))."\t$organism\t$taxon\t$taxonId\t$seq\tSwiss-Prot\n" ;
		  	}
		  	else{#présence a la fois de ec et de go
		  		print SQL "$id\t".join(',',sort(@Ec))."\t$organism\t$taxon\t$taxonId\t$seq\tSwiss-Prot\t".join(';',sort(@GO))."\n" ;
		  		print SQLGO "$id\t".join(',',sort(@Ec))."\t$organism\t$taxon\t$taxonId\t$seq\tSwiss-Prot\t".join(';',sort(@GO))."\n" ;
				}
			}
		}
		$id = '';
		$taxonId = '';
		$seq = '';
		@Ec = ();
		@Organism = ();
		@Taxon = ();  
		@GO=();
	}      
    }
    close(IN);
    close(SQL);
    close(SQLGO);
}
# Extrait les sequences annotees avec un EC number complet
# et les sequences annotees avec des GOterms
sub ExtractSeqSwissProt{
    print STDERR "ExtractSeqSwissProt\n";
    my ($fileIn,$fileOut) = @_;
    my $id = '';
    my $taxonId = 0;
    my $de = '';
    my $taxon = '';
    my $organism = '';
    my $i = '';
    my $seq = '';
    my @Tmp = ();
    my @Ec = ();
    my @Organism = ();
    my @Taxon = ();
    my @GO=();#Goterm,qualité
    my @tmpGO=();
#print STDERR "ExtractSeqSwissProt\n";
# Creation du fichier SQL des sequences annotees avec un EC
    open SQL,'>'.$fileOut;

    #ajouter des infos pour les go sur la meme ligne
    #mais attention prendre aussi la qualité de l'annotation GO
    #lignes du type:
    #DR   GO; GO:0006355; P:regulation of transcription, DNA-dependent; IEA:UniProtKB-KW.
    
    # Lecture du fichier 
    open IN,$fileIn;    
    while (<IN>) {
        chomp($_);
	# Recuperation de l'ID de la sequence
        if ($_ =~ /ID +([^_ ]*_[^ ]*) /) {
            $id = $1;
        }
	# Recuperation des annotations enzymatiques completes
        elsif ($_ =~ /^DE *([^ ].*)$/) {
		$de = $1 ;
		if ($de =~ /EC=/){
			if ($de =~ /EC=([0-9\.]*);/){
				push(@Ec,$1);
			}
			else{
				die("Problem with EC $_ ($id)\n") unless ($de =~ /EC=([0-9\.\-n]*);/);
			}
		}

	}
	# Recuperation de l'organisme
	elsif($_ =~ /^OS *([^ ].*)$/){
		push (@Organism, $1) ;
	}
	# Recuperation du taxon ID
	elsif($_ =~ /^OX *NCBI_TaxID=([0-9]*);/){
		$taxonId = $1 ;
	}
	# Recuperation du taxon
	elsif($_ =~ /^OC *([^ ].*)$/){
		push (@Taxon, $1);
	}
	elsif ($_ =~ /^SQ /)  {
            $_ = <IN>;
            while ($_ =~ /^ /) {
		chomp($_);
                $seq .= $_;
                $_ = <IN>;                
            }
            $seq =~ s/ //g; # Suppression des espaces dans la sequence
        }
        elsif ($_ =~/^DR\s*GO;/){
	    @tmpGO=split(/;/,$_);
	    $tmpenr="$tmpGO[1],$tmpGO[3]";
	    push(@GO,$tmpenr);
        }
	if ($_ =~ /^\/\//) { 
		# On conserve la sequence si au moinq un EC complet ou un GO
		if (($#Ec >= 0)||($#GO>=0)){
		  $organism = join(' ', @Organism);
		  $taxon = join(' ', @Taxon);
		  # $organism = $1 if ($organism =~ /^([^\(]*) \(/);
		  $taxon =~ s/\.$//;
		  $organism =~ s/\.$//;
		  chomp($id);
		  if($#Ec==-1){
		    $Ec[0]="-";
		  }
		  if($#GO==-1){
		    $GO[0]="-";
		  }
		  print SQL "$id\t".join(',',sort(@Ec))."\t$organism\t$taxon\t$taxonId\t$seq\tSwiss-Prot\t".join(';',sort(@GO))."\n" ;
		}
		$id = '';
		$taxonId = '';
		$seq = '';
		@Ec = ();
		@Organism = ();
		@Taxon = ();  
		@GO=();
	}      
    }
    close(IN);
    close(SQL);
}

#extraire les seq de metacyc
sub ExtractSeqMetaCyc
{
    print STDERR "ExtractSeqMetaCyc\n";
    my ($dirIn,$fileOutMetaSql) = @_;
    # INPUT : Fichiers de MetaCyc utilises
    my $fileProt = $dirIn.'data/proteins.dat';
    my $fileReac = $dirIn.'data/reactions.dat';
    my $fileEnzReac = $dirIn.'data/enzrxns.dat';
    #print STDERR "ExtractSeqMetaCyc\n";
    #  INPUT : Liste des Bases de donnees dont les sequences devront etre extraites
    # Car les sequences proteiques ne ficgurent pas dans MetaCyc, juste les identifiants
    my %Database = ();
    # Ordre de recherche dans les bases
    # Si on a lien vers differentes bases, on va par exemple privilegier la base SwissProt par rapport a ref seq
    $Database{'SWISSPROT'} = 1; $Database{'UNIPROT'}= 2 ; $Database{'REFSEQ'}= 3; $Database{'PID'}=4; $Database{'PDB'}=5;

    # Autres variables	
    my $ec = '';
    my $id = '';
    my $enzyme = '';
    my $component = '';
    my $i = 0;
    my $reaction = '';
    my $tmp = '';
    my $idTmp = '';
    my $complex = '';
    my %ReactionEc = ();
    my %EnzymaticReaction = ();
    my %ProteinReaction = ();
    my %ProteinComplex = ();
    my %ComplexProtein = ();    
    my %ProteinLink = ();
    my %SearchSeq = ();
    my %Specie = ();
    my @Reaction = ();
    my @Database = ();
	
    # Tri des bases de donnees dans l'ordre de preference
    my @Database = ();
    foreach my $database (keys(%Database))
    {
		$Database[$Database{$database} - 1] = $database;
    }
	
    # Lecture du fichier decrivant les reactions chimiques
    if (-e $fileReac) 
    {
		print STDERR "lecture du fichier decrivant les reactions chimiques\n";
        open IN,$fileReac or die ("impossible d'ouvrir le fichier $fileReac\n");
        while (<IN>) 
        {
        	chomp($_);
			if ($_ =~ /EC-NUMBER - EC-([0-9\.]*)$/)
			{
				if ($ec eq '')
				{
					$ec = $1 ;
				}
				else
				{
					$ec.=" ".$1;
					print STDERR "Several EC for the same enzymatic reaction ($ec and $1)\n";
				}
			}
			push(@Reaction, $1) if ($_ =~ /ENZYMATIC-REACTION - (.*)$/); 
			if ($_ eq '//')
			{
				# Si un EC complet a ete attribue a la reaction, on enregistre la reaction avec son activite
				if ($ec ne '')
				{
					foreach $i (@Reaction)
					{
						$ReactionEc{$i}{$ec} = '' ;
					}
				}
				# Reinitialisation des varaibles
				$ec = '' ;
				@Reaction = ();
			}
        }
        close(IN);
    }
    else { die ("No file $fileReac\n"); }

	# Lecture du fichier decrivant les reactions enzymatiques
    if (-e $fileEnzReac) 
    {
    	print STDERR "lecture du fichier decrivant les reactions enzymatiques\n";
        open IN,$fileEnzReac or die("impossible d'ouvrir le fichier $fileEnzReac\n");
        while (<IN>) 
        {
            chomp($_);
            $id = $1 if ($_ =~ /UNIQUE-ID - (.*)$/); # ID de la reaction
            $enzyme = $1 if ($_ =~ /ENZYME - (.*)$/); # Nom de l'enzyme
            if ($_ eq '//') 
            {
                $EnzymaticReaction {$id} {$enzyme} = '' if ($enzyme ne '');
                $id = '';
                $enzyme = '';
            }
        }
        close(IN);
    }
    else { die ("No file $fileEnzReac\n"); }

	# Lecture du fichier decrivant toutes les proteines
    if (-e $fileProt) 
    {
        open IN,$fileProt or die ("pb avec l'ouverture du fichier $fileProt\n");
        print STDERR "lecture du fichier decrivant toutes les proteines\n";
        while (<IN>) 
        {
            chomp($_);
            $id = $1 if ($_ =~ /UNIQUE-ID - (.*)$/); # ID de la proteines
            $ProteinReaction {$id}{$1} = '' if ($_ =~ /^CATALYZES - (.*)$/); # ID de la reaction enzymatique qu'elle catalyse
			if ($_ =~ /^DBLINKS/)
			{
				chomp($_);
		#print STDERR "dans DBLINKS $_\n";
		#DBLINKS - (PDB "1AJR" NIL |fulcher| 3381769294 NIL NIL)
		if ($_ =~ /^DBLINKS - *\(([^ ]*) "([^\"]*)".*$/){# ligne du type DBLINKS - (qqchose "qqchose"qqchose
			$database = $1;
			$idTmp = $2;
		}
		else{
			#DBLINKS - (KEGG "
			#/Synpcc7942_0535" NIL |brito| 3513696542 NIL NIL)
			if($_ =~ /^DBLINKS - *\(([^ ]*) ".$/){#cas pas géré par sandrine, ajout cecile
				#print STDERR "if cecile\n";
				$database=$1;
				#print STDERR $database."\n";
				while(!(/^\/[^ ]*" .*\)$/)){
					$_ = <IN>;#on regarde la ligne suivante
					if($_ =~/^\/([^ ])*" .*\)$/){
					  $idTmp=$1;
					  die ("Problem with $id in RefSeq ($idTmp in $database)\n") unless ($_ =~ /^\//);
					}
				}
									
			}
			else{
			  #DBLINKS - (PID "AAB70707
			  #/" NIL |christ| 3319379217 NIL NIL)
			  if ($_ =~ /^DBLINKS - *\(([^ ]*) "([^\"]*)$/){# ligne du type :DBLINKS - (qqchose "qqchoseqqchose =>pas de fermeture du "
				$database = $1;
				$idTmp = $2;
				while (!($_ =~ /^\/"[^\"]*\)$/)) {#tant que la ligne n'est pas du type : /"qqchose) ) => ne gere pas le cas 				
					$_ = <IN>;#on regarde la ligne suivante
					#print STDERR "ligne suivante regardée $_\n";
					#print STDERR "$id\n $idTmp \n $database\n";
					die ("Problem with $id in RefSeq ($idTmp in $database)\n") unless ($_ =~ /^\//);
				}
			  }	
			  else{
				print "Warning : Ligne DBLINKS dans un format different que celui attendu\n$_\n";
				exit();
			  }
			}
		}
		# Si lien vers une des bases que l'on a selectionne, on conserve le lien
		# Base de donnnees selectionnees si possible de telecharger la sequence automatiquent avec un script
		$ProteinLink{$id} {$database} = $idTmp if (exists($Database{$database})); 
	    }			
	# Si la proteine appartient a un complexe
    if ($_ =~ /COMPONENT-OF - (.*)$/) {
                $component = $1 ; #Nom du complexe
                $ComplexProtein {$component} {$id} = ''; # Pour connaitre proteine qui constituent le complexe
		$ProteinComplex {$id} {$component} = '';
	    }
	  # Si la proteine est un complexe
            if ($_ =~ /COMPONENTS - (.*)$/) {
                $component = $1 ; #Nom du complexe
                $ComplexProtein  {$id}{$component} = ''; # Pour connaitre proteine qui constituent le complexe
		$ProteinComplex  {$component}{$id} = '';
	    }
            $id = '' if ($_ eq '//');
        }
        close(IN);
    }
    else {
        die "No file $fileProt\n";
    }

    # Pour chaque reaction enzymatique
    foreach my $enzReacn (keys(%EnzymaticReaction)) {
    	print STDERR "dans enz reacn\n";
        #si reaction a un ec 
        if (exists($EnzymaticReaction {$enzReacn} )) {
	# Pour chaque enzyme associee a la reaction
            foreach $enzyme (keys(% {$EnzymaticReaction {$enzReacn} })) {
		print STDERR "dans le foreach\n";
                $enzyme =~ s/\|//g unless (exists($ProteinReaction {$enzyme}));
		 # Si l'enzyme correspond a une proteine
                if (exists($ProteinReaction{$enzyme})) {
		    print "existe protein reaction\n";
		# Si la reaction a un EC complet
                    if (exists($ReactionEc{$enzReacn})) {
			print "existe reaction ec\n";
			# Reherche des proteines catalysant cette reaction
			# Si la proteine possede des liens vers d'autres bases
                            if (exists($ProteinLink{$enzyme})) {
				print "existe protein link\n";
				# On recherche la base de donnees dans laquelle la sequence sera recherchee
				DATABASE1: for ($i = 0; $i <= $#Database; $i ++){
					print "$Database[$i]\n";
					if (exists($ProteinLink{$enzyme}{$Database[$i]})) {
						print "dans protein link\n";
						# 1 : {Base de donnees dans laquelle la sequence doit etre recherchee}
						# 2 : {ID a rechercher}
						# valeur : liste des ECs complets associes a l'ID
						$SearchSeq {$Database[$i]} {$ProteinLink {$enzyme}{$Database[$i]}}  = join(',',sort(keys(% {$ReactionEc {$enzReacn}})));
						last DATABASE1;
					}
				}
			}
			  # Si l'enzyme est un complexe
			if (exists($ComplexProtein {$enzyme})) {
			    # Pour chaque proteine du complexe
			  foreach $id (keys(% {$ComplexProtein {$enzyme}})) {																	
				    if (exists($ProteinLink {$id})) {
					    # On recherche la base de donnees dans laquelle la sequence sera recherchee
					    DATABASE2: for ($i = 0; $i <= $#Database; $i ++){
					      if (exists($ProteinLink {$id}{$Database[$i]})) {
						    # 1 : {Base de donnees dans laquelle la sequence doit etre recherchee}
						    # 2 : {ID a rechercher}
						    # valeur : liste des ECs complets associes a l'ID
						    $SearchSeq {$Database[$i]} {$ProteinLink {$id}{$Database[$i]}}  = join(',',sort(keys(% {$ReactionEc {$enzReacn}})));
						    last DATABASE2;
					      }
					    }
				    }
			    }
			}
		      }
                }
                else {
                    print "Warning: No sequence for $enzyme in $fileProt\n";
                }
            }
        }
    }
	
    $i = 0;
	# Pour lister les especes afin de recuperer ensuite leur taxon
    $Specie{'-'} = '-';
	
	# Recuperation des sequences dans la base de donnnees SWISSPROT

    foreach $id (keys(% {$SearchSeq {'SWISSPROT'}})) {
    	print STDERR "recuperation des sequences dans la base de données SWISSPROT\n";
        Extract_UP_seq($id,$SearchSeq {'SWISSPROT'} {$id},$fileOutMetaSql,\%Specie);
	$i ++;
	# On fait une pause toutes les 10 secondes pour ne pas saturer le site internet
        if ($i > 10){
            $i = 0;
            sleep(2);
        }
    }
	# Recuperation des sequences dans la base de donnnees REFSEQ
    foreach $id (keys(% {$SearchSeq {'REFSEQ'}})) {
    		print STDERR "recuperation des sequences dans la base de données Refseq\n";
        Extract_pid_refseq_seq($id,$SearchSeq {'REFSEQ'} {$id},$fileOutMetaSql,\%Specie);
        $i ++;
		# On fait une pause toutes les 10 secondes pour ne pas saturer le site internet
        if ($i > 10){
            $i = 0;
            sleep(2);
        }
    }
	# Recuperation des sequences dans la base de donnnees PID
    foreach $id (keys(% {$SearchSeq {'PID'}})) {
    	print STDERR "recuperation des sequences dans la base de données PID\n";
        Extract_pid_refseq_seq($id,$SearchSeq {'PID'} {$id},$fileOutMetaSql,\%Specie);
        $i ++;
		# On fait une pause toutes les 10 secondes pour ne pas saturer le site internet
        if ($i > 10){
            $i = 0;
            sleep(2);
        }
    }
	# Recuperation des sequences dans la base de donnnees UNIPROT
    foreach $id (keys(% {$SearchSeq {'UNIPROT'}})) {
    		print STDERR "recuperation des sequences dans la base de données UNIPROT\n";
        Extract_UP_seq($id,$SearchSeq {'UNIPROT'} {$id},$fileOutMetaSql,\%Specie);
        $i ++;
	# On fait une pause toutes les 10 secondes pour ne pas saturer le site internet
        if ($i > 10){
            $i = 0;
            sleep(2);
        }
    }
    #recuperation des sequences dans la base de données PDB
    foreach $id (keys(% {$SearchSeq {'PDB'}})) {
   	print STDERR "recuperation des sequences dans la base de données PDB\n";
        Extract_PDB_seq($id,$SearchSeq {'PDB'} {$id},$fileOutMetaSql,\%Specie);
        $i ++;
		# On fait une pause toutes les 10 secondes pour ne pas saturer le site internet
        if ($i > 10){
            $i = 0;
            sleep(2);
        }
    }
}

# Extraction de la sequence $id de la base de donnnees Uniprot
sub Extract_UP_seq{
	my ($id,$ec,$fileOut,$refSpecie) = @_;
	my $page_html = '';
	my $i = 0;
	my $taxonId = 0;
    my $seq = '';
	my @Line = (); 
    my @Organism = ();	
    my @Taxon = ();
    #print STDERR "dans extract up seq \n";
	#Si recuperation de la sequence sans erreur
    if (($page_html = get('http://www.uniprot.org/uniprot/'.$id.'.txt')) || ($page_html = get('http://www.uniprot.org/uniprot/?query=replaces:'.$id.'&format=txt'))) {
		#On recupere chaque ligne dans un tableau
		@Line = split(/\n/,$page_html);
		for ($i = 0; $i <= $#Line; $i ++){
			# Recuperation de l'organisme
			if($Line[$i] =~ /^OS *([^ ].*)$/){
				push (@Organism, $1) ;
			}
			# Recuperation du taxon ID
			elsif($Line[$i] =~ /^OX *NCBI_TaxID=([0-9]*);/){
				$taxonId = $1 ;
			}
			# Recuperation du taxon
			elsif($Line[$i] =~ /^OC *([^ ].*)$/){
				push (@Taxon, $1);
			}
			elsif ($Line[$i] =~ /^SQ /)  {
				$i ++;
				while ($Line[$i] =~ /^ /) {
					$seq .= $Line[$i];
					$i ++;;                
				}
				$seq =~ s/ //g; # Suppression des espaces dans la sequence
			}
		}
		print "Warning no organism for $id in Uniprot\n" if ($#Organism < 0);
		print "Warning no taxonomy for $id in Uniprot\n" if ($#Taxon < 0);
		print "Warning no taxonomy ID for $id in Uniprot\n" if ($taxonId == 0);
		print "Warning no sequence for $id in Uniprot\n" if ($seq eq '');
		$$refSpecie{join(' ', @Organism)} = join(' ', @Taxon)."\t$taxonId" unless (exists($$refSpecie{join(' ', @Organism)}));
		open OUT,'>>'.$fileOut;
		chomp($id);
		print OUT "$id\t$ec\t".join(' ', @Organism)."\t".$$refSpecie{join(' ', @Organism)}."\t$seq\tMetaCyc\n" ;
		close(OUT);
	}
	else{
		print "Warning : Page not found for $id :\nhttp://www.uniprot.org/uniprot/$id\.txt\n'http://www.uniprot.org/uniprot/?query=replaces:$id&format=txt'\n";
	}
}


# Extraction de la sequence $id de la base de donnnees RefSeq
sub Extract_pid_refseq_seq{
		#print STDERR "dans extract pid refseq seq \n";
    my ($id,$ec,$fileOut,$refSpecie) = @_;
    my @Line = ();
    my @Seq = ();
    my $seq = '';
    my $tmpSeq = '';
    my $specie = '-';
    my $tmpSpecie = '-';
    my $taxon = '';
    my $i = 0;
	
    my $ua2 = new LWP::UserAgent;
    my $url2 ="http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&id=$id&retmode=text&rettype=gp&email=sandrine.grossetete\@igmors.u-psud.fr";
		#print STDERR $url2."\n";die;
	# Recuperation de la page associee a la sequence
    my $request2 = new HTTP::Request GET=> $url2;
    my $page2 = $ua2->request($request2);
    while (!($page2->is_success) && ($i < 5)) {
        sleep(1);
        $request2 = new HTTP::Request GET=> $url2;
        $page2 = $ua2->request($request2);
        $i ++;
    }
	# Si il y a eu une erreur lors de la recuperation de la sequence id
    if (!($page2->is_success))    {
            print "Erreur de recuperation dans la base RefSeq  : $url2 !\n" . $page2->error_as_HTML . "\n\n";
    }
	# Si la page a ete recuperee avec succes
    else    {
		# Recuperation du contenu de la page dans $i
        $i = $page2->content ;

		# Recuperation des informations qui nous interessent dans la page
        @Line = split(/LOCUS/,$i);
		if ($#Line >= 0){
			for ($i = 1; $i <= $#Line; $i ++){
				$tmpSpecie = '-';
				$tmpSeq = '';
				if ($Line[$i] =~ /\nSOURCE *([^ ][^\n]*)\n/){
					$tmpSpecie = $1 ; # Recuperation de l'organisme					
				}
				else{
					print "No specie found for $id in RefSeq\n";
				}
				if ($Line[$i] =~ /\nORIGIN([^\/]*)\//){
					$tmpSeq = $1 ; # Recuperation de la sequence
					$tmpSeq =~ s/[\n 0-9]//g; # Suppresion de la position et des retours a la ligne
					if (($tmpSeq =~ /^[ATGCNnatgc]*$/) && ($Line[$i] =~ /\/translation="([ \na-zA-Z]*)"/)){
						$tmpSeq = $1;
						$tmpSeq =~ s/[\n ]//g; # Suppresion des espaces et des retours a la ligne
					}
				}
				#Enregistrement de la sequence
				if (($tmpSeq ne '') && !($tmpSeq =~ /^[ATGCNnatgc]*$/)) {
					if ($seq eq ''){
						$seq = $tmpSeq;
						$specie = $tmpSpecie;
						$$refSpecie{$specie} = searchTaxon($specie) unless (exists($$refSpecie{$specie}));  # Recuperation du taxon
					}
					else{
						$seq = "Several result";
						$specie .= "\nSeveral result";
					}
				}
			}
			if ($seq eq '') {
				print "Warning : No amino acis sequences for $id in RefSeq\n".join("\n",@Line)."\n\n" ;
			}
			elsif($seq eq 'Several result'){
				print "Warning : Several results for $id in RefSeq (none result saved)\n".join("\n",@Line)."\n\n" ;
			}
			else{
				if (exists($$refSpecie{$specie})){
					open OUT,'>>'.$fileOut;
					chomp($id);
					print OUT "$id\t$ec\t$specie\t$$refSpecie{$specie}\t$seq\tMetaCyc\n" ;
					close(OUT);
					if ($specie eq '-'){
						print "Warning : No specie for $id in Ref Seq\n".join("\n",@Line)."\n\n" ;
					}
				}
				else{
					print "Warning : none taxonomy found for $specie : \n$id\t$ec\t$specie\t$seq\tMetaCyc\n" ;
				}
			}
		}
		else{
			print "Warning : No results for $id in RefSeq\n" ;
		}
     }
}

# Extraction de la sequence $id de la base de donnnees PDB
sub Extract_PDB_seq{
		#print STDERR "dans extract pdb seq\n";
    my ($id,$ec,$fileOut,$refSpecie) = @_;
    my @Line = ();
    my %Seq = ();
    my $header = '';
    my $specie = '-';
    my $seq = '';
		my $i = '';
	# Recuperation des informations reliees a $id
    $i = get('http://www.rcsb.org/pdb/explore/explore.do?structureId='.$id) ;
    $specie = $1 if ($i =~ /TreeEntity[^>]*>([^>]*)<\/a>/);
    $specie =~ s/^[^A-Za-z]*//g;
    $specie =~ s/[^A-Za-z]*$//g;
    $$refSpecie{$specie} = searchTaxon($specie) unless (exists($$refSpecie{$specie}));
	# Recuperation de la sequence $id
    $seq = get('http://www.rcsb.org/pdb/files/fasta.txt?structureIdList='.$id) ;
    #$seq du type ">2F9Y:A|PDBID|CHAIN|SEQUENCE
		#MGSSHHHHHHSSGLVPRGSHMSLNFLDFEQPIAELEAKIDSLTAGSRQDEKLDINIDEEVHRLREKSVELTRKIFADLGA
		#WQIAQLARHPQRPYTLDYVRLAFDEFDELAGDRAYADDKAIVGGIARLDGRPVMIIGHQKGRETKEKIRRNFGMPAPEGY
		#RKALRLMQMAERFKMPIITFIDTPGAYPGVGAEERGQSEAIARNLREMSRLGVPVVCTVIGEGGSGGALAIGVGDKVNML
		#QYSTYSVISPEGCASILWKSADKAPLAAEAMGIIRPRLKELKLIDSIIPEPLGGAHRNPEAMAASLKAQLLADLADLDVL
		#STEDLKNRRYQRLMSYGYA
		#>2F9Y:B|PDBID|CHAIN|SEQUENCE
		#MSWIERIKSNITPTRKASIPEGVWTKCDSCGQVLYRAELERNLEVCPKCDHHMRMTARNRLHSLLDEGSLVELGSELEPK
		#DVLKFRDSKKYKDRLASAQKETGEKDALVVMKGTLYGMPVVAAAFEFAFMGGSMGSVVGARFVRAVEQALEDNCPLICFS
		#ASGGARMQEALMSLMQMAKTSAALAKMQERGLPYISVLTDPTMGGVSASFAMLGDLNIAEPKALIGFAGPRVIEQTVREK
		#LPPGFQRSEFLIEKGAIDMIVRRPEMRLKLASILAKLMNLPAPNPEAPREGVVVPPVPDQEPEA"
    if ($seq ne ''){#si sequence trouvée
        @Line = split(/\n/,$seq);
        $seq = '';
				# Pour chaque ligne retournee (peut avoir plusieurs sequences)
				print "scalar ligne=".scalar(@Line)."\n";
        for ($i = 0; $i < scalar(@Line); $i ++){
						# Si nouvelle sequence trouvee
            if ($Line[$i] =~ /^>/){
								# On enregistre la precedente
                $Seq{$seq}{$id} = '' if ($seq ne '');
								# Recuperation du nouvel ID
                $id  = $1 if ($Line[$i] =~ />([^\|]*)\|/);
                $seq = '';
            }
						# Recuperation de la suite de la sequence
            else{
                $seq .= $Line[$i];
            }
        }
				# Enregistrement de la derniere sequence
        $Seq{$seq}{$id} = '' if ($seq ne '');
				# Enregistrement des sequences retournees (peut avoir plusieurs IDs avec meme sequence
        open OUT,'>>'.$fileOut;
				# Pour chaque sequence proteique
        foreach $seq (keys(%Seq)){
						# Tableau des IDs partageant la meme sequence
            @Line = sort(keys(%{$Seq{$seq}})); 
						# On enregistre uniquement le premier ID
						my $prem=shift(@Line);
						chomp($prem);
            print OUT $prem."\t$ec\t$specie\t$$refSpecie{$specie}\t$seq\tMetaCyc\n";
        }
        close(OUT);
    }
    else {
        print "Pb with $id (PDB): No res\n";
    }
   
}

# Recherche des informations taxonomiques reliees a une espece $specie donnee
sub searchTaxon  {
#print STDERR "searchtaxon\n";
    my ($specie) = @_;
    my $res = '';
		my $taxon = '';
		my $idTaxon = '';
    my @Line = ();
		my $tmpSpecie = '';
		my $otherInfo = '';
		if ($specie =~ /^([A-Z][^ ]* [^ ]*) (.*)$/){
			$tmpSpecie = $1;
			$otherInfo = $2;
		}
		elsif ($specie =~ /^[a-z ]* ([A-Z][^ ]* [^ ]*) (.*)$/){
			$tmpSpecie = $1;
			$otherInfo = $2;
			$tmpSpecie =~s/ /+/g;
			$otherInfo =~s/ /+/g;
			$otherInfo =~s/[\(\)]//g;
		}
		else{
			$tmpSpecie = $specie;
			$tmpSpecie =~s/ /+/g;
		}
		$tmpSpecie=~s/\'//g;#ajout cecile
		if ($otherInfo ne ''){
			($res = get('http://www.uniprot.org/taxonomy/?query=%22'.$otherInfo.'%22+AND+(common%3A%22'.$tmpSpecie.'%22+OR+scientific%3A%22'.$tmpSpecie.'%22)&sort=score'))  || die ("Problem with specie $specie\n".'http://www.uniprot.org/taxonomy/?query=%22'.$otherInfo.'%22+AND+(common%3A%22'.$tmpSpecie.'%22+OR+scientific%3A%22'.$tmpSpecie.'%22)&sort=score');
			unless ($res =~ /<p class="summary">/){
				($res = get('http://www.uniprot.org/taxonomy/?query=.'.$tmpSpecie.'+'.$otherInfo.'&sort=score'))  || die ("Problem with specie $specie\n".'http://www.uniprot.org/taxonomy/?query=.'.$tmpSpecie.'+'.$otherInfo.'&sort=score');
			}
		}
		else{
			($res = get('http://www.uniprot.org/taxonomy/?query=(common%3A%22'.$tmpSpecie.'%22+OR+scientific%3A%22'.$tmpSpecie.'%22)&sort=score'))  || die ("Problem with specie $specie\n".'http://www.uniprot.org/taxonomy/?query=(common%3A%22'.$tmpSpecie.'%22+OR+scientific%3A%22'.$tmpSpecie.'%22)&sort=score');
		}
		unless ($res =~ /<p class="summary">/){
			($res = get('http://www.uniprot.org/taxonomy/?query=.'.$tmpSpecie.'&sort=score'))  || die ("Problem with specie $specie\n");
		}
							
    if ($res ne '')    {
		@Line = split(/<tr>/,$res);
		$res = '';
		for (my $i = 0; $i <= $#Line; $i ++){
			# if ($Line[$i] =~ /<p class="summary">([^\<]*)<\/p>.*uniprot\/\?query=organism:([0-9]*)\+/){
				# $taxon = $1;
				# $idTaxon = $2;
			# }
			if ($Line[$i] =~ /<a href=".\/([0-9]*)">([^<]*)<\/a>.*<p class="summary">([^\<]*)<\/p>/){
				$idTaxon = $1;
				$tmpSpecie = $2;
				$taxon = $3;
				$res = "$taxon\t$idTaxon" ;

				unless ($res =~ /related/){
					$res =~ s/&rsaquo;/;/g;
					$res =~ s/ ; /_/g;
					return $res;
				}
			}
			elsif (($res eq '') && ($Line[$i] =~ /<p class="summary">([^\<]*)<\/p>/)){
				$res = "$1\t-" ;
				$res =~ s/&rsaquo;/;/g;
				$res =~ s/ ; /_/g;
			}
		}
		if ($res eq ''){
			print "No taxonomy for $specie\n".'http://www.uniprot.org/taxonomy/?query=%22'.$otherInfo.'%22+AND+(common%3A%22'.$tmpSpecie.'%22+OR+scientific%3A%22'.$tmpSpecie.'%22)&sort=score'."\n";
			print join("\n",@Line)."\n";
		}
		else{
			print "No taxonomy ID for $specie ($res)\n";
			print join("\n",@Line)."\n";
		}
		return "-\t-";
    }
}

__END__

