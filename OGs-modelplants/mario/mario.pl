#!/usr/bin/perl
use warnings;
use Cwd;
use Math::Combinatorics;
use Getopt::Long;
use Parallel::ForkManager;

# Cecile PEREIRA
# 25/06/2014
# Meta-approach of homologous groups detection
# you can find a complete description in the paper : [nom du papier]

#------------------------------------------
# PARAMETERS
#------------------------------------------
#default parameters
$RepCourant = cwd();
%arguments=();
$arguments{'inter'}=4;
$arguments{'e'}=1e-10;
$arguments{'pa'}=40;
$arguments{'cpu'}=1;
my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
my $dateUs =  ($mon+1)."_".$mday."_".(1900+$year);
$arguments{'res'}="$RepCourant/RESULTS_metOG_$dateUs/";
$arguments{'tmp'}="$RepCourant/TMP_METAPP_metOG_$dateUs/";
$version=1.1;

#Save parameters and test parameters format
GetOptions(\%arguments,"f=s","i=s","inter:i","e:f","pa:f","cpu:i","res:s","tmp:s");

#print help
if(defined $arguments{'h'}){
  help_parameters();
  die;
}

#test input parameters
if((! defined $arguments{'i'}) or (! defined $arguments{'f'})){
  print "DIE: the folder containing result of initial methods and the result folder are needed\n";
  help_parameters();
  die;
}

#-------------------------------------------
#print start date and parameters
print "Meta-approach, version $version\n";
print "Started :\t";
system("date");
print "Parameters : \n";
foreach $t(keys(%arguments)){
  print "$t\t$arguments{$t}\n";
}
print "\n";

$DScript="$RepCourant/SCRIPTS/";
$RES=$arguments{'res'};
$TMP=$arguments{'tmp'};
print "Tempory files : $TMP\nResults files : $RES\n";

#------------------------------------------
# INPUT ARE TESTED
#------------------------------------------

#test the number of input files (initial methods)
chdir($arguments{'i'});
@homo_init=glob("*");
$nb_meth=@homo_init;
if($nb_meth<=1){
  print "At least the result of two methods have to be compared\n";
  help_parameters();
  die;
}

# convert input file if necessary
# creation of the output folder
mkdir($RES);
mkdir($TMP);
mkdir("$TMP/ERR/");
mkdir("$TMP/OUT/");
mkdir("$TMP/Intersections");
system("cp * $TMP/Intersections/") and die("ERROR: can not copy input files");
chdir("$TMP/Intersections/");
# parser for orthoXML format:
foreach $hi (@homo_init){
	if($hi=~/\.xml$/){ 
	  system("perl $DScript/parser_orthoXML.pl $hi") and die("ERROR in orthoXML parser or orthoXML files");
	  system("rm $hi");
	}
}

#convert seqXML to fasta format
chdir($arguments{'f'});
@filesxml=glob('*.xml');
$nbfxml=@filesxml;
if($nbfxml>0){
  mkdir("$TMP/FASTA/");
  foreach $tmpf(@filesxml){
    system("perl seqxml2any --format fasta $tmpf") and die("ERROR in conversion of seqxml: seqxml2any");
    $nom=$tmpf;
    $nom=~s/xml//;
    system("mv $tmpf\.fasta $TMP/FASTA/$tmpf\.fa") and die("ERROR: mv .fa to .fasta");
  }
  system("cp *.fa $TMP/FASTA/") and die("ERROR copy .fa to temporary FASTA folder");#copy of input files in fasta format (.fa)
  $arguments{'f'}="$TMP/FASTA/";
}

#-------------------------------------------
# MAIN
#-------------------------------------------

# 1) creation of all possible intersections
for($combination=2;$combination<=$nb_meth;$combination++){
  @tab=combine($combination,@homo_init);#tableau de tableau
  foreach $t(@tab){
    #do the intersection of the methods in the table $t
    @methodsinter=sort @$t;#this methods have to be in the intersection
    $mfinal=pop(@methodsinter);
    $mstart=join("%",@methodsinter);
    unless(-e "$TMP/Intersections/$mstart\%$mfinal"){
      system("perl $DScript/intersection_2_methods.pl $TMP/Intersections/$mstart $TMP/Intersections/$mfinal > $TMP/Intersections/$mstart\%$mfinal") and die("ERROR script intersection");
    }
  }
}

# 2) Tests if intersections are not empty
chdir("$TMP/Intersections/");
@intertmp=glob("*%*");
$nb=0;
foreach $it(@intertmp){
	if(-z $it){
		print "WARNING: interesection file $it empty\n";
		$nb++;
	}
}
if($nb==($#intertmp+1)){
	die("NO intersections");
}


# 3) Name of species
chdir($arguments{'f'});
@species=glob("*.fa");

# 4) Seeds and HMM
print "$arguments{'pa'}";
$pm=Parallel::ForkManager->new($arguments{'cpu'});
for($i=($nb_meth-1);$i>=1;$i--){
  $nbinter=$i+1;
  print "\nIntersection of $nbinter methods\ndate: ";
  system("date");
  # seed creation, seed fasta, profile HMM for each seed and one file by specie for ungrouped sequences
  system("perl -I $DScript/perl_library/ $DScript/build_graines_inter_multithreads.pl $TMP/Intersections/ $TMP/FINAL_GROUPS/ $arguments{'inter'} $i $arguments{'f'} $TMP $arguments{'cpu'} $DScript $arguments{'cpu'}") and die("ERROR creation intermediate groups, Muscle and HMM and unassigned sequences. Step : intersection of  $nbinter methods");
  # HMM profile database creation
  system("perl $DScript/database_creation.pl $TMP $nbinter $DScript") and die("ERROR database creation step $nbinter");
  # Foreach specie, comparison ungrouped sequences and HMM profile and analyse of the result
  mkdir("$TMP/hmmscan_out$nbinter");
  mkdir("$TMP/hmmscan$nbinter/");
  mkdir("$TMP/hmmanalysis$nbinter/");
  #@jobs=();
  foreach $s(@species){
	my $pid=$pm->start and next;
    	mkdir("$TMP/hmmscan$nbinter/$s/");
	system("perl $DScript/jobhmmscan_hmmanalyse.pl $TMP $nbinter $s $DScript $arguments{'e'} $arguments{'pa'}");
	$pm->finish;
  }
  $pm->wait_all_children;
  system("perl $DScript/fastafinalgroups.pl $RepCourant $nbinter $arguments{'f'} $arguments{'e'} $arguments{'pa'} $TMP") and die("ERROR final group creation for the step intersection of $nbinter methods");
  system("cp $TMP/group$nbinter/fastagroup$nbinter/* $TMP/FINAL_GROUPS/") and die ("ERROR copy group$nbinter/fastagroup$nbinter to FINAL_GROUPS");
  system("rm $TMP/BD$nbinter -rf") and die("ERROR removing BD$nbinter to temporary files");
  system("rm $TMP/HMMSeeds$nbinter -rf") and die("ERROR removing HMMSeeds$nbinter to temporary files");
}

print "Add intersection of 2 methods with size < minimum size threshold\ndate\n";
system("date");

# 5) Add intersection of 2 methods with size < minimum size threshold and no sequences already presents in final groups
system("perl $DScript/select_groups_size_inf_minsize.pl $TMP/Intersections/ $TMP/Build_seeds_inter2T/ $TMP/FINAL_GROUPS/ $arguments{'inter'}") and die("ERROR in the add of groups with size < minimum size parameter");
system("perl $DScript/fasta_groups_size_inf_minsize.pl $arguments{'f'} $TMP/Build_seeds_inter2T/ $TMP/FINAL_GROUPS/") and die('ERROR in the add of groups with size < minimum size parameter');

# 6) final results in RES folder 
system("tar -cvjf $RES/FINAL_GROUPS.bz2 $TMP/FINAL_GROUPS/");
system("perl $DScript/saveOrthoXML.pl $TMP/FINAL_GROUPS/ $RES/meta-approach.xml") and die('ERROR: results just in fasta format, not produce in orthoXML format');

print "finished:\t";
system("date");

#----------
#functions
#----------
sub help_parameters{
  my $text=<<"TEX";

Meta-approach script version $version
---------------------------------------------------------------
Fallowing library and programs have to be installed:
	perl library (use cpan for the installation): 	Math::Combinatorics, 
							Cwd,  
							Bio::SeqIO,  
							File::Basename, 
							File::Spec, 
							XML::Simple 
							Parallel::ForkManager;
	programs : muscle, HMMER in SCRIPTS folder
---------------------------------------------------------------

Parameters :
-h	print the help
-i	folder where are stored homolog groups to combine
		one file by initial method (put at least the result of 2 methods)
		two files format allowed
		1) orthoXML format : description at http://orthoxml.org/
		   file name have to be ended by '.xml'
		   Noted that if the file contain paires and not groups, groups will be made as 
		   ensemble of proteins with a relation of orthology with all other proteins in the same group
		2) "groups" format : one line by group
				     protein names separed by a ';'
-f	folder containing proteomes sequences files both fasta files or seqXML files are allowed, 
		warning: use one fasta file by specie
		sequences files in seqXML: file name have to be ended by '.xml'
		sequences files in fasta format: fine name have to be ended by '.fa'
-res    result folder, default parameter : RESULTS_metOG_date/
-tmp    tempory folder, default parameter : TMP_metOG_date/
-inter	minimum size of the selected intersections for the creation of intermediate groups
		default parameter : 4
-e	e-value threshold, default parameter : 1e-10
-pa	aligment percentage threshold, default parameter : 40 
-cpu	number of used cpu by hmmer, default parameter : 1

example: 
perl MaRiO.pl -i ~/OrthologGroups/ -f ~/fastaProteomes/ -cpu 2
TEX
  print "$text\n";
}
