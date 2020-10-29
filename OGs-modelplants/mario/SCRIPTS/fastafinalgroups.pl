#!/usr/bin/perl -w
use strict;

# Description
# ortholog groups fasta files
# cecile pereira, 06 march 2013
#######################################################################
# Declaration
###########

# Parametres
my $dir=$ARGV[0];#INPUT initial folder
my $methode=$ARGV[1]; # INPUT number of methods in the intersection
my $genomes=$ARGV[2]; # INPUT : folder with proteomes in fasta format
my $evalSeuil = $ARGV[3]; # INPUT : E-value thershold
my $pcalign = $ARGV[4];# INPUT : alignment thershold
my $dirOut1=$ARGV[5];# INPUT : temporary output folder
my $dirFastaGroup = "$dirOut1/fasta_group$methode/"; # INPUT : forder with intermediate groups in fasta format
my $dirHmmAnlaysis = "$dirOut1/hmmanalysis$methode/"; # INPUT : forder with result keeped (comparison seq/profile HMM validate)
my $dirOut = $dirOut1.'/group'."$methode/"; # OUTPUT : folder with ortholog group resulting of intermediate groups for $methode methods taking into account for the intersection
my $dirOutNewFasta = $dirOut.'fastagroup'."$methode/"; # OUTPUT : fasta 
my $fileOutIdGroupeInter4 = $dirOut."/inter_et_seq_id_seq_$methode.sql";# id associed with the group

# Autres variables
my $file = '';
my $nb = 0;
my $geno = '';
my $id = '';
my %IdInfo = ();
my %Group = ();
my %IdGroup = ();
my @Geno = ();
my @Line = ();
my @File = ();
my @IdGroup = ();

#######################################################################
# Programme
###########

# genome name and sequences
chdir($genomes);
my $idtmp="";
@Geno=glob("*.fa");
for(my $i=0;$i<=$#Geno;$i++){
  open(IN,$Geno[$i]) or die("ERROR can't open $Geno[$i]");
  $Geno[$i]=~s/\.fa$//;
  while(<IN>){
    chomp;
    if(/^>(.+)/){
      $idtmp=$1;
      my @tmpidtmp=split(/\s/,$idtmp);
      $idtmp=$tmpidtmp[0];
      $IdInfo{$idtmp}{'genome'}=$Geno[$i];
    }
    else{
      $IdInfo{$idtmp}{'aa'}.=$_;
      $IdInfo{$idtmp}{'aa'} =~  s/\*$//;
    }
  }
  close(IN);
}

# Verification que les repertoires INPUT existent
die("ERROR: No INPUT directory : dir $dir\n") unless (-e $dir);
die("ERROR: No INPUT directory : dirFastaGroup $dirFastaGroup\n") unless (-e $dirFastaGroup);
die("ERROR: No INPUT directory : dirHmmAnalaysis $dirHmmAnlaysis\n") unless (-e $dirHmmAnlaysis);

foreach $geno (@Geno){
	die("ERROR: No result for $geno ($dirHmmAnlaysis)\n") unless (-e $dirHmmAnlaysis.$geno.'.fa.out');
}

# OUTPUT
mkdir ($dirOut) unless (-e $dirOut);
mkdir ($dirOutNewFasta) unless (-e $dirOutNewFasta);

# foreach intermediate group take each sequences 
chdir($dirFastaGroup);
@File = glob("*fa");
foreach $file (@File) {
    $nb = $1 if ($file =~ /^(.*)\.fa/);
    open IN,$dirFastaGroup.$file or die ("ERROR: can't open  $file");
    while(<IN>) {
        chomp($_);
        if ($_ =~ /^\>(.*)$/) {
	  my $nomid=$1; 
          $IdGroup{$nomid} = '';
          $Group{$nb}{$nomid} = '';
        }
    }
    close(IN);
}

# read result and put sequences into the group
chdir($dirHmmAnlaysis);
@File = glob("*out");
foreach $file (@File) {
    if ($file =~ /^(.*)\.out/) {
        $geno = $1;
        open IN,$dirHmmAnlaysis.$file;
        while (<IN>) {
            chomp($_);
            @Line = split(/\t/,$_);
            unless($Line[1]=~/\|/){
	      $Group{$Line[1]}{$Line[0]} = '';
	    }
        }
        close(IN);
    }
}

# Group creation

# OUTPUT files
open IDSEQ,'>'.$fileOutIdGroupeInter4 or die("ERROR: impossible d'ouvrir $fileOutIdGroupeInter4");
# Foreach group 
foreach $nb (keys(%Group)) {
	# Initialisation of variables for the nb group
		
	@IdGroup = keys(%{$Group{$nb}});
	# Write fasta file of the group
	open FASTA,'>'.$dirOutNewFasta.$nb.'.fa';
	foreach $id (@IdGroup) {
	  if (exists($IdInfo{$id}{'aa'})) {
            print FASTA ">$id\t$IdInfo{$id}{'genome'}\n$IdInfo{$id}{'aa'}\n";
            print IDSEQ "$id\t$nb\n";
	  }
	  else {
            die("ERROR: No amino sequences for $id ($nb)\n");
	  }
      }
      close(FASTA);
}
close(IDSEQ);

# tar fasta file
if ($dirOutNewFasta =~ /\/([^\/]*)\/$/){
	$dirOutNewFasta = $1 ;
	chdir($dirOut);
	system("tar -cjf ".$dirOutNewFasta.".bz2 $dirOutNewFasta/");
}
