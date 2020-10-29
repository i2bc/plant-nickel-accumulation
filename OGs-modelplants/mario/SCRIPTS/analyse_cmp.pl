#!/usr/bin/perl -w
use strict;

# Description
# Analyse of the comparison result of intermediate group/unassigned sequences

#cecile pereira
#04 march 2014
#######################################################################
# Declaration
###########

#print OUT $dirWork.'6_intersection_jobs.pl '.$dirOutHmmScan."/$dirGeno/ $evalSeuil $dirOutHmmScanout/$dirGeno.out_o $pcalign > $dirOutHmmAnalysis/$dirGeno.out\n";
	
# Parametres 
die ("ERROR : the script need 5 inputs\n") unless ($#ARGV == 4);
my $dirIn = $ARGV[0]; # Repertory with result files  
# INPUT : repertoire contenant les fichiers generes par hmmer lors de la comparaison
my $evalSeuil = $ARGV[1]; # evalue parameter
my $OutHmmScanout=$ARGV[2]; # out_o file
my $pcalign=$ARGV[3];# alignment parameter 
my $HMM_all=$ARGV[4];# record of HMM profils

# Autres parametres
my @File = ();
my @Line = ();
my $file = '';
my $geno = '';
my $dir = '';
my %dejavu = ();
my @tmp=();


# Recuperation de la date
 my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
$year = 1900 + $year;
$mon ++;
$hour = '0'.$hour if ($hour =~ /^[0-9]$/);
$mday = '0'.$mday if ($mday =~ /^[0-9]$/);
$mon = '0'.$mon if ($mon =~ /^[0-9]$/);
my $date =  $year.$mday.$mon;

#######################################################################
# Script
###########

if ($dirIn =~ /^(.*\/)([^\/]*)\/$/){
    $dir = $1; # record with all results
    $geno = $2; # genome name
}
else{
	die ("ERROR with $dirIn\n");
}

chdir($dirIn);
@File = glob("*out");# files create by hmmer for the genome $geno
%dejavu=();#sequence;groupes;evalue

foreach $file (@File) {
#read result for sequences present in $file (all unassigned proteins of the genome $geno)
    open IN,$dirIn.$file or die("ERROR : can't open the file $dirIn$file");
    while (<IN>) {
        unless ($_ =~ /^#/) {
            chomp($_);
            @Line = split(/ +/,$_);
		#$Line[7]: evalue parameter on the best domain
		if($Line[7] < $evalSeuil){
			$dejavu{"$Line[0]\t$Line[2]"}=$Line[7];#groupe \t prot = evalue BEST DOMAINE
		}
        }
    }
    close(IN);
}

#-------------------------------------
# percentage alignment
#-------------------------------------
my %HMMlength=();#profile HMMs length, key number of the profile, value length of the profil
my %querylength=();
my %querymodellength=();
my %querymodelpourcent=();
my %querymodellengthHMM=();
my %querymodelpourcentHMM=();
my $groupe="";
my $query='';
my $tmplen=0;
my $pourcent=0;
my $long=0;
my $longhmm=0;
my $pourcentHMM=0;

#1) association intermediate groupe / profile length
open(HMM,$HMM_all)or die("ERROR can't read hmm_all: $HMM_all\n");
my $namehmm="";
my $lengthhmm=0;
while(<HMM>){
  if(/^NAME/){
    chomp;
    @tmp=split(/\s+/,$_);
    $namehmm=$tmp[1];
  }
  else{
    if(/^LENG/){
      @tmp=split(/\s+/,$_);
      $lengthhmm=$tmp[1];
      $HMMlength{$namehmm}=$lengthhmm;
    }
  }
}
close(HMM);


#2) read result comparison
open(IN,$OutHmmScanout) or die("ERROR can't read $OutHmmScanout");
while(<IN>){
    unless(/^#/){#commentary
      if(/^Query:\s+(.+)\s+\[L=(\d+)\]/){
		$query=$1;
		$tmplen=$2;#sequence length
		$query=~s/\s//g;
		$groupe='';
		$querylength{$query}=$tmplen;#sequence id / sequence length
      }
      else{
	if(/>> (\d+)/){
	    unless($groupe eq ""){#unless first group
			if(exists($dejavu{"$groupe\t$query"})){
				if(exists($querymodellength{$query}{$groupe})){#sequence 
					if(exists($querymodellengthHMM{$query}{$groupe})){
						$pourcent=$querymodellength{$query}{$groupe}{$dejavu{"$groupe\t$query"}}/$querylength{$query}*100;
						$querymodelpourcent{$query}{$groupe}=$pourcent;
						$pourcentHMM=$querymodellengthHMM{$query}{$groupe}{$dejavu{"$groupe\t$query"}}/$HMMlength{$groupe}*100;#ajout deja vu 17/12/2013
						$querymodelpourcentHMM{$query}{$groupe}=$pourcentHMM;#result for the BEST DOMAIN!
					}
				}
			}
	    }

	    $groupe=$1;
	}
	if(/!/){
	
#example of result line :
# >> 3686
#    #    score  bias  c-Evalue  i-Evalue hmmfrom  hmm to    alifrom  ali to    envfrom  env to     acc
#  ---   ------ ----- --------- --------- ------- -------    ------- -------    ------- -------    ----
#    1 !  167.6   0.2   1.7e-53   7.6e-50      72     240 ..      14     182 ..       4     183 .] 0.95

	  #keep % alignment of the best domain

	  #alignment length on sequence
	  my @tmp=split(/\s+/,$_);
	  $long=$tmp[11]-$tmp[10]+1;
	  $querymodellength{$query}{$groupe}{$tmp[6]}=$long;#query group evalue sequence alignment
	  #alignment length on profile
	  $longhmm=$tmp[8]-$tmp[7]+1;
	  $querymodellengthHMM{$query}{$groupe}{$tmp[6]}=$longhmm;
	}
	if(/^\/\//){
	  unless($groupe eq ""){
	    if(exists($dejavu{"$groupe\t$query"})){#evalue ok
	      if(exists($querymodellength{$query}{$groupe}{$dejavu{"$groupe\t$query"}})){
			if(exists($querymodellengthHMM{$query}{$groupe}{$dejavu{"$groupe\t$query"}})){
				$pourcent=$querymodellength{$query}{$groupe}{$dejavu{"$groupe\t$query"}}/$querylength{$query}*100;
				$querymodelpourcent{$query}{$groupe}=$pourcent;
				$pourcentHMM=$querymodellengthHMM{$query}{$groupe}{$dejavu{"$groupe\t$query"}}/$HMMlength{$groupe}*100;
				$querymodelpourcentHMM{$query}{$groupe}=$pourcentHMM;
			}
	      }
	    }
	  }
	}
      }
    }
}
close(IN);

#evalue sort
my @dejavutab=sort{$dejavu{$a}<=>$dejavu{$b}} keys(%dejavu);#lower first

my %dejares=();#hash keys proteins with result already find
foreach my $dv(@dejavutab){
  $dv=~s/ +//g;
  my @tmp=split(/\t/,$dv);
  if(exists($querymodelpourcent{$tmp[1]}{$tmp[0]})){#proteine groupe 
    if(exists($querymodelpourcentHMM{$tmp[1]}{$tmp[0]})){
      if($querymodelpourcent{$tmp[1]}{$tmp[0]}>=$pcalign){
	if($querymodelpourcentHMM{$tmp[1]}{$tmp[0]}>=$pcalign){
	  unless(exists($dejares{$tmp[1]})){
		  print "$tmp[1]\t$tmp[0]\t$dejavu{$dv}\t$querymodelpourcent{$tmp[1]}{$tmp[0]}\t$querymodelpourcentHMM{$tmp[1]}{$tmp[0]}\n";
		  $dejares{$tmp[1]}="";
	  }
	}
      }
    }
  }
}
