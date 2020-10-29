#!/usr/bin/perl -w
# description
# add groups comming from intersection of 2 methods with size lower than $inter
# cecile pereira
# 06 march 2014

$DOSINT=$ARGV[0];#Intersection folder
$DOSRES=$ARGV[1];#seed folder
$DOSMETH=$ARGV[2];#folder with final groups : FASTA_FINAL
$inter=$ARGV[3];
mkdir($DOSRES);
%idvus=();#keys proteins already in a final group

#saver protein  already in a final group in %idvus
chdir($DOSMETH);
@fichiers4m=glob("*fa");
foreach $f(@fichiers4m){
	open(IN,$f)or die("ERROR: can't open $f");
	while(<IN>){
		chomp;
		@tmp=split(/\t/,$_);
		if($tmp[0]=~/^>(.+)/){
			@tmpidautre=split(/\s/,$1);
			$idtmp=$tmpidautre[0];
			$idvus{$idtmp}='';	
		}
	}
	close(IN);
}

chdir($DOSINT);

#read intersection of 2 methods
%groupes2=();
$nbg=0;
chdir("$DOSINT");
@t=glob("*");
foreach $te(@t){
  if($te=~/^\w+\%\w+$/){
    $inter2="$DOSINT/$te";
    open(IN,$inter2) or die("ERROR: can't open $inter2");
    while(<IN>){
	    chomp;
	    @tmp=split(";",$_);
	    $groupegarde=1;
	    foreach $val(@tmp){
		    if(exists($idvus{$val})){
			    $groupegarde=0;
			    last;
		    }
	    }
	    if($groupegarde){
		    $groupes2{$nbg}=[@tmp];
		    $nbg++;
	    }
    }
    close(IN);
  }
}

#if two groups with the same proteins, only keep the larger
conservergroupes(\%groupes2,2,$inter);

#############
#fonctions
#############
sub conservergroupes{
	($refhash,$nbmeth,$mininter)=@_;
	@cleestries=sort{$#{$$refhash{$b}}<=>$#{$$refhash{$a}}}keys(%$refhash);
	foreach $cle(@cleestries){
	   $taillegroupe=$#{$$refhash{$cle}}+1;
	   if($taillegroupe<$mininter){
		$groupegarde=1;
		foreach $id(@{$$refhash{$cle}}){
			if(exists($idvus{$id})){
				$groupegarde=0;
				last;
			}
		}
		if($groupegarde){
			open(OUT,">$DOSRES/$$refhash{$cle}[0].out")or die("ERROR: can't open the file $DOSRES/$$refhash{$cle}[0].out");
			$groupRef=join(",",@{$$refhash{$cle}});
			print OUT "SELECTED 2\t$groupRef\t$nbmeth\n";
			close(OUT);
			foreach $id(@{$$refhash{$cle}}){
				$idvus{$id}="";
			}
		}
		delete($$refhash{$cle});
	   }
	}
}
