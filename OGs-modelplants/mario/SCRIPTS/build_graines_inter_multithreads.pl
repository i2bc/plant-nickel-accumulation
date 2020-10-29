#!/usr/bin/perl -w
#use threads;
use Parallel::ForkManager;
#use threads::shared;

#use Thread::Semaphore;

# seed creation
# Write FASTA FILES:
# - seeds
# - sequences without in a seed (one file by species)
# Create profile HMM for each seed
# 16 june 2014
 
my $DOSINT=$ARGV[0];#intersection folder
my $DOS4METH=$ARGV[1];#folder with final groups
my $inter=$ARGV[2];#intersection minimum size
my $nbmeth=$ARGV[3];
my $dirFASTA = $ARGV[4];#modif 07/03
my $dirWork=$ARGV[5];
#my $nbcpu=$ARGV[6];
my $Dscript=$ARGV[7];
my $methode=$nbmeth+1;
my $dirGroupFasta = $dirWork."/fasta_group$methode/"; # OUTPUT : folder with seeds in fasta
my $dirSeqWithoutGroupFasta = $dirWork."/fasta_seq_without_groups$methode/"; # OUTPUT : folder with protein ungrouped
my $DO="$dirWork/ALIGNEMENT$methode/";#FOLDER FOR ALIGNED SEQUENCES
my $DH="$dirWork/HMMSeeds$methode/";
my $tmpinter="$dirWork/tmpinterchoose/";
mkdir($tmpinter);
#my @jobs=();

#my $threadsAllowed=$ARGV[8];
$pm=Parallel::ForkManager->new($ARGV[8]);
#my $sem=Thread::Semaphore->new($ARGV[8]);#number of cpu allowed

# Folder creation
mkdir($dirGroupFasta) unless (-e $dirGroupFasta);
mkdir($dirSeqWithoutGroupFasta) unless (-e $dirSeqWithoutGroupFasta);
mkdir($DO);
mkdir($DH);
my %idvus=();#keys: protein already in a group

#####################################################
#each sequence is associed to an id and a genome
#####################################################
chdir($dirFASTA);
my @fichfasta=glob("*.fa");
my %idseq=();
foreach my $f(@fichfasta){
	my $genome=$f;
	$genome=~s/\.fa$//;
	open(IN,$f);
	my $id="";
	while(<IN>){
		chomp;
		if(/>(.*)/){
			$id=$1;
			my @tmp=split(/\s/,$id);
			$id=$tmp[0];
			$idseq{$id}{$genome}="";#idseq specie sequence
		}
		else{
			$idseq{$id}{$genome}.=$_;
		}
	}
	close(IN);
}


###########################################################################
# READING FINAL GROUPS IN ORDER TO KNOW PROTEINS ALLREADY IN A INTERSECTION
###########################################################################

#number of groups already mades and proteins already in a FINAL GROUP
mkdir($DOS4METH);
chdir($DOS4METH);
my @idGrouptmp=glob("*.fa");
my $maxig=0;#max id of FINAL GROUP
my $ig="";
foreach $ig(@idGrouptmp){
	$ig=~s/\.fa$//;
	if($ig>$maxig){
		$maxig=$ig;
	}
	#DELETE ID already in a final group
	open(IN,"$ig.fa")or die("ERROR: can't open $ig.fa");
	while(<IN>){
	  if(/^>/){
	    chomp;
	    my @tmpline=split(/\s/,$_);
	    my $idtmp=$tmpline[0];
	    $idtmp=~s/^>//;
	    $idvus{$idtmp}='';
	  }
	}
	close(IN);
}
$idGroup=$maxig;
$idGroup=$idGroup+10;

###########################################################################
# READING INTERSECTION OF nbmeth+1 METHODS
###########################################################################
chdir($DOSINT);
%groupes3=();
$nbg=0;
# KEEP GROUPS WITHOUT PROTEINS IN FINAL GROUPS
@files=glob("*");
foreach $f(@files){
  @tmp=split(/%/,$f);
  if($#tmp==$nbmeth){
    open(IN,$f)or die("ERROR: can not open the file $f");
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
		    $groupes3{$nbg}=[@tmp];
		    $nbg++;
	    }
    }
    close(IN);
  }
}
open(OUT,">$tmpinter/grintertmp");
@cleestries=sort{$#{$groupes3{$b}}<=>$#{$groupes3{$a}}}keys(%groupes3);
foreach $c(@cleestries){
	if(($#{$groupes3{$c}}+1)<$inter){
		last;	
	}
	print OUT join("\t",@{$groupes3{$c}})."\n";
}
close(OUT);
#if two intersections of the same size, keep the bigger.
conservergroupes("$tmpinter/grintertmp",($nbmeth+1),$idGroup);

chdir($dirGroupFasta);
@fastamuscle=glob("*\.fa");
foreach $f(@fastamuscle){
	#Fork
	my $pid=$pm->start and next;
	$f=~s/\.fa//;
	system("perl $Dscript/build_graines_inter_multithreads_job.pl $DH $f $dirGroupFasta $DO $dirWork $Dscript");
	$pm->finish;
}
$pm->wait_all_children;
#############
#fonctions
#############

# selection of intersections (intermediate groups)
# alignment multiple of intermediate groups muscle
# creation of profile HMMs
# creation of fasta file of unassigned sequences
sub conservergroupes{
	($figroup,$nbmeth,$nbgc)=@_;
	my @jobs=();
	open(IN,$figroup);
	while(<IN>){
	        chomp;
	        my @groupetmp=split(/\t/,$_);
		$groupegarde=1;
		foreach $id(@groupetmp){
			if(exists($idvus{$id})){
				$groupegarde=0;
				last;
			}
		}
		if($groupegarde){
			open(OUT,">$dirGroupFasta/$nbgc.fa")or die("ERROR: Can not open the file $dirGroupFasta/$nbgc.fa");
			foreach $id (@groupetmp){
				$idvus{$id}="";
				my @tmp=keys(%{$idseq{$id}});
				print OUT ">$id\n$idseq{$id}{$tmp[0]}\n";
				delete($idseq{$id});
			}
			close(OUT);
			#$sem->down;
			#push(@jobs,threads->create("lethread",$Dscript,$DH,$nbgc,$dirGroupFasta,$DO,$dirWork));
			$nbgc++;
			#sleep(1) while(scalar threads->list(threads::running)>=$ta);
			#$count = threads->list(threads::running);#a retirer, juste pour voir combien de threads tournent en meme temps;
			#print "nbthreads : $count\n";
		}
		#delete($$refhash{$cle});
	  # }
	   #delete($$refhash{$cle});
	}
	# CREATION OF UNGROUPED PROTEINS FASTA FILES
	foreach my $idpseul (keys(%idseq)){
		unless(exists($idvus{$idpseul})){
		  my @tmpgenome=keys(%{$idseq{$idpseul}});
		  open(OUT,">>$dirSeqWithoutGroupFasta/$tmpgenome[0].fa") or die ("impossible d'ouvrir le fichier $dirSeqWithoutGroupFasta/$tmpgenome[0].fa");
		  print OUT ">$idpseul\n$idseq{$idpseul}{$tmpgenome[0]}\n";
		  close(OUT);
		}
	}
	close(IN);
	$_->join for @jobs;
	
}

sub lethread{
	my ($Dscript,$DH,$nbgc,$dirGroupFasta,$DO,$dirWork)=@_;
        system("perl $Dscript/build_graines_inter_multithreads_job.pl $DH $nbgc $dirGroupFasta $DO $dirWork $Dscript");
}
