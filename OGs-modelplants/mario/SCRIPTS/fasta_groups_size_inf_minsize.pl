#!/usr/bin/perl -w
use strict;

# Description
# Write fasta file :
# - intersection of 2 methods with size < size min
# 
# cecile pereira
# 06 march 2014
#######################################################################
# Parameters
###########

my $dirFASTA = $ARGV[0];
my $dirInGroup=$ARGV[1];
my $RES=$ARGV[2];

# Autres variables
my @File = ();
my @Line = ();
my @Id = ();
my $file = '';
my $id = '';
my $idGroup = 0;
my $genome = '';
#######################################################################
# MAIN
###########

################################################
# For protein, save id and seq
################################################
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
			my @tmpidautre=split(/\s/,$1);
			$id=$tmpidautre[0];
			$idseq{$id}{$genome}="";#idseq specie seq
		}
		else{
			$idseq{$id}{$genome}.=$_;
		}
	}
	close(IN);
}

# id of last final group
chdir($RES);
my @idGrouptmp=glob("*.fa");
my $maxig=0;
my $ig="";
foreach $ig(@idGrouptmp){
        $ig=~s/\.fa$//;
        if($ig>$maxig){
                $maxig=$ig;
        }
}
$idGroup=$maxig;
$idGroup=$idGroup+10;# number step between groups obtain by intersection of several number of methods

# read seeds files
chdir($dirInGroup);
@File = glob("*out");
foreach $file (@File) {
    open (IN,$dirInGroup.'/'.$file) or die("ERROR : can't open the file ".$dirInGroup.'/'.$file);
    while (<IN>) {#lignes du type: SELECTED 0	canbrvCABR.00004G_202,cancavCACA.00033-j268,canglgnlGLVCAGL0I00308g,cannivCANI.00020-j150,klubavKLBA.00037-j32,kludevKLDE.00015G_3753,sacceSCRT_01122	4
            chomp($_);
            #2nd col, the list of proteins on groups (id separed by ',')
            @Line = split("\t",$_);
            $idGroup ++;
	    # Create fasta file
	    @Id = split(/,/,$Line[1]); 
	    open (OUT,'>'.$RES.$idGroup.'.fa') or die("ERROR : can't open the file $RES/$idGroup.fa\n");
	    # foreach protein, take sequence in the hash table and write it in the fasta file
	    foreach $id (@Id) {
		    my @tmp=keys(%{$idseq{$id}});
		    print OUT ">$id\t$tmp[0]\n".$idseq{$id}{$tmp[0]}."\n";
		    delete($idseq{$id});
	    }
	    close(OUT);
    }
    close(IN);
    die ("ERROR file $RES$idGroup.fa empty\n") if(-z $RES.$idGroup.'.fa');
}
