#!/usr/bin/perl -w
use strict;

# Description
# HMM database creation
# perl database_creation.pl [Result folder] [methods]
# cecile pereira, 2/05/2012
#######################################################################
# Declaration
###########

# Parametres
my $dirWork = $ARGV[0]; # INPUT : OUTPUT FOLDER
my $methode=$ARGV[1];
my $dscript=$ARGV[2];
my $dirHmmDb = $dirWork."/BD$methode/"; #OUTPUT FOLDER OF THE HMM DATABASE
my $dirHmm = $dirWork."/HMMSeeds$methode/"; # INPUT : FOLDER WITH HMM PROFILS

# Autres variables
my @File = ();
my $file = '';

#######################################################################
# Programme
###########

# FOLDER CREATION
mkdir ($dirHmmDb) unless (-e $dirHmmDb);

# COPY all hmm profils on the same file
chdir($dirHmm);
@File = glob("*msf");
foreach $file (@File) {
	#print "$file\n";
    	system('cat '.$file.' >> '.$dirHmmDb.'hmm_all');
}

#HMM database creation
chdir($dirHmmDb);
system("hmmpress hmm_all >> $dirWork/ERR/hmmpress$methode.err");
