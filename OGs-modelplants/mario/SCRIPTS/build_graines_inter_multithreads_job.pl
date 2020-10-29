


$DH=$ARGV[0];
$nbgc=$ARGV[1];
$dirGroupFasta=$ARGV[2];
$DO=$ARGV[3];
$dirWork=$ARGV[4];
$Dscript=$ARGV[5];





$nbtest=0;

while((! -e "$DH/$nbgc.msf")&&($nbtest<5)){
	system("muscle -in $dirGroupFasta$nbgc.fa -quiet -out $DO/$nbgc.fa -loga $dirWork/ERR/muscle$nbgc.err");
        if(-e "$DO$nbgc.fa"){
        	#fasta to stockholm format
                system("$Dscript/fasta2stockholm.pl $DO/$nbgc.fa > $DO/$nbgc.sto");
                #BUILD HMM PROFILE
                #system("hmmbuild --cpu $nbcpu $DH/$nbgc.msf $DO/$nbgc.sto > $dirWork/OUT/hmmbuild$nbg.out >> $dirWork/ERR/hmmbuild$nbgc.err");
                system("hmmbuild --cpu 1 $DH/$nbgc.msf $DO/$nbgc.sto > $dirWork/OUT/hmmbuild$nbgc.out >> $dirWork/ERR/hmmbuild$nbgc.err");
        }
        $nbtest++;
}
if(! -e "$DH/$nbgc.msf"){
        print STDERR "Problem the intersection $nbgc not tacking into account\n";
}

