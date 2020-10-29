package Fonctions;
use strict;

# Construit un Hash qui fait le lien entre les pids "Vrai"
# et les pids temporaires:
# Clés: Pids "Vrai"
# Valeurs: Pids "temporaires".
# Les données sont obtenues du fichier en argument
sub lirePidTemp{
    my ($file,$log) = @_;
    my %hash;

    open(LOG,">>$log");
    print  "Lecture des Pids Temporaires: $file\n";
    print LOG "Lecture des Pids Temporaires: $file\n";
    open(IN,$file) or die "Fichier Inconnu";
    while(<IN>){
	chomp;
	my ($a,$b) = split(/\t/);
	$hash{$a} = $b;
    }
    close(IN);
    close(LOG);

    return %hash;
}

# Cette fonction est très semblable à la précédente
# Sauf qu'elle renvoie un hash avec comme:
# Clés: Pids temporaires
# Valeurs: Pid "réèls"
sub lirePid{
    my ($file,$log) = @_;
    my %hash;

    open(LOG,">>$log");
    print  "Lecture des Pids Temporaires: $file\n";
    print LOG "Lecture des Pids Temporaires: $file\n";
    open(IN,$file) or die "Fichier Inconnu";
    while(<IN>){
	chomp;
	my ($a,$b) = split(/\t/);
	$hash{$b} = $a;
    }
    close(IN);
    close(LOG);

    return %hash;
}

# Ecris le fichier en argument 1
# Contenant: le lien entre les Pids "Vrai" et les pids temporaires
sub ecrirePidTemp{
    # Repertoire en argument 1
    # Fichier de sortie en argument 2
    # et Fichier de Log en argument 3
    my ($dir,$pidOut,$logFile) = @_;
    # Liste des fichiers du repertoire en argument
    my @files;
    # Abréviation des espèces
    my $abrv;
    my $chr;
    # Pid temporaire courant
    my $pid_temp=0;
    
    opendir(DIR,$dir);
    my @temp = readdir(DIR);
    foreach(@temp)
    {
	push @files,$_ if(/\.sgml/);
    }
    closedir(DIR);
    
    open(LOG,">>$logFile");
    open(OUT,">pidtemp");
    print "Ecriture $logFile\n";
    print LOG "Ecriture $logFile\n";
    foreach my $f (@files)
    {
	print "$f\n";
	open(IN,$ARGV[0]."/".$f);
	# tableau contenant les pids reels des espèces comme clé
	# Et comme valeur la position de début du pid dans le génome
	my %pid;
	while(<IN>)
	{
	    # Si c'est une ligne de description de l'espèce
	    if(/\A\<DESCR\>/)
	    {
		
		# On récupère l'abbréviation de l'espèce
		my ($abrv) = /\<ABRV\>.*?\<\/ABRV\>/;
		# On récupère le chromosome en cours, et le cas échéant le plasmide
		if(/\<CHR\>(.*?)\<\/CHR\>/)
		{
		    $chr = $1;
		}
		else
		{
		    ($chr) = /\<PLAS\>(.*?)\<\/PLAS\>/;
		}
	    }
	    # Sinon, on récupère les pids
	    else
	    {
		# PID
		my ($pid) = /\<PID\>(.*?)\<\/PID\>/;
		
		# Start
		my($start)=/\<START\>(.*?)\<\/START\>/;
		my ($end) = /\<END\>(.*?)\<\/END\>/;
		# On remplie la valeur start pour le pid "pid"
		$pid{$chr}{$pid}=min($start,$end);
	    }
	}
	close(IN);
	
	# On trie maintenant les clés du hash selon la valeur de début!
	foreach my $sAbrv(keys %pid)
	{
	    my %hashTemp = %{$pid{$sAbrv}};
	    my @cles = sort{$hashTemp{$a}<=>$hashTemp{$b}} keys %hashTemp;
	    foreach(@cles)
	    {
		print OUT "$_\t$pid_temp\n";
		$pid_temp++;
	    }
	    $pid_temp+=1000;
	}
    }
    close(OUT);
    close(LOG);
}


sub min{
    my ($a,$b)=@_;
    return $a if($a<$b);
    return $b;
}
1;
