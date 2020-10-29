package tab2mario;
#fonction permettant de parser les fichiers d'entree contenant les groupes d'orthologues et les inparalogues
#entrees:
#le dossier où sont stockés les groupes
#l'extension du fichier contenant les orthologues
#l'extension du fichier contenant les paralogues
#la methode pour lequelles on lit les fichiers resultats
#sorties:
#hash ayant pour clé les groupes faits par la méthode
sub convert
{  
  my ($meth, $dirGroup,$DOSSORT)=@_;
  my $fileOrthologs = $meth.'_orthologous_groups.sql';
  my $fileInparalogs= $meth.'_inparalogs.sql';
  my %Group=(); 
  # Lecture des inparalogues
  print $dirGroup.$fileInparalogs."\n";
  open IN, $dirGroup.$fileInparalogs or die ("imposssible douvrir le fichier inparalogs");
  while (<IN>) {#lecture du fichier de paralogue
    chomp($_);
    @Res = split(/\s/,$_); # colonne 1 : ID de la proteine, colonne 2 : numero du groupe#modif cecile avt split("\t",$_)
    if (exists($Group{$Res[1]})){
      $Group{$Res[1]}.= ';'.$Res[0];		
#	print "$Res[1] $Group{$Res[1]}\n";				
    }
    else{
      $Group{$Res[1]}=$Res[0];
    }
  }
  close(IN);
  # Lecture des groupes d'orthologues
  #my %groupes=();
  open IN, $dirGroup.$fileOrthologs or die("impossible d'ouvrir $dirGroup$fileOrthologs");
  print "$DOSSORT/$meth\_groups\n";
  open OUT , ">".$DOSSORT."/".$meth."_groups" or die ("impossible d'ouvrir le ficher $DOSSORT/$meth\_groups");
  while (<IN>) {
  #exemple de ligne du fichier:canbrvCABR.00029-j161	cancavCACA.00015-j229	canglgnlGLVCAGL0G06666g	cannivCANI.00004-j76	klubavKLBA.00010-j91	kludevKLDE.00012-j221	sacceSCRT_03239	1
    chomp($_);
    # Recuperation des differentes colonnes dans un tableau
    @Res = split(/\s/,$_);#modif cecile avant split("\t",$_)
    $nb = pop(@Res); # numero du groupe
    print "$nb\n";
    $groupetmp=join(";",@Res);
    $groupetmp=~s/\s/\;/g;
    $groupetmp=~s/;-//g;
    $groupetmp=~s/^-;//;
    if(exists($Group{$nb})){
    	print OUT $groupetmp.";".$Group{$nb}."\n";
    }
    else{
    	print OUT $groupetmp."\n";
    }
  }
  close(IN);
  close(OUT);
}

1;
