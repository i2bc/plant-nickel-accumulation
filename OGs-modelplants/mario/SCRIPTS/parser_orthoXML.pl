#!/usr/bin/perl
use XML::Simple;
use re 'eval';
#AIM: read orthoXML format and save like : 
#one group by clique of orthologs, 
#one line by group
#sequences separed by a ';'
#cecile pereira



my $FichierXML = $ARGV[0];
my $FichierResulat = $FichierXML;
$FichierResulat=~s/\.xml$/_groups/;

my $parser = XML::Simple->new( KeepRoot => 1 );
my $doc = $parser->XMLin($FichierXML);# Tout le fichier XML est dans $doc sous forme d'arbre


%proteines=();
%protgroup=();
%groupprot=();
@edges=();
# Création du fichier résultat
open( my $FhResultat, '>', $FichierResulat )
  or die("ERROR: can not opent the file $FichierResulat\n$!");
  if(exists $doc->{orthoXML}->{origin}){
    print "parsed methode : $doc->{orthoXML}->{origin}\n";
  }
  foreach $specie (keys(%{$doc->{orthoXML}->{species}})){
#     print "$specie\n";
    foreach $idg (keys(%{$doc->{orthoXML}{species}{$specie}{database}{genes}{gene}})){
#       print "idg $idg\n";
      if(exists($proteines{$idg})){
	die("ERROR: several proteins with the same id");
      }
      $proteines{$idg}=$doc->{orthoXML}->{species}->{$specie}->{database}->{genes}->{gene}->{$idg}->{protId};
#       print "$proteines{$idg}\n";
    }
  }
  foreach $group (keys(%{$doc->{orthoXML}{groups}{orthologGroup}})){
#     print "$group\n";
    foreach $p (keys(%{$doc->{orthoXML}{groups}{orthologGroup}{$group}{geneRef}})){
      $protgroup{$p}{$group}="";#storage of all groups
      $groupprot{$group}{$p}='';
#       print "$p \t $group\n";
    }
  }


#recherche de cliques si une proteine peut etre dans plusieurs groupes
$cliquesearch=0;
foreach $p(keys(%protgroup)){
  $nbg=keys(%{$protgroup{$p}});
  if($nbg>1){
    #recherche de clique
    $cliquesearch=1;
  }

}
#-----------------------
$cliquesearch=1;#aretirer

#-----------------------
@allcliques=();
if($cliquesearch){
  #1) creer toutes les paires de proteines
  foreach $g(keys(%groupprot)){
    @protg=(keys(%{$groupprot{$g}}));
    for($i=0;$i<$#protg;$i++){
      for($j=$i+1;$j<=$#protg;$j++){
	if($protg[$i]>$protg[$j]){
	  @tmp=($protg[$j],$protg[$i]);
	}
	else{
	  @tmp=($protg[$i],$protg[$j]);
	}
	push(@edges,[@tmp]);
      }
    }
  }
  #2) recherche de cliques
  $flag=1;
  $k=2;
  while($flag==1){
#     print "k $k\n";
    @cliques=getcliques($k,\@edges);
#     print join("\n",@cliques)."\n";
    $nbcliques=@cliques;
#     print "nbcliques $nbcliques\n";
    @allcliques=(@allcliques,@cliques);
    if($nbcliques==0){
      $flag=0;
    }
    $k++;
  }
  #selection des cliques a ecrire:
  %vu=();
  for($i=$#allcliques;$i>=0;$i--){
    @tmp=split(';',$allcliques[$i]);
    $garde=1;
    foreach $t(@tmp){
      if(exists($vu{$t})){
	$garde=0;
      }
    }
    if($garde==1){
#       print  $FhResultat "$allcliques[$i]\n";
      @prottmp=split(';',$allcliques[$i]);
      for($i=0;$i<=$#prottmp;$i++){
	$prottmp[$i]=$proteines{$prottmp[$i]};
      }
      print $FhResultat join(';',@prottmp),"\n";
      foreach $t(@tmp){
	$vu{$t}="";
      }
    }
  }
}
else{
  foreach $gr(keys(%groupprot)){
#     print $FhResultat join(";",keys(%{$groupprot{$gr}})),"\n";
    @prottmp=keys(%{$groupprot{$gr}});
    for($i=0;$i<=$#prottmp;$i++){
      $prottmp[$i]=$proteines{$prottmp[$i]};
    }
    print $FhResultat join(';',@prottmp),"\n";
  }
}


close($FhResultat);


sub getcliques {

     my ($k,$edges) = @_;
     my @cliques = ();
     my @vertices = ();
     
      @vertices = edges2vertices(@{$edges});

     my   $string = (join ',' => @vertices) .  ';'
                . (join ',' => map "$_->[0]-$_->[1]", @{$edges});

     my  $regex = '^ .*\b '
               . join(' , .*\b ' => ('(\d+)') x $k)
               . '\b .* ;'
               . "\n";

    for (my $i = 1; $i < $k; $i++) {
            for (my $j = $i+1; $j <= $k; $j++) {
                $regex .= '(?= .* \b ' . "\\$i-\\$j" . ' \b)' . "\n";
            }
        }

     # Backtrack to regain all the identified k-cliques (Credit Mike Mikero)
     $regex .= '(?{ push (@cliques, join(";", map $$_, 1..$k) ) })(?!)';
     $string =~ /$regex/x; 
     
     return sort @cliques;
}

#----Subroutines -------------------
sub edges2vertices {
  my @edges = @_;
  my %hTemp;
  my @vertices;
  
 my  @aTemp = map{@{$_}} @edges;
      @hTemp{@aTemp}  = ();
  @vertices = sort keys %hTemp;   
  return @vertices;  
}

sub edges2vertices_slow {
  #AoA to uniq array;

  my @edges = @_;
  my @vertices;
  my @uniqv; 
  
   for my $i ( 0 .. $#edges ) {
               for my $j ( 0 .. $#{$edges[$i]} ) {
                   push @vertices, $edges[$i][$j];
               }
           }

       @uniqv = sort keys %{{map {$_,1} @vertices}};
    return @uniqv;
}
