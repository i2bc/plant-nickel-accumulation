#!/usr/bin/perl
use warnings;
#save the groups in orthoXML format
#Cécile Pereira
#07 march 2014

$FG=$ARGV[0];#folder with groups fasta files
$OXF=$ARGV[1];#orthoxml final file
#----------------------------------
#1) header of the orthoXML file
#----------------------------------
orthoxmlheader($OXF);

#----------------------------------
#2)foreach specie print sequences at least in one group and ortholog groups
#----------------------------------
chdir($FG);
@files=glob("*");
%genomeid=();#specie idprotein geneID
%groupido=();#orthologGroupID geneID
$nb=1;#gene id for the orthoXML file
foreach $f(@files){
  $nbf=$f;
  $nbf=~s/\.fa$//;
  open(IN,$f)or die("ERROR: can't open the file $f");
  while(<IN>){
    chomp;
    if(/>(.+)/){
      @tmp=split(/\t/,$1);
#      print join(';',@tmp);
      $genomeid{$tmp[1]}{$tmp[0]}=$nb;
      $groupido{$nbf}{$nb}="";
      $nb++;
    }
  }
  close(IN);
}
orthoxmlspecies($OXF);
orthoxmlgroups($OXF);

#----------------
#functions
#----------------
#create the header of the orthoxml file
sub orthoxmlheader{
  ($fio)=@_;
  open(OUT,">$fio") or die("ERROR: can't open the output file for orthoXML $fio");
  $text=<<"TEX";
<orthoXML xmlns="http://orthoXML.org/2011/" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" version="0.3" origin="meta-approach, Pereira C. et al" originVersion="1.0" xsi:schemaLocation="http://orthoXML.org/2011/ http://www.orthoxml.org/0.3/orthoxml.xsd">
TEX
  print OUT $text;
  close(OUT);
}

#create species part
sub orthoxmlspecies{
  ($fio)=@_;
  open(OUT,">>$fio");
  foreach $g(keys(%genomeid)){
    print OUT "\t<species name=".'"'.$g.'"'.">\n\t\t<database>\n\t\t\t<genes>\n";
    foreach $p(keys(%{$genomeid{$g}})){
      $text=<<"TEX";
\t\t\t\t<gene id="$genomeid{$g}{$p}" protId="$p"/>
TEX
      print OUT $text;
    }
    $text=<<"TEX";
\t\t\t</genes>
\t\t</database>
\t</species>
TEX
    print OUT $text;
  }
  close(OUT);
}


#print groups
sub orthoxmlgroups{
  ($fio)=@_;
  open(OUT,">>$fio");
  print OUT "\t<groups>\n";
  foreach $g(keys(%groupido)){
    print OUT "\t\t<orthologGroup id=\"$g\">\n";
    foreach $p(keys(%{$groupido{$g}})){
      $text=<<"TEX";
\t\t\t<geneRef id="$p">
\t\t\t</geneRef>
TEX
      print OUT $text;
    }
    print OUT "\t\t</orthologGroup>\n";
  }
  print OUT "\t".'</groups>'."\n";
  print OUT '</orthoXML>'."\n";
  close(OUT);
}
