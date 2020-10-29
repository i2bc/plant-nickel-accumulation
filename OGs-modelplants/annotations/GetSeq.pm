package GetSeq;
use LWP::Simple;
use LWP::UserAgent;

# Extraction de la sequence $id de la base Uniprot
sub Extract_UP_seq
{
	my ($id,$ec,$fileOut,$refSpecie) = @_;
    if (($page_html = get('http://www.uniprot.org/uniprot/'.$id.'.txt')) || ($page_html = get('http://www.uniprot.org/uniprot/?query=replaces:'.$id.'&format=txt'))) 
    {
		#On recupere chaque ligne dans un tableau
		@Line = split(/\n/,$page_html);
		for ($i = 0; $i <= $#Line; $i ++)
		{
			# Recuperation de l'organisme
			if($Line[$i] =~ /^OS *([^ ].*)$/)
			{
				push (@Organism, $1) ;
			}
			# Recuperation du taxon ID
			elsif($Line[$i] =~ /^OX *NCBI_TaxID=([0-9]*);/)
			{
				$taxonId = $1 ;
			}
			# Recuperation du taxon
			elsif($Line[$i] =~ /^OC *([^ ].*)$/)
			{
				push (@Taxon, $1);
			}
			elsif ($Line[$i] =~ /^SQ /)  
			{
				$i ++;
				while ($Line[$i] =~ /^ /) 
				{
					$seq .= $Line[$i];
					$i ++;;                
				}
				$seq =~ s/ //g; # Suppression des espaces dans la sequence
			}
		}
		print "Warning no organism for $id in Uniprot\n" if ($#Organism < 0);
		print "Warning no taxonomy for $id in Uniprot\n" if ($#Taxon < 0);
		print "Warning no taxonomy ID for $id in Uniprot\n" if ($taxonId == 0);
		print "Warning no sequence for $id in Uniprot\n" if ($seq eq '');
		$$refSpecie{join(' ', @Organism)} = join(' ', @Taxon)."\t$taxonId" unless (exists($$refSpecie{join(' ', @Organism)}));
		open OUT,'>>'.$fileOut;
		chomp($id);
		print OUT "$id\t$ec\t".join(' ', @Organism)."\t".$$refSpecie{join(' ', @Organism)}."\t$seq\tMetaCyc\n" ;
		close(OUT);
	}
	else
	{
		print "Warning : Page not found for $id :\nhttp://www.uniprot.org/uniprot/$id\.txt\n'http://www.uniprot.org/uniprot/?query=replaces:$id&format=txt'\n";
	}
}


# Extraction de la sequence $id de la base de donnnees RefSeq
sub Extract_pid_refseq_seq
{
	#print STDERR "dans extract pid refseq seq \n";
    my ($id,$ec,$fileOut,$refSpecie) = @_;
    
	
    my $ua2 = new LWP::UserAgent;
    my $url2 ="http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&id=$id&retmode=text&rettype=gp&email=sandrine.grossetete\@igmors.u-psud.fr";
		#print STDERR $url2."\n";die;
	# Recuperation de la page associee a la sequence
    my $request2 = new HTTP::Request GET=> $url2;
    my $page2 = $ua2->request($request2);
    while (!($page2->is_success) && ($i < 5)) 
    {
        sleep(1);
        $request2 = new HTTP::Request GET=> $url2;
        $page2 = $ua2->request($request2);
        $i ++;
    }
	# Si il y a eu une erreur lors de la recuperation de la sequence id
    if (!($page2->is_success))    
    {
            print "Erreur de recuperation dans la base RefSeq  : $url2 !\n" . $page2->error_as_HTML . "\n\n";
    }
	# Si la page a ete recuperee avec succes
    else    
    {
		# Recuperation du contenu de la page dans $i
        $i = $page2->content ;

		# Recuperation des informations qui nous interessent dans la page
        @Line = split(/LOCUS/,$i);
		if ($#Line >= 0)
		{
			for ($i = 1; $i <= $#Line; $i ++)
			{
				$tmpSpecie = '-';
				$tmpSeq = '';
				if ($Line[$i] =~ /\nSOURCE *([^ ][^\n]*)\n/)
				{
					$tmpSpecie = $1 ; # Recuperation de l'organisme					
				}
				else
				{
					print "No specie found for $id in RefSeq\n";
				}
				if ($Line[$i] =~ /\nORIGIN([^\/]*)\//)
				{
					$tmpSeq = $1 ; # Recuperation de la sequence
					$tmpSeq =~ s/[\n 0-9]//g; # Suppresion de la position et des retours a la ligne
					if (($tmpSeq =~ /^[ATGCNnatgc]*$/) && ($Line[$i] =~ /\/translation="([ \na-zA-Z]*)"/))
					{
						$tmpSeq = $1;
						$tmpSeq =~ s/[\n ]//g; # Suppresion des espaces et des retours a la ligne
					}
				}
				#Enregistrement de la sequence
				if (($tmpSeq ne '') && !($tmpSeq =~ /^[ATGCNnatgc]*$/)) 
				{
					if ($seq eq '')
					{
						$seq = $tmpSeq;
						$specie = $tmpSpecie;
						$$refSpecie{$specie} = searchTaxon($specie) unless (exists($$refSpecie{$specie}));  # Recuperation du taxon
					}
					else
					{
						$seq = "Several result";
						$specie .= "\nSeveral result";
					}
				}
			}
			if ($seq eq '') 
			{
				print "Warning : No amino acis sequences for $id in RefSeq\n".join("\n",@Line)."\n\n" ;
			}
			elsif($seq eq 'Several result')
			{
				print "Warning : Several results for $id in RefSeq (none result saved)\n".join("\n",@Line)."\n\n" ;
			}
			else
			{
				if (exists($$refSpecie{$specie}))
				{
					open OUT,'>>'.$fileOut;
					chomp($id);
					print OUT "$id\t$ec\t$specie\t$$refSpecie{$specie}\t$seq\tMetaCyc\n" ;
					close(OUT);
					if ($specie eq '-')
					{
						print "Warning : No specie for $id in Ref Seq\n".join("\n",@Line)."\n\n" ;
					}
				}
				else
				{
					print "Warning : none taxonomy found for $specie : \n$id\t$ec\t$specie\t$seq\tMetaCyc\n" ;
				}
			}
		}
		else
		{
			print "Warning : No results for $id in RefSeq\n" ;
		}
     }
}

# Extraction de la sequence $id de la base de donnnees PDB
sub Extract_PDB_seq
{
	#print STDERR "dans extract pdb seq\n";
    my ($id,$ec,$fileOut,$refSpecie) = @_;
    my @Line = ();
    my %Seq = ();
    my $header = '';
    my $specie = '-';
    my $seq = '';
		my $i = '';
	# Recuperation des informations reliees a $id
    $i = get('http://www.rcsb.org/pdb/explore/explore.do?structureId='.$id) ;
    $specie = $1 if ($i =~ /TreeEntity[^>]*>([^>]*)<\/a>/);
    $specie =~ s/^[^A-Za-z]*//g;
    $specie =~ s/[^A-Za-z]*$//g;
    $$refSpecie{$specie} = searchTaxon($specie) unless (exists($$refSpecie{$specie}));
	# Recuperation de la sequence $id
    $seq = get('http://www.rcsb.org/pdb/files/fasta.txt?structureIdList='.$id) ;
    #$seq du type ">2F9Y:A|PDBID|CHAIN|SEQUENCE
		#MGSSHHHHHHSSGLVPRGSHMSLNFLDFEQPIAELEAKIDSLTAGSRQDEKLDINIDEEVHRLREKSVELTRKIFADLGA
		#WQIAQLARHPQRPYTLDYVRLAFDEFDELAGDRAYADDKAIVGGIARLDGRPVMIIGHQKGRETKEKIRRNFGMPAPEGY
		#RKALRLMQMAERFKMPIITFIDTPGAYPGVGAEERGQSEAIARNLREMSRLGVPVVCTVIGEGGSGGALAIGVGDKVNML
		#QYSTYSVISPEGCASILWKSADKAPLAAEAMGIIRPRLKELKLIDSIIPEPLGGAHRNPEAMAASLKAQLLADLADLDVL
		#STEDLKNRRYQRLMSYGYA
		#>2F9Y:B|PDBID|CHAIN|SEQUENCE
		#MSWIERIKSNITPTRKASIPEGVWTKCDSCGQVLYRAELERNLEVCPKCDHHMRMTARNRLHSLLDEGSLVELGSELEPK
		#DVLKFRDSKKYKDRLASAQKETGEKDALVVMKGTLYGMPVVAAAFEFAFMGGSMGSVVGARFVRAVEQALEDNCPLICFS
		#ASGGARMQEALMSLMQMAKTSAALAKMQERGLPYISVLTDPTMGGVSASFAMLGDLNIAEPKALIGFAGPRVIEQTVREK
		#LPPGFQRSEFLIEKGAIDMIVRRPEMRLKLASILAKLMNLPAPNPEAPREGVVVPPVPDQEPEA"
    if ($seq ne '')
    {#si sequence trouvée
        @Line = split(/\n/,$seq);
        $seq = '';
				# Pour chaque ligne retournee (peut avoir plusieurs sequences)
				print "scalar ligne=".scalar(@Line)."\n";
        for ($i = 0; $i < scalar(@Line); $i ++)
        {
						# Si nouvelle sequence trouvee
            if ($Line[$i] =~ /^>/)
            {
								# On enregistre la precedente
                $Seq{$seq}{$id} = '' if ($seq ne '');
								# Recuperation du nouvel ID
                $id  = $1 if ($Line[$i] =~ />([^\|]*)\|/);
                $seq = '';
            }
						# Recuperation de la suite de la sequence
            else
            {
                $seq .= $Line[$i];
            }
        }
				# Enregistrement de la derniere sequence
        $Seq{$seq}{$id} = '' if ($seq ne '');
				# Enregistrement des sequences retournees (peut avoir plusieurs IDs avec meme sequence
        open OUT,'>>'.$fileOut;
				# Pour chaque sequence proteique
        foreach $seq (keys(%Seq))
        {
						# Tableau des IDs partageant la meme sequence
            @Line = sort(keys(%{$Seq{$seq}})); 
						# On enregistre uniquement le premier ID
						my $prem=shift(@Line);
						chomp($prem);
            print OUT $prem."\t$ec\t$specie\t$$refSpecie{$specie}\t$seq\tMetaCyc\n";
        }
        close(OUT);
    }
    else 
    {
        print "Pb with $id (PDB): No res\n";
    }  
}

# Recherche des informations taxonomiques reliees a une espece $specie donnee
sub searchTaxon  {
#print STDERR "searchtaxon\n";
    my ($specie) = @_;
    my $res = '';
		my $taxon = '';
		my $idTaxon = '';
    my @Line = ();
		my $tmpSpecie = '';
		my $otherInfo = '';
		if ($specie =~ /^([A-Z][^ ]* [^ ]*) (.*)$/){
			$tmpSpecie = $1;
			$otherInfo = $2;
		}
		elsif ($specie =~ /^[a-z ]* ([A-Z][^ ]* [^ ]*) (.*)$/){
			$tmpSpecie = $1;
			$otherInfo = $2;
			$tmpSpecie =~s/ /+/g;
			$otherInfo =~s/ /+/g;
			$otherInfo =~s/[\(\)]//g;
		}
		else{
			$tmpSpecie = $specie;
			$tmpSpecie =~s/ /+/g;
		}
		$tmpSpecie=~s/\'//g;#ajout cecile
		if ($otherInfo ne ''){
			($res = get('http://www.uniprot.org/taxonomy/?query=%22'.$otherInfo.'%22+AND+(common%3A%22'.$tmpSpecie.'%22+OR+scientific%3A%22'.$tmpSpecie.'%22)&sort=score'))  || die ("Problem with specie $specie\n".'http://www.uniprot.org/taxonomy/?query=%22'.$otherInfo.'%22+AND+(common%3A%22'.$tmpSpecie.'%22+OR+scientific%3A%22'.$tmpSpecie.'%22)&sort=score');
			unless ($res =~ /<p class="summary">/){
				($res = get('http://www.uniprot.org/taxonomy/?query=.'.$tmpSpecie.'+'.$otherInfo.'&sort=score'))  || die ("Problem with specie $specie\n".'http://www.uniprot.org/taxonomy/?query=.'.$tmpSpecie.'+'.$otherInfo.'&sort=score');
			}
		}
		else{
			($res = get('http://www.uniprot.org/taxonomy/?query=(common%3A%22'.$tmpSpecie.'%22+OR+scientific%3A%22'.$tmpSpecie.'%22)&sort=score'))  || die ("Problem with specie $specie\n".'http://www.uniprot.org/taxonomy/?query=(common%3A%22'.$tmpSpecie.'%22+OR+scientific%3A%22'.$tmpSpecie.'%22)&sort=score');
		}
		unless ($res =~ /<p class="summary">/){
			($res = get('http://www.uniprot.org/taxonomy/?query=.'.$tmpSpecie.'&sort=score'))  || die ("Problem with specie $specie\n");
		}
							
    if ($res ne '')    {
		@Line = split(/<tr>/,$res);
		$res = '';
		for (my $i = 0; $i <= $#Line; $i ++){
			# if ($Line[$i] =~ /<p class="summary">([^\<]*)<\/p>.*uniprot\/\?query=organism:([0-9]*)\+/){
				# $taxon = $1;
				# $idTaxon = $2;
			# }
			if ($Line[$i] =~ /<a href=".\/([0-9]*)">([^<]*)<\/a>.*<p class="summary">([^\<]*)<\/p>/){
				$idTaxon = $1;
				$tmpSpecie = $2;
				$taxon = $3;
				$res = "$taxon\t$idTaxon" ;

				unless ($res =~ /related/){
					$res =~ s/&rsaquo;/;/g;
					$res =~ s/ ; /_/g;
					return $res;
				}
			}
			elsif (($res eq '') && ($Line[$i] =~ /<p class="summary">([^\<]*)<\/p>/)){
				$res = "$1\t-" ;
				$res =~ s/&rsaquo;/;/g;
				$res =~ s/ ; /_/g;
			}
		}
		if ($res eq ''){
			print "No taxonomy for $specie\n".'http://www.uniprot.org/taxonomy/?query=%22'.$otherInfo.'%22+AND+(common%3A%22'.$tmpSpecie.'%22+OR+scientific%3A%22'.$tmpSpecie.'%22)&sort=score'."\n";
			print join("\n",@Line)."\n";
		}
		else{
			print "No taxonomy ID for $specie ($res)\n";
			print join("\n",@Line)."\n";
		}
		return "-\t-";
    }
}

__END__

