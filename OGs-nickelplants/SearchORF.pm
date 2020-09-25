package SearchORF;

sub search_one 
{#1	
	($sequence, $a) = @_;
	#print $sequence."\n";
	#Definition des codon start et stop
	$codon_start="ATG";
    $codon_stop="TAG|TGA|TAA";
    # recherche des regions sans codon STOP dans la sequence
	my $regexpr= "((?(?!$codon_stop)[ATCG]{3}){$a,})";
	
	#$x = "foobar";
    #$x =~ /foo(?!bar)/;  # pas de correspondance, 'bar' suit 'foo'
    #$x =~ /foo(?!baz)/;  # reconnue, 'baz' ne suit pas 'foo'
    #$x =~ /(?<!\s)foo/;  # reconnue, il n'y a pas de \s avant 'foo'
    
	$l= length($sequence);
	my $maxlength = 0;				
	pos $sequence = 0;	
	while ($sequence =~ m/$regexpr/g)				
	{ 			
		if (length($1) > $maxlength)
		{
			$maxlength = length($1);			
			$ORF = $1;
			$start = $-[1] +1;
			$end = $+[1];		
		}	
		pos($sequence) = length($`) +1;		
	}
	#print "start: ".$start."\n";
	#print "end: ".$end."\n";
	if ($sequence =~ m/$ORF($codon_stop)/)
	{
		$ORF .= $1;
		$end += 3;
	}
	if ($ORF =~ /((?(?!$codon_start)[ATCG]{3})*)($codon_start)/)
	{
		#print "first ATG pos: ",pos($ORF)."\n";
		$meth = $-[2];
		#print length($ORF)," ORF: ".$ORF."\n";
		#print "first methionine position: $meth\n";
	}
	return ($ORF, $start, $end, $meth); 			
}								
1;

