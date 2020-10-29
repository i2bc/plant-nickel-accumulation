package metacyc2fasta;
use strict;
use warnings;
use XML::LibXML;
sub uniprot2fasta
{
	# extraction des données a partir d'un ID Uniprot de MetaCyc
	# La  methode basique Load_XML produit des erreurs (deprecated IDs)
	# Le script exit quand le fichier XML est vide alors que l'URL demandée reste accessible
	# A revoir pour gerer l'exception (ajout d'un test avec eval actuellement)


	my ($id) = @_; 
	#get data from UNIPROT using Protein IDs (XML formatted)
	use WWW::Curl::Easy;
	my $curl = WWW::Curl::Easy->new;        
	$curl->setopt(CURLOPT_HEADER,0); # boolean: 0 no header/1 header
	$curl->setopt(CURLOPT_URL, "http://www.uniprot.org/uniprot/$id.xml"); 

	# A filehandle, reference to a scalar or reference to a typeglob can be used here.
	my $response;
	$curl->setopt(CURLOPT_WRITEDATA,\$response);
	# Starts the actual request
	my $retcode = $curl->perform;
	# Looking at the results...
	if ($retcode != 0) 
	{
		# Error code, type of error, error message
		print("An error happened: $retcode ".$curl->strerror($retcode)." ".$curl->errbuf."\n");
		return -1;
	}
	if (!$response) { return -1;}	
	else
	{
		my @data = ();
		my $dataset = '';
		my $org = '';
		my $ec_list = '';
		eval
		{
			my $parser = XML::LibXML->new();
			my $struct = $parser->parse_string($response);
			my $root = $struct -> getDocumentElement();
			
			my ($entry) = $root->getElementsByTagName('entry');			
			if($entry->hasAttribute('dataset'))
			{  
				my @atts = $entry->getAttributes();
				foreach my $at (@atts) 
				{
					my $db = $at -> getName();
					if ($db eq 'dataset')
					{
						$dataset = $at -> getValue(); # recherche le l'origine de la sequence (Swiss-Prot ou TrEMBL)
						push (@data, $dataset);
					}
				}
			}
			my ($name) = $root->getElementsByTagName('name');
			my $n = $name->textContent; 
			push(@data,$n);	# extraction de l'ID de la sequence (ID specifique utilise dans le fichier uniprot_sprot.xml) 
			
			my ($organism) = $root->getElementsByTagName('organism');
			my @org_names = $organism->getElementsByTagName('name');
			foreach my $name (@org_names)
			{  
				my @name_at = $name->getAttributes();
				foreach my $at (@name_at) 
				{
					if ($at->getName() eq 'type' and $at->getValue eq "scientific")
					{
						$org = $name->textContent;
						push (@data, $org);
					}
				}
			}		
			my ($recname) = $root->getElementsByTagName('recommendedName'); # unique tag ou no tag
			if ($recname)
			{
				my ($fullName) = $recname->getElementsByTagName('fullName'); # unique tag		
				my $fn = $fullName-> textContent; 
				push (@data, $fn); # extraction de la fonction de la sequence 
				my @ecnodes = $recname->getElementsByTagName('ecNumber'); # multiple tag  or no tag
				foreach my $n (@ecnodes)
				{  
					$ec_list .= $n->textContent.' '; # recherche des EC numbers a l'interieur du tag recommendedName
				}
			}
			else # recherche du tag fullname si pas de tag recommendedName
			{
				my ($fullName) = $root->getElementsByTagName('fullName'); # unique tag
				my $fn = $fullName-> textContent;
				push (@data, $fn);
			}
			push(@data, $ec_list);	# extraction des EC numbers 		
			my @sequence = $root->getElementsByTagName('sequence'); # multiple tag			
			foreach my $s (@sequence)
			{  
				if($s->hasAttribute('length')) # attribut caracterisant le tag contenant la séquence protéique
				{
					my $seq = $s->textContent; #obtenir la séquence formatée
					$seq =~ s/\n//g;
					push (@data, $seq) or die ("no sequence found!");# extraction de la séquence
				} 
			}
		}; warn $@ if $@;
		return @data;	
	}
}
sub ncbi2fasta
{
		my ($id) = @_; 
		my @data = ();
		#get data from NCBI using Protein IDs (XML formatted)
		use WWW::Curl::Easy;
		my $curl = WWW::Curl::Easy->new;
        
		$curl->setopt(CURLOPT_HEADER,0); # boolean: 0 no header/1 header
		$curl->setopt(CURLOPT_URL, "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&id=$id&retmode=XML"); 

		# A filehandle, reference to a scalar or reference to a typeglob can be used here.
		my $response;
		$curl->setopt(CURLOPT_WRITEDATA,\$response);
		# Starts the actual request
		my $retcode = $curl->perform;
		# Looking at the results...
		if ($retcode != 0) 
		{
			#Error code, type of error, error message
			print("An error happened: $retcode ".$curl->strerror($retcode)." ".$curl->errbuf."\n");
			next;
		}
		if (!$response) {}	
		else
		{
			
			
			my $ncbi_id = '';
			my $org = '';
			my $function = '';
			my $product = '';
			my $ec_list = '';
			my $seq = '';
			eval
			{
				my $parser = XML::LibXML->new();
				my $struct = $parser->parse_string($response);
				my $root = $struct -> getDocumentElement();
				my ($moltype) = $root->getElementsByTagName('GBSeq_moltype');
				my $type= $moltype->textContent;
				if ($type eq "AA")
				{
					my ($locus) = $root->getElementsByTagName('GBSeq_locus');
					$ncbi_id = $locus->textContent; 
					my ($def) = $root->getElementsByTagName('GBSeq_definition');
					my $definition = $def->textContent;
					my @definition = split(/ \[|\]/ , $definition);
					$function = $definition[0];
					$org = $definition[1];
					my ($sequence) = $root->getElementsByTagName('GBSeq_sequence'); # unique tag		
					$seq = $sequence->textContent; #obtenir le texte entre les balises
					my @qualifiers = $root->getElementsByTagName('GBQualifier'); 
					foreach my $q (@qualifiers)
					{  
						my ($n)= $q->getElementsByTagName('GBQualifier_name');
						my $qname = $n->textContent;
						if ($qname eq "organism")
						{
							my ($qval) = $q->getElementsByTagName('GBQualifier_value');	
							if (!$org) { $org = $qval->textContent;}
							 			
						} 
						if ($qname eq 'EC_number')
						{
							my ($qval) = $q->getElementsByTagName('GBQualifier_value');	
							$ec_list .= $qval->textContent.' '; 			
						}
						if ($qname eq "product")
						{
							my ($qval) = $q->getElementsByTagName('GBQualifier_value');	
							$product = $qval->textContent; 			
						}
						if ($qname eq "function")
						{
							my ($qval) = $q->getElementsByTagName('GBQualifier_value');	
							if (!$function) {$function = $qval->textContent;} 
						}
					}	
				}	
				@data = ($ncbi_id, $org, $function, $product, $ec_list, $seq);		
			}; warn $@ if $@;
			return @data;	
		}
	}

1;
