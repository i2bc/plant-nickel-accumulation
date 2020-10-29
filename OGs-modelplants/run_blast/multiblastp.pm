package multiblastp;
use File::Basename;

sub run_formatdb
{	
	my $fastafile = shift(@_);
	unless (-e $fastafile.'.phr')
	{
		system("makeblastdb -dbtype prot -in $fastafile");
		
	}
	return  $fastafile." formatted \n";
}

sub run_blastp
{
	my ($query, $hit, $out_dir) = @_;
	print "running blast beetwen ",$query," and ",$hit,"\n";
	
	@suffix = (".fa", ".fasta");
	$queryname = basename($query,@suffix); 
	$hitname = basename($hit, @suffix);
	unless ( -e $out_dir.'/'.$queryname."_vs_".$hitname.".blast") 
	{
		# blasts format de sortie:  qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore query length subject length
		system("blastp -query $query -db $hit  -outfmt \"7 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen\" -seg yes -out  ".$out_dir.'/'.$queryname."_vs_".$hitname.".blast");
	}
	return $queryname."_versus_".$hitname." done\n";
}
1;
