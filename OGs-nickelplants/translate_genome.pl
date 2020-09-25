#!/usr/bin/perl -w

# Recherche des ORFs dans les contigs obtenus par RNA-seq (marioplants project)
use Bio::SeqIO;
use File::Basename;
#use lib "/home/chris/tmarioplants";

use SearchORF;
	$min_aa = $ARGV[0];
	$file = $ARGV[1];
	@suffix = (".fa");
	
	($name,$path,$suffix) = fileparse($file,@suffix);
	#print "*".$name,"|","*",$path,"|","*",$suffix, "\n";

	# read contig
	$in  = Bio::SeqIO->new(-file => $file ,
                           -format => 'fasta');
    # store longer ORF (3 forward phases) 
	$outnt = Bio::SeqIO->new(-file => ">COGs/orf/$name.fa" ,
                           -format => 'fasta');
    # store longer ORF translation
    $outp = Bio::SeqIO->new(-file => ">COGs/proteome/$name.fa" ,
                          -format => 'fasta'); 
    # create a file with oriented RNA
    $RNA = Bio::SeqIO->new(-file => ">COGs/oriented_contigs/$name.fa" ,
                          -format => 'fasta');      
    # store infos (tabular format)
    open OUT, ">COGs/orf_tab/$name.tab";
    $count = 0; 
    $countrev =0;                
	while (my $seq = $in->next_seq()) 
	{
		$count++;
		my $id = $seq->display_id(); 
		#print $id,"\n"; 
		my $desc = $seq->desc();        
		my $sequence = $seq->seq(); 
		#print $sequence."\n";
		my $revcomp = $seq->revcom();
		my $revcompseq = $revcomp->seq;
		#print $revcompseq."\n";
		print OUT $id."\t".length($sequence)."\t";	
		($s, $start, $end, $meth) =  SearchORF::search_one($sequence, $min_aa);
		($srev, $startrev, $endrev, $methrev) =  SearchORF::search_one($revcompseq, $min_aa);
			
		if(length($srev) > length($s) )
		{ 
			$countrev++;
			$s = $srev; 
			$start = $startrev; 
			$end = $endrev; 
			$meth = $methrev;
			$revcomp->desc($revcomp->desc."[revcomp]");
			$RNA->write_seq($revcomp);
			$strand = "minus";
		}
		else 
		{ 
			$RNA->write_seq($seq);
			$strand = "plus";
		}
		print OUT $strand."\t".length($s) ."\t".sprintf ("%0.0f", length($s)/length($sequence)*100)."\t".$start."\t".$end."\t".$meth."\n";
		$os= Bio::Seq->new (-display_id => $id, -seq => $s, -format=> 'fasta');
		$outnt->write_seq ($os);
		$op = $os->translate();
		$outp->write_seq($op);								
	}
	close OUT;	
	print "$file: $count sequences ; $countrev reversed sequences; (".sprintf("%0.0f",$countrev/$count*100).")\n";
	


