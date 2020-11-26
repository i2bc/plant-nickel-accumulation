#!/usr/bin/perl -w
# Coded by C. Drevet
use File::Basename;
if ($#ARGV <6) 
{ 
	print "back_to_name-list.pl <inputfile> <valeur de split inputfile> <fichier references numero/genome> <valeur de split reference genomes> <dossier de reference numero/proteines> <valeur de split numero/proteines> <tag>\n";
	exit();
}
my $f = $ARGV[0]; # the file to process 
my $split_f = $ARGV[1];
my $file_numero_genomes = $ARGV[2];
my $split_g = $ARGV[3];
my $dir_numero_prot = $ARGV[4].'/';
my $split_p = $ARGV[5];
my $tag = $ARGV[6];
# build a table of genome codes
open (REF, $file_numero_genomes);
my @genome;
while($lg = (<REF>))
{
	chomp $lg;
	@lg = split(/$split_g/,$lg);
	$lg[1] =~ s/\.fa//;
	$genome[$lg[0]] = $lg[1];
}
close REF;
#open the file to process 

open(OUT,">named_".basename($f).$tag);        
open (IN, $f);
while (<IN>) { $fcontent .= $_; }
print $fcontent."\n\n\n";
@f = split (/\n/, $fcontent);
print $f[0]."\n";
print $f[1]."\n";
foreach $cluster (@f)
{
		chomp $cluster;	
		print  $cluster."\n";
		@group = split(/$split_f/, $cluster);
		foreach $p (@group)
		{
			print $p."\n";
			@p = split (/_/, $p);
			$genome_nb = $p[0]; 
			print $genome_nb."\t".$p[1]."\n";
			$genome_code = $genome[$genome_nb];
			print "CODE".$genome_code."\n";
			print $dir_numero_prot.$genome_nb.".name\n"; 
			open (GEN,$dir_numero_prot.$genome_nb.".name");
			while ($n = (<GEN>))
			{
				chomp $n;
				@n = split(/$split_p/,$n);
				if($n[0] eq $p) 
				{ 
					$hp = $n[1];
					$fcontent =~ s/^$p;/$hp;/; # premiere proteine du fichier (plusieurs proteines sur la premiere ligne)
					$fcontent =~ s/^$p\n/$hp\n/; # premiere proteine du fichier (une seule proteine sur la premiere ligne du fichier, groupes de paralogs)
					$fcontent =~ s/\n$p;/\n$hp;/; # proteines des debuts de ligne (2-fin)
					$fcontent =~ s/\n$p\n/\n$hp\n/; # proteines en debut et fin de ligne (groups de parologs)
					$fcontent =~ s/;$p/;$hp/; # autres cas					
					last;
				}
			}
			close GEN;	
		}		
}		
print OUT $fcontent;
close OUT;
