#!/usr/bin/perl -w
if ($#ARGV <0) 
{ 
	print "perl orthofinder2mario.pl <repertoire de travail>\n";
	exit();
}
my $work_dir = $ARGV[0] ;

system("mkdir $work_dir/mario_input") unless ( -e "$work_dir/mario_input") ;
$orthofinderfile = $work_dir.'/OrthoFinder/input/OrthologousGroups.txt';
print $orthofinderfile."\n";
$mariofile = $work_dir.'/mario_input/orthofinder_groups';

open (IN, $orthofinderfile);
open (OUT, ">$mariofile");
while ($l= <IN>) 
{ 
	chomp $l;
	@l = split (/:\s+/, $l);
	$group = $l[1];
	@prot = split (' ', $group);
	if(scalar(@prot) >1)
	{
		$group =~ s/ /;/g ;
		print OUT $group."\n";
	}
}
close IN;
close OUT;

