#!/usr/bin/perl -w
#make the intersection
#all groups with a size greater than 2 are conserved.

#how use this script: 
#perl test_inter_v2.pl [file groups 1] [file groups 2] >[result file]

#[file groups 1] and [file groups 2] :
#each file contain all groups for one method.
#each group is on a line 
#protein of groups are sepered by ;

#the file is like that (don't put #):
#prot1;prot2;prot3
#prot5;prot6;prot7;prot8

#here prot1;prot2;prot3 is one group of 3 proteins
#and prot5;prot6;prot7;prot8 is one other group

#the result file will be in the same format but will contain only the intersection between [file groups 1] and [file groups 2]

#attention : in each file be care of each protein is not in more than 1 group. Is a protein could be in two group this script will not give you the good result

open IN1,$ARGV[0];
$nb_line=0;
%hash1=();
while (<IN1>)#travel of the first file pass through the groups found by the first method 
{
  $nb_line++;#the number of the line will be the ID of the group
  chomp;
  @tmp=split(/;/,$_);
  foreach $val (@tmp)
  {
    $hash1{$val}=$nb_line;#association between protein name and potential group number (line number in A)
  }
}
close IN1;

open IN2,$ARGV[1];
$nb_line=0;
while (<IN2>) #pass through the groups found with the second method
{
  $nb_line++;
  chomp;
  @tmp=split(/;/,$_);

  %hash_res=();
  foreach $val (@tmp) #for each protein of the current group (method 2)   
  {
	
    push @{$hash_res{$hash1{$val}}},$val if (exists $hash1{$val});#association between number of the protein in the first file and the protein name
    #it's a table, so all protein with are in the same time in one group in the second method (the line $_) and in the first method (ĥash1{$val} = ID group in first method) will be in the same table
  }
  foreach $val (keys %hash_res) #one key by table so one key by intersection
  {
    $nb2=@{$hash_res{$val}};
    if ($nb2 >= 2) #we conserve the intersection only is they are 2 or more proteins
    {
      $res_group=join(';',@{$hash_res{$val}});
    	print "$res_group\n";
    }
  }
}
close IN2;
