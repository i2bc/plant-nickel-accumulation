#!/usr/bin/perl -w
use Graph;
use Graph::Undirected;
use strict;

# Nouveau graphe
my $g = Graph::Undirected->new;

$g->add_edges("1","2");
$g->add_edges("2","1");
$g->add_edges("1","2");
$g->add_edges("2","3");
$g->add_edges("3","2");
$g->add_edges("8","9");
$g->add_edges("3","8");
$g->add_edges("1","3");
$g->add_edges("9","3");
$g->add_edges("3","9");

my @b = $g->bridges();

for my $br (@b){
    print $br->[0]."<->".$br->[1]."\n";
}

my @e = $g->edges();
for my $e (@e){
    print $e->[0]."-".$e->[1]."\n";
}
