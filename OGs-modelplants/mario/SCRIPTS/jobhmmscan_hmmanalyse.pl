#!/usr/bin/perl -w


$TMP=$ARGV[0];
$nbinter=$ARGV[1];
$s=$ARGV[2];
$DScript=$ARGV[3];
$arge=$ARGV[4];
$argpa=$ARGV[5];







system("hmmscan --cpu 1 -E 0.1 -o $TMP/hmmscan_out$nbinter/$s.out_o --tblout $TMP/hmmscan$nbinter/$s/$s.out $TMP/BD$nbinter/hmm_all $TMP/fasta_seq_without_groups$nbinter/$s > $TMP/OUT/hmmscan$nbinter$s.out >> $TMP/ERR/hmmscan$nbinter$s.err") and die("ERROR HMMSCAN comparison");
system("perl $DScript/analyse_cmp.pl $TMP/hmmscan$nbinter/$s/ $arge $TMP/hmmscan_out$nbinter/$s.out_o $argpa $TMP/BD$nbinter/hmm_all > $TMP/hmmanalysis$nbinter/$s.out ") and die("ERROR filter HMM comparison results");

        
