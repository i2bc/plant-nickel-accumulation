#!/usr/bin/env perl

# seqxml2any
# Copyright Dave Messina <dmessina@cpan.org>
# First version: 2010-01-15
#
# For questions, help, support, or suggestions:
# please email info@orthoxml.org
#
# For more information on SeqXML:
# please visit http://seqxml.org

use strict;
use warnings;
use Bio::SeqIO;
use Getopt::Long;
use File::Basename;
use File::Spec;

my ( $format );
GetOptions(
    "format=s"        => \$format,
);

my $usage = "
seqxml2any - convert SeqXML to (almost) any sequence format

Usage: seqxml2any --format <format> seqfile1.xml seqfile2.xml ... seqfileN.xml

Required Parameters:
--format        <format>   the format of the output sequence

Examples:
any2seqxml --format fasta mynewdata.xml
any2seqxml --format genbank first_run.xml second_run.xml

More Information:
The output filename will be the same as the input filename
except the suffix will be the format type. For example,
my_data.xml will become my_data.fasta if --format fasta.

Supported output formats are anything BioPerl supports.
Here are a few common ones:
    fasta
    genbank
    embl
    interpro
    pir
    swiss

For the complete list of supported formats,
    see http://www.bioperl.org/wiki/HOWTO:SeqIO
For more info on SeqXML, see http://seqxml.org
";
die "ERROR:\n$usage" unless ( $format and @ARGV );

foreach my $infile (@ARGV) {

    # input stream
    my $in = Bio::SeqIO->new(
        '-format' => 'seqxml',
        '-file'   => $infile
    );

    # create output filename
    # will be the same as the input filename except the suffix will
    # be the format name
    my ( $base, $path, $suffix ) = fileparse( $infile, qr/\.[^.]*/ );
    my $newbase = $base . '.' . $format;
    my $outfile = File::Spec->catfile( $path, $newbase );

    # output stream
    my $out = Bio::SeqIO->new(
        '-format' => $format,
        '-file'   => ">$outfile",
    );
    
    # main loop
    # each time, one seq record is read from the input
    # and written to the output
    while ( my $inseq_obj = $in->next_seq ) {
        $out->write_seq($inseq_obj);
    }
}
