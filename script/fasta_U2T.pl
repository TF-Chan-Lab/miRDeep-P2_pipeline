#!/bin/perl
#Author:ALin
#Purpose: Input a transcript fasta and output a DNA fasta.
#Usage: perl fasta_U2T.pl <input> <output>

if(@ARGV != 2){
	print "Usage: perl fasta2fasta.pl <input> <output>\n";
	exit;
}

use strict;
use Bio::Seq;
use Bio::SeqIO;


my $in = shift;
my $out = shift;

my $seq_in = Bio::SeqIO->new(-file => $in, -format => 'fasta',);
my $seq_out = Bio::SeqIO->new(-file => ">$out", -format => 'fasta');

while(my $seq_obj = $seq_in->next_seq){
	my $new_seq = $seq_obj->seq;
	$new_seq =~ s/U/T/g;
	$seq_obj->seq($new_seq);
	$seq_out->write_seq($seq_obj);
}





