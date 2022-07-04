#!/bin/perl
#Author: ALin
#Purpose: To combine fasta entries with the same sequence.
#Usage: perl unique_fasta_${version}.pl <in file> <output> <prefix>
#Change log:
#	v1.2	2020-07-06	Change the delineator of the description from "," to "&".
#	v1.3	2022-07-04	Force all letters in a sequence to be capital.

use strict;
use Bio::SeqIO;
use Bio::Seq;

my $version = "v1.3";

if(@ARGV != 3){
	print "Usage: unique_fasta_${version}.pl <in file> <output> <prefix>\n";
	exit;
}

my $in = shift @ARGV;
my $out = shift @ARGV;
my $prefix = shift @ARGV;
my %seq = ();

#open(IN, "<", $in) or die"Cannot open $in!\n";
#open(OUT, ">", $out) or die"Cannot create $out!\n";

my $seq_in = Bio::SeqIO->new( -format => 'fasta', -file => $in,);
my $seq_out =  Bio::SeqIO->new( -format => 'fasta', -file => ">$out",);
while(my $seq_obj = $seq_in->next_seq){
	my $temp_seq = uc($seq_obj->seq);
	my $temp_id = $seq_obj->id;
	unless(exists $seq{$temp_seq}){
		%{$seq{$temp_seq}} = ();
	}
	unless(exists $seq{$temp_seq}{$temp_id}){
		$seq{$temp_seq}{$temp_id} = 1;
	}
}

my $count = 1;

foreach my $seq (keys %seq){
	my $temp_id = $prefix . "_$count";
	my $temp_desc = "";
	foreach my $id (keys %{$seq{$seq}}){
		$temp_desc .= "$id&";
	}
	my $seq_obj = Bio::Seq->new();
	$seq_obj->id($temp_id);
	$seq_obj->seq($seq);
	$seq_obj->desc($temp_desc);
	$seq_out->write_seq($seq_obj);
	$count++;
}

#close IN;
#close OUT;


















