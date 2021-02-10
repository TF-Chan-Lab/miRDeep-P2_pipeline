#!/bin/perl
#Author: ALin
#Purpose: To parse a blast result.
#Usage: perl general_blast_parser.pl <in blast result> <output flie>

use strict;
use Bio::SearchIO; 
use Bio::DB::Fasta;
use Bio::DB::GenBank;
use Bio::PrimarySeqI;
use Bio::SeqIO::genbank;
use Bio::LocationI;
use Bio::SeqFeature::Generic;
use Bio::Seq;
use Bio::Tools::CodonTable;

if(@ARGV != 2){
	print "#Usage: perl general_blast_parser.pl <in blast result> <output flie>\n";
	exit;
}


my $in = shift @ARGV;
my $out = shift @ARGV;

open(IN, "<", $in) or die "Cannot open $in!\n";
open(OUT, ">", $out) or die "Cannot create $out!\n";

my $in_obj = new Bio::SearchIO(-format => 'blast', -file => $in);

print OUT "Query_name\tQuery_acc\tQuery_length\tHit name\tHit acc\tDescription\tHit length\tLocus\tRank\tNum hsps\tQuery start\tQuery end\tQuery strand\tHit start\tHit end\tHit strand\tFrame\tHit score\tQuery cover\tHit cover\tIdentity\tE-value\tQuery_seq\tHit_seq";

while(my $result = $in_obj->next_result){
	my $query_name = $result->query_description;
        my $query_acc = $result->query_accession;
       	my $query_length = $result->query_length;
	while(my $hit = $result->next_hit){
		my $hit_name = $hit->name;
              	my $hit_acc = $hit->accession;
          	my $hit_description = $hit->description;
       		my $hit_length = $hit->length;
       		my $locus = $hit->locus;
		my $num_hsps = $hit->num_hsps;
		while(my $hsp = $hit->next_hsp){
			my $rank = $hsp->rank;
			my $query_start = $hsp->start('query');
			my $query_end = $hsp->end('query');
			my $query_strand = $hsp->strand('query');
			my $hit_start = $hsp->start('hit');
			my $hit_end = $hsp->end('hit');
			my $hit_strand = $hsp->strand('hit');
			my $frame = ($hsp->query->frame + 1) * $query_strand;
			my $hit_score = $hsp->score;
			my $identity = $hsp->percent_identity;
			my $query_cover = ($query_end - $query_start + 1) / $query_length;
			my $hit_cover = ($hit_end - $hit_start + 1) / $hit_length;
			my $query_seq = $hsp->query_string;
			my $hit_seq = $hsp->hit_string;
			my $evalue = $hsp->evalue;
			print OUT "\n$query_name\t$query_acc\t$query_length\t$hit_name\t$hit_acc\t$hit_description\t$hit_length\t$locus\t$rank\t$num_hsps\t$query_start\t$query_end\t$query_strand\t$hit_start\t$hit_end\t$hit_strand\t$frame\t$hit_score\t$query_cover\t$hit_cover\t$identity\t$evalue\t$query_seq\t$hit_seq";
    		}
  	}
}

close IN;
close OUT;
