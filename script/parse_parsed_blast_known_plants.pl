#!/bin/perl
#Author: ALin
#Purpose: To parsed a parsed blast result and get the entries of identical query sequence. Criteria for a known miRNA are as following: 1. plus strand; 2. exact match. Criteria for a variant of a known miRNA are as following: 1. plus strand;  2. No indel; 3. number of mismatch <= 2.
#Usage: perl parse_parsed_blast_known_plant.pl <in parsed blast file> <prefix>

use strict;

if(@ARGV != 2){
        print "Usage: perl parse_parsed_blast_known_plant.pl <in parsed blast file> <prefix>\n";
        exit;
}

my $in = shift @ARGV;
my $prefix = shift @ARGV;
my $known = $prefix . "_known.txt";
my $variant = $prefix . "_variant.txt";
my $new_line_flag_known = 0;
my $new_line_flag_variant = 0;
my %ids = ();

#my $hit_cov_cutoff = 0.9;
my $max_num_mismatch = 2;

open(IN, "<", $in) or die "Cannot open $in!\n";
open(KNOWN, ">", $known) or die "Cannot create $known!\n";
open(VAR, ">", $variant) or die "Cannot create $variant!\n";

while(<IN>){
	chomp;
	my $strand_flag = 0;
	my $exact_match_flag = 0;
	#my $hit_cov_flag = 0;
	my $indel_flag = 0;
	#my $seed_flag = 0;
	my $mismatch_flag = 0;

	my $line = $_;
	my @line = split('\t', $line);
	
	
	#Criterion 1: strand
	if($line[15] == 1){
		$strand_flag = 1;
	}
	else{
		next;
	}

	#Criterion 2: exact match

	if(($line[18] == 1) && ($line[19] == 1) && ($line[20] == 100)){
		$exact_match_flag = 1;
	}
	
	if(($strand_flag == 1) && ($exact_match_flag  == 1)){
		if($new_line_flag_known == 0){
			$new_line_flag_known = 1;
		}
		else{
			print KNOWN "\n";
		}
		print KNOWN $line;
		next;
	}
	
	#Criterion 3: indel
	if(($line[22] =~ /\-/) || ($line[23] =~ /\-/)){
		next;
	}
	else{
		$indel_flag = 1;
	}
	
	#Criterion 4: mismatch
	my $length_aligned = $line[14] - $line[13] + 1;
	my $num_mismatch = $line[2] - $length_aligned;
	for(my $i = 0; $i < $length_aligned; $i++){
		my $temp_query_nt = substr($line[22], $i, 1);
		my $temp_hit_nt = substr($line[23], $i, 1);
		if($temp_query_nt ne $temp_hit_nt){
			$num_mismatch++;
		}
	}
	if($num_mismatch <= $max_num_mismatch){
		$mismatch_flag = 1;
	}
	else{
		next;
	}
	if(($strand_flag == 1) && ($exact_match_flag == 0) && ($indel_flag == 1) && ($mismatch_flag == 1)){
		if($new_line_flag_variant == 0){
			$new_line_flag_variant = 1;
		}
		else{
			print VAR "\n";
		}
		print VAR "$line";
	}
}

#close IN;

#foreach my $id (keys %ids){
#	if($new_line_flag == 0){
#		$new_line_flag = 1;
#	}
#	else{
#		print OUT "\n";
#	}
#	print OUT $id;
#}


close KNOWN;
close VAR;







