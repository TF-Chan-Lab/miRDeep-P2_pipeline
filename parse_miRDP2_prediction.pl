#!/bin/perl
#Author: ALin
#Purpose: Parse miRDP2 prediction file and merge the files indicated in a list.
#Usage: parse_miRDP2_prediction.pl <input list> <prefix>

if(@ARGV != 2){
	print "Usage: parse_miRDP2_prediction.pl <input list> <prefix>\n";
	exit;
}

use Bio::SeqIO;
use Bio::Seq;
use Switch;
use strict;

my $in = shift;
my $prefix = shift;
my $mature = $prefix . "_mature.fa";
my $arms = $prefix . "_arms.txt";

open(IN, "<", $in) or die "Cannot open $in!\n";
open(ARMS, ">$arms") or die "Cannot open $arms!\n";
my $mature_out = Bio::SeqIO->new(-file => ">$mature", -format => 'fasta', );

my %mature_seq = ();
my %mature_pre = ();
my %mature_pri = ();
my %mature_locus = ();
my $mature_count = 1;
my $temp_pos = "";
my $temp_locus = "";	
my $temp_5pi_seq = "";
my $temp_3pi_seq = "";
my $temp_pre_seq = "";
my $temp_pri_seq = "";

my %pair = ();

while(<IN>){
	chomp;
	my $file = $_;
	open(FILE, "<$file") or die "Cannot open $file!\n";
	while(<FILE>){
		chomp;
		my $line = $_;
		my @line = split('\t', $line);
		$line[0] =~ s/\s//g;
		$line[1] =~ s/\s//g;
	

		switch ($line[0]){
		
			#case "score_star"	{$desc .= "SCORE_STAR=$line[1] ";}
			#case "score_randfold"	{$desc .= "SCORE_RANDFOLD=$line[1] ";}
			#case "score_mfe"	{$desc .= "SCORE_MEF=$line[1] ";}
			#case "score_freq"	{$desc .= "SCORE_FREQ=$line[1] ";}
			#case "score"	{$desc .= "SCORE=$line[1] ";}
			case "mature_arm"	{$temp_pos = $line[1];}
			case "mature_seq"	{
							if($temp_pos eq "first"){
								$temp_5pi_seq = $line[1];
							}
							else{
								$temp_3pi_seq = $line[1];
							}
						}
			case "pre_seq"	{$temp_pre_seq = $line[1];}
			case "pri_id"	{$temp_locus = $line[1];}
			case "pri_seq"	{$temp_pri_seq = $line[1];}
			case "star_seq"	{
						if($temp_pos eq "first"){
							$temp_3pi_seq = $line[1];
						}
						else{
							$temp_5pi_seq = $line[1];
						}
						##################################
						unless(exists $mature_seq{$temp_5pi_seq}){
							$mature_seq{$temp_5pi_seq} = "";
						}
						unless(exists $mature_seq{$temp_3pi_seq}){
							$mature_seq{$temp_3pi_seq} = "";
						}
						#Precursors
						unless(exists $mature_pre{$temp_5pi_seq}){
							%{$mature_pre{$temp_5pi_seq}} = ();
						}
						$mature_pre{$temp_5pi_seq}{$temp_pre_seq} = 1;
						unless(exists $mature_pre{$temp_3pi_seq}){
							%{$mature_pre{$temp_3pi_seq}} = ();
						}
						$mature_pre{$temp_3pi_seq}{$temp_pre_seq} = 1;
						#Primary
						unless(exists $mature_pri{$temp_5pi_seq}){
                                                        %{$mature_pri{$temp_5pi_seq}} = ();
                                                }
                                                $mature_pri{$temp_5pi_seq}{$temp_pri_seq} = 1;
                                                unless(exists $mature_pri{$temp_3pi_seq}){
                                                        %{$mature_pri{$temp_3pi_seq}} = ();
                                                }
                                                $mature_pri{$temp_3pi_seq}{$temp_pri_seq} = 1;
						#Locus
						unless(exists $mature_locus{$temp_5pi_seq}){
							%{$mature_locus{$temp_5pi_seq}} = ();
						}
						$mature_locus{$temp_5pi_seq}{$temp_locus} = 1;
						unless(exists $mature_locus{$temp_3pi_seq}){
							%{$mature_locus{$temp_3pi_seq}} = ();
						}
						$mature_locus{$temp_3pi_seq}{$temp_locus} = 1;
						#Pair
						unless(exists $pair{$temp_5pi_seq}){
							%{$pair{$temp_5pi_seq}} = ();
						}
						unless(exists $pair{$temp_5pi_seq}{$temp_3pi_seq}){
							$pair{$temp_5pi_seq}{$temp_3pi_seq} = "";
						}
						$pair{$temp_5pi_seq}{$temp_3pi_seq} .= "$temp_locus:$temp_pre_seq:$temp_pri_seq;";
						###################################
						$temp_pos = "";
						$temp_locus = "";
						$temp_5pi_seq = "";
						$temp_3pi_seq = "";
						$temp_pre_seq = "";
						$temp_pri_seq = "";
					}
		
		}

	}
	close FILE;
}

close IN;

foreach my $mseq (keys %mature_seq){
	my $seq_obj = Bio::Seq->new();
	my $id = $prefix . "_mature_" . $mature_count;
	my $desc = "COOR=";
	foreach my $coor (keys %{$mature_locus{$mseq}}){
		$desc .= "$coor,";
	}
	$desc .= "; Precursor=";
	foreach my $pre (keys %{$mature_pre{$mseq}}){
		$desc .= "$pre,";
	}
	$desc .= "; Primary=";
	foreach my $pri (keys %{$mature_pri{$mseq}}){
		$desc .= "$pri,";
	}
	$mature_seq{$mseq} = $id;
	$seq_obj->id($id);
	$seq_obj->seq($mseq);
	$seq_obj->desc($desc);
	$mature_out->write_seq($seq_obj);
	$mature_count++;
}

my $new_line_flag = 0;

foreach my $five_pi_seq (keys %pair){
	foreach my $three_pi_seq (keys %{$pair{$five_pi_seq}}){
		if($new_line_flag == 0){
				$new_line_flag = 1;
		}
		else{
			print ARMS "\n";
		}
		print ARMS "$mature_seq{$five_pi_seq}\t$mature_seq{$three_pi_seq}\t$pair{$five_pi_seq}{$three_pi_seq}";
	}
}

close ARMS;




