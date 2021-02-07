#!/bin/perl
#Author: ALin
#Purpose: To count reads in a bam file that map to certain reference with appropriate mapping quality.
#Usage: bam2ref_counts.pl [option]
#		-bam <String> Input bam file
#		-aln <String> Output aligned reads names in a list
#		-un <String> Output unaligned reads names in a list
#		-q <Integer> Minimal mapping quality (Inclusive; default: 0)
#		-p <Float> Minimal percentage of mismatch. (Exclusive; default: off)
#		-m <Integer> Minimal number of mismatch. (Exclusive; default: off)
#		-f <String> Fasta file of the reference for exact length match between the read and reference

use strict;
use Getopt::Long;
use Bio::SeqIO;
use Bio::Seq;

my $in = "";
my $aln;
my $un;
my $mapq = 0;
my $percentage;
my $mis;
my $fasta;
my $help = 0;


GetOptions(
	'bam=s'	=>	\$in,
	'aln=s' =>	\$aln,
        'un=s'   =>      \$un,
        'q=i'   =>      \$mapq,
	'p=f'	=>	\$percentage,
	'm=i' 	=>	\$mis,
	'f=s'	=>	\$fasta,
        'h!'    =>      \$help,
);

if($help){
        print_usage();
        exit;
}

unless($in){
	print_usage();
	exit;
}

if($percentage && $mis){
	print "-p and -m are exclusive!\n";
	print_usage();
	exit;
}

open(BAM, "samtools view $in |") or die "Cannot open $in!\n";

if($aln){
	open(ALN, ">$aln") or die "Cannot create $aln!\n";
}

if($un){
	open(UN, ">$un") or die "Cannot create $un!\n";
}

my %ref_length = ();

my %ref_counts = ();
$ref_counts{"Unaligned"} = 0;
$ref_counts{"Non-unique"} = 0;
my $aln_new_line_flag = 0;
my $un_new_line_flag = 0;

if($fasta){
	my $fasta_in = Bio::SeqIO->new(-file => "$fasta", -format => 'fasta',);
	while(my $seq_obj = $fasta_in->next_seq){
		my $temp_id = $seq_obj->id;
		my $temp_length = $seq_obj->length;
		$ref_length{$temp_id} = $temp_length;
		if(exists $ref_counts{$temp_id}){
			print "Non-unique reference ID!\n";
			exit;
		}
		else{
			$ref_counts{$temp_id} = 0;
		}
	}
}

my %aligned = ();
my %unaligned = ();

OUTTER:while(<BAM>){
	chomp;
	my $line = $_;
	my @line = split('\t', $line);
	#print "$line[0]\n";
	if($line[4] < $mapq || $line[2] eq "*"){
		$unaligned{$line[0]} = 1;
	}
	else{
		#unless(exists $ref_counts{$line[2]}){
		#	$ref_counts{$line[2]} = 0;
		#}
		my @cigar = split('', $line[5]);
		my $num = @cigar;
		#foreach my $tag (@cigar){
		#	print "$tag\t";
		#}
		#print "|\t";
		my $length = length($line[9]);
		$line =~ /NM:i:(\d+)/;
		my $mismatch = $1;
		my $match = $length - $mismatch;
#		my $match = 0;
#		for(my $i = 1; $i < $num; $i++){
#			if($cigar[$i] eq "M"){
#				my $temp_num = 0;
#				BACK: for(my $j = $i - 1; $j >= 0; $j--){
#					if($cigar[$j] =~ /\d/){
#						$temp_num += $cigar[$j] * (10 ** ($i - $j - 1));
#						#print "$temp_num\t";
#					}
#					else{
#						last BACK;
#					}
#				}
#				$match += $temp_num;
#				print "$match\n";
#			}
#		}
		if($fasta && $match != $ref_length{$line[2]}){
			$unaligned{$line[0]} = 1;
			next OUTTER;
		}
		my $mis_per = $mismatch / $length;
		if((($percentage) && ($mis_per > $percentage)) || (($mis) && ($mismatch > $mis))){
			$unaligned{$line[0]} = 1;
		}
		else{
			unless(exists $aligned{$line[0]}){
				%{$aligned{$line[0]}} = ();
			}
			$aligned{$line[0]}{$line[2]} = 1;
		}
	}
}
close BAM;

foreach my $read (keys %aligned){
	my $num_refs = keys %{$aligned{$read}};
	if($num_refs > 1){
		$ref_counts{"Non-unique"}++;
	}
	else{
		foreach my $ref (keys %{$aligned{$read}}){
			$ref_counts{$ref}++;
		}
	}
	if($aln){
		if($aln_new_line_flag == 0){
			$aln_new_line_flag = 1;
		}
		else{
			print ALN "\n";
		}
		print ALN "$read";
	}
}

foreach my $read (keys %unaligned){
	if(exists $aligned{$read}){
		next;
	}
	$ref_counts{"Unaligned"}++;
	if($un){
		if($un_new_line_flag == 0){
			$un_new_line_flag = 1;
		}
		else{
			print UN "\n";
		}
		print UN "$read";
	}	
}


if($aln){
	close ALN;
}

if($un){
	close UN;
}


foreach my $ref (keys %ref_counts){
	print "$ref\t$ref_counts{$ref}\n";
}



sub print_usage{
	print "Usage: bam2genus_counts.pl [option]\n\t-bam <String> Input bam file\n\t-aln <String> Output aligned reads names in a list\n\t-un <String> Output unaligned reads in a list\n\t-q <Integer> Minimal mapping quality (Inclusive; default: 0)\n\t-p <Float> Minimal percentage of mismatch (Exclusive; default: off)\n\t-m <Integer> Minimal number of mismatch. (Exclusive; default: off)\n\t-f <String> Fasta file of the reference for exact length match between the read and reference\n";
}
