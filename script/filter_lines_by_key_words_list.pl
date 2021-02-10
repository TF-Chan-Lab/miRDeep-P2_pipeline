#!/bin/perl
#Author: ALin
#Purpose: To filter lines in any files by looking at a list of words in specific column (0-based).
#Usage: filter_lines_by_key_words_list.pl <in file> <output> <list of key words> <column A> [column B]...

use strict;

if(@ARGV < 4){
	print "#Usage: filter_lines_by_key_words_list.pl <in file> <output> <list of key words> <column A (0-based)> [column B]...\n";
	exit;
}

my $in = shift @ARGV;
my $out = shift @ARGV;
my $list = shift @ARGV;
my %key_words = ();
my @columns = ();

while(@ARGV > 0){
	my $temp_column = shift @ARGV;
	print "Column $temp_column (+1) will be checked...\n";
        push(@columns, $temp_column);
}

open(IN, "<", $in) or die"Cannot open $in!\n";
open(LIST, "<", $list) or die"Cannot open $list!\n";
open(OUT, ">", $out) or die"Cannot create $out!\n";

while(<LIST>){
	chomp;
	my $key_word = $_;
	$key_words{$key_word} = 1;
}

close LIST;

my $num = (keys %key_words);
my $num_c = @columns;
print "Total entries in the list: $num\n$num_c column(s)\n";

while(<IN>){
	my $flag = 0;
	my $line = $_;
	chomp $line;
	my @line = split('\t', $line);
	foreach my $i (@columns){
		if($key_words{$line[$i]} == 1){
			#print "$key_words{$line[$i]}\n";
			$flag++;
		}
	}
	#print "$flag\n";
	if($flag == 0){
		print OUT "$line\n";
	}
}

close IN;
close OUT;


















