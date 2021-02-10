#!/bin/perl
#Author: ALin
#Purpose: To combine the HTSeq output as indicated by a list. The list should contain the file name only. Option can be indicated by whether to normalize or not. Normalization is sample-based.
#Usage: perl combine_cufflinks_genes.pl [option] <input list> <output>

if(@ARGV < 2){
	print "#Usage: perl combine_cufflinks_genes.pl [option] <input list> <output>\n";
	print "Options:\n";
	print "\t-n String Integer\tInput of a normlization file, where the column (zero-based) used for normalization is indicated by the file name and an integer\n";
	exit;
}


use strict;
use Switch;

my $nor_file = "";
my $NORMALIZATION = "FALSE";
my $NOR_COL = 0;

foreach my $argv (@ARGV){
	switch($argv){
		
		case "-n"	{shift @ARGV; $nor_file = shift @ARGV; $NOR_COL = shift @ARGV; $NORMALIZATION = "TRUE";}
	}
}

my $list = shift @ARGV;
my $out = shift @ARGV;

open(LIST, "<", $list) or die "Cannot open $list!\n";
open(OUT, ">", $out) or die "Cannot create $out!\n";

my %genes = ();
my @samples = ();
my @nor_factors = ();

if($NORMALIZATION eq "TRUE"){
	open(NOR, "<", $nor_file) or die "Cannot open $nor_file!\n";
	while(<NOR>){
		my $line = $_;
		my @line = ();
		chomp $line;
		@line = split('\t', $line);
		push(@nor_factors, $line[$NOR_COL]);
	}
close NOR;
}

my $samples_count = 0;
my $num_samples = 0;

while(<LIST>){
	$num_samples++;
}
close LIST;

open(LIST, "<", $list) or die "Cannot open $list!\n";

while(<LIST>){
	my $file = $_;
	chomp $file;
	$file =~ m/(.+)\.[^\.]+/;
	push(@samples, $1);
	open(IN, "<", $file) or die "Cannot open $file!\n";
	my $flag = 1;
	while(<IN>){
		if($flag == 0){
			$flag = 1;
		}
		else{
			my $line = $_;
			my @line = ();
			chomp $line;
			@line = split('\t', $line);
			unless(exists $genes{$line[0]}){
				@{$genes{$line[0]}} = (0) x $num_samples;
			}
			if($NORMALIZATION eq "TRUE"){
				$line[1] *= $nor_factors[$samples_count];
			}
			$genes{$line[0]}[$samples_count] = $line[1];
		}
	}
	$samples_count++;
	close IN;
}

close LIST;

print OUT "Gene_id";

foreach my $i (@samples){
	print OUT "\t$i";
}

print OUT "\tMean\tMedian\tVariance\tCV2";

foreach my $gene (keys %genes){
	print OUT "\n$gene";
	foreach my $i (@{$genes{$gene}}){
		print OUT "\t$i";
	}
	my $mean = mean($genes{$gene});
	my $median = median($genes{$gene});
	my $var = variance($genes{$gene});
	my $cv2 = "NA";
	if($mean > 0){
		$cv2 = $var / ($mean ** 2);
	}
	print OUT "\t$mean\t$median\t$var\t$cv2";
}

close OUT;

sub mean{
	my @array = @{$_[0]};
	my $sum = 0;
	my $n = scalar @array;
	foreach my $i (@array){
		$sum += $i;
	}
	my $m = $sum / $n;
	#print "$m\n";
	return $m;
}

sub variance{
	my @array = @{$_[0]};
	my $ss = 0;
	my $n = @array;
	my $mean = mean(\@array);
	foreach my $i (@array){
		$ss += ($i - $mean) ** 2;
	}
	my $var = $ss / ($n - 1);
	#print "$var\n";
	return $var;
}

sub median{
	my @array = @{$_[0]};
	my @sorted_array = sort {$a <=> $b} @array;
	my $n = @sorted_array;
	if($n % 2 == 0){
		my $temp = $n / 2;
		my $m = ($sorted_array[$temp - 1] + $sorted_array[$temp]);
		#print "$m\n";
		return $m;
	}
	else{
		my $pos = ($n - 1) / 2;
		#print "$sorted_array[$pos]\n";
		return $sorted_array[$pos];
	}
}


















