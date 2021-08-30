#! /usr/bin/perl -w
################################################################################
# generate_matrix.pl version 1.0
# Copyright (C) Runmao Lin
# Contact (E-mail): linrunmao@caas.cn
# 
# This program is provided under the terms of a personal license to the recipient and 
# may only be used for the recipient's own research at an academic insititution.
# 
# For using this program in a company or for commercial purposes, a commercial license 
# is required.
# 
# The purpose of this program is used to generate three files of 'matrix.mtx', 
# 'features.tsv' and 'barcodes.tsv' for scRNA-seq data analysis.
################################################################################

use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;

my $reformat_matrix_byCellID; # such as 'reformat_matrix_MG200_MC3_byCellID.xls'
my $features; # gene ID list, such as 'features.tsv'
my $output_prefix; # the prefix for output files
my $help;
my $version=sprintf("%.1f",1.0);

GetOptions
(
	"reformat_matrix_byCellID=s" => \$reformat_matrix_byCellID,            # string
	"features=s" => \$features,                                            # string
	"output_prefix=s" => \$output_prefix,                                  # string
	"help" => \$help                                                       # flag
);

################################################################################
# usage
############################# usage begin ######################################
my $usage= << "USAGE";

Example: perl $0 -reformat_matrix_byCellID  reformat_matrix_MG200_MC3_byCellID.xls  -features  features.tsv  -output_prefix  step02 
version: 1.0
Options:
	-reformat_matrix_byCellID <file>      the 'reformat_matrix_MG200_MC3_byCellID.xls' file for scRNA-seq analysis
	-features <file>                      gene ID list, such as 'features.tsv' file for scRNA-seq data
	-output_prefix <prefix>               the prefix for output files, such as 'step02'
	-help                                 print help information

USAGE
############################## usage end #######################################

if ($help || !(defined $reformat_matrix_byCellID) || !(defined $features) || !(defined $output_prefix))
{
	print $usage;
	exit;
}

my @months = qw(Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec);
my ($timesecond, $timeminute, $timehour, $timedaymonths, $timemonth, $timeyear, $timedayweek, $timedayYear, $timedaylightsavings) = localtime();
$timeyear+=1900;

print "[$timehour\:$timeminute\:$timesecond\, $months[$timemonth] $timedaymonths\, $timeyear] ";
print "Start to run the \'$0\' program ...\n";

## 'reformat_matrix_MG200_MC3_byCellID.xls' file format
# CellID	BrChl_00001	BrChl_00003
# AAACCTGCAAAGGTGC-1	0	0
# AAACCTGCACAGGAGT-1	0	0

################################################################################
# Reading file containing gene ID list
################################################################################
print "Reading '$features' ...\n";

my %genelist;
my $genenum=0;
if($features=~/\.gz$/)
{
	open L,"gzip -dc $features | " || die "Cannot open the file '$features'.\n";
}
else
{
	open L,"$features" || die "Cannot open the file '$features'.\n";
}
while(<L>)
{
	chomp;
	my @sp=split(/\s+/,$_);
	if($_ ne "")
	{
		$genelist{$sp[0]}=1;
	}
}
close L;

################################################################################
# Reading the file, such as 'reformat_matrix_MG200_MC3_byCellID.xls'
################################################################################
print "Reading '$reformat_matrix_byCellID' ...\n";

my %genepos;
my %expgene;
my %record;
my %cellnum2id;
my $cellnum=0;
my $expcount=0;
if($reformat_matrix_byCellID=~/\.gz$/)
{
	open F,"gzip -dc $reformat_matrix_byCellID | " || die "Cannot open the file '$reformat_matrix_byCellID'.\n";
}
else
{
	open F,"$reformat_matrix_byCellID" || die "Cannot open the file '$reformat_matrix_byCellID'.\n";
}
while(<F>)
{
	chomp;
	my @sp=split(/\t+/,$_);
	if($_=~/^CellID\t/)
	{
		for(my $i=1;$i<@sp;$i++)
		{
			if(defined $genelist{$sp[$i]})
			{
				$genepos{$i}=$sp[$i];
			}
		}
	}
	elsif($_ ne "")
	{
		$cellnum++;
		$cellnum2id{$cellnum}=$sp[0];
		my $expnum=0;
		my %expgenelist; # in one cell
		my %checkexpgene;
		foreach my $pos(keys %genepos)
		{
			if($sp[$pos]>0)
			{
				$expnum++;
				$checkexpgene{$pos}=1;
				$expgenelist{$pos}=$sp[$pos];
			}
		}
		if($expnum>0)
		{
			foreach my $pos(keys %expgenelist)
			{
				$record{$cellnum}{$pos}=$expgenelist{$pos};
				$expcount++;
				$expgene{$pos}=$genepos{$pos};
			}
		}
	}
}
close F;

################################################################################
# Writing results to '*_features.tsv' file
################################################################################
print "Writing results to '$output_prefix\_features.tsv' ...\n";

my $checkexpgenenum=0;
my %genenum2id;
open OB,">$output_prefix\_features.tsv";
foreach my $pos(sort {$a<=>$b} keys %expgene)
{
	$checkexpgenenum++;
	$genenum2id{$pos}=$checkexpgenenum;
	print OB "$expgene{$pos}\t$expgene{$pos}\n";
}
close OB;

################################################################################
# Writing results to '*_matrix.mtx' file
################################################################################
print "Writing results to '$output_prefix\_matrix.mtx' ...\n";

open OA,">$output_prefix\_matrix.mtx";
print OA "%%MatrixMarket matrix coordinate integer general\n";
print OA "%\n";
print OA "$checkexpgenenum $cellnum $expcount\n";
foreach my $cid (sort {$a<=>$b} keys %record)
{
	foreach my $nid (sort {$b<=>$a} keys %{$record{$cid}})
	{
		print OA "$genenum2id{$nid} $cid $record{$cid}{$nid}\n";
	}
}
close OA;

################################################################################
# Writing results to '*_barcodes.tsv' file
################################################################################
print "Writing results to '$output_prefix\_barcodes.tsv' ...\n";

open OC,">$output_prefix\_barcodes.tsv";
foreach my $cid(sort {$a<=>$b} keys %cellnum2id)
{
	print OC "$cellnum2id{$cid}\n";
}
close OC;

($timesecond, $timeminute, $timehour, $timedaymonths, $timemonth, $timeyear, $timedayweek, $timedayYear, $timedaylightsavings) = localtime();
$timeyear+=1900;
print "[$timehour\:$timeminute\:$timesecond\, $months[$timemonth] $timedaymonths\, $timeyear] ";
print "End of running the \'$0\' program.\n";

__END__
