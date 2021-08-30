#! /usr/bin/perl -w
################################################################################
# select_cells_by_cellID.pl version 1.0
# Copyright (C) Runmao Lin
# Contact (E-mail): linrunmao@caas.cn
# 
# This program is provided under the terms of a personal license to the recipient and 
# may only be used for the recipient's own research at an academic insititution.
# 
# For using this program in a company or for commercial purposes, a commercial license 
# is required.
# 
# The purpose of this program is used to select information from 
# 'reformat_MG200_MC4_byCellID.xls'.
################################################################################

use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;

my $reformat_byCellID; # file 'reformat_MG200_MC4_byCellID.xls'
my $cell_ID_list; # the file contains cell IDs
my $output; # the output file
my $help;
my $version=sprintf("%.1f",1.0);

GetOptions
(
	"reformat_byCellID=s" => \$reformat_byCellID,     # string
	"cell_ID_list=s" => \$cell_ID_list,               # string
	"output=s" => \$output,                           # string
	"help" => \$help                                                       # flag
);

################################################################################
# usage
############################# usage begin ######################################
my $usage= << "USAGE";

Example: perl $0  -reformat_byCellID  reformat_MG200_MC4_byCellID.xls  -cell_ID_list  cell_ID.list  -output  cell_ID.txt 
version: 1.0
Options:
	-reformat_byCellID <file>             the 'reformat_MG200_MC4_byCellID.xls' file for scRNA-seq data
	-cell_ID_list <file>                  list of cell IDs
	-output <file>                        the output file
	-help                                 print help information

USAGE
############################## usage end #######################################

if ($help || !(defined $reformat_byCellID) || !(defined $cell_ID_list) || !(defined $output))
{
	print $usage;
	exit;
}

my @months = qw(Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec);
my ($timesecond, $timeminute, $timehour, $timedaymonths, $timemonth, $timeyear, $timedayweek, $timedayYear, $timedaylightsavings) = localtime();
$timeyear+=1900;

print "[$timehour\:$timeminute\:$timesecond\, $months[$timemonth] $timedaymonths\, $timeyear] ";
print "Start to run the \'$0\' program ...\n";

# the 'reformat_MG200_MC4_byCellID.xls' file format
# CellID	BrChl_00001
# AAACCTGCAAAGGTGC-1	0

################################################################################
# Reading the file containing cell IDs
################################################################################
print "Reading '$cell_ID_list' ...\n";

my %cellid;
open I,"$cell_ID_list" || die "Cannot open the file '$cell_ID_list'.\n";
while(<I>)
{
	chomp;
	my @sp=split(/\s+/,$_);
	if($_ ne "")
	{
		$cellid{$sp[0]}=1;
	}
}
close I;

################################################################################
# Reading the file, such as 'reformat_MG200_MC4_byCellID.xls'
################################################################################
print "Reading '$reformat_byCellID' ...\n";

open OUT,">$output";
if($reformat_byCellID=~/\.gz$/)
{
	open F,"gzip -dc $reformat_byCellID | " || die "Cannot open the file '$reformat_byCellID'.\n";
}
else
{
	open F,"$reformat_byCellID" || die "Cannot open the file '$reformat_byCellID'.\n";
}
while(<F>)
{
	chomp;
	my @sp=split(/\t/,$_);
	if($_=~/^CellID/)
	{
		print OUT "$_\n";
	}
	elsif($_ ne "")
	{
		if(defined $cellid{$sp[0]})
		{
			print OUT "$_\n";
		}
	}
}
close F;
close OUT;

($timesecond, $timeminute, $timehour, $timedaymonths, $timemonth, $timeyear, $timedayweek, $timedayYear, $timedaylightsavings) = localtime();
$timeyear+=1900;
print "[$timehour\:$timeminute\:$timesecond\, $months[$timemonth] $timedaymonths\, $timeyear] ";
print "End of running the \'$0\' program.\n";

__END__
