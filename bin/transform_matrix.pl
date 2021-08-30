#! /usr/bin/perl -w
################################################################################
# transform_matrix.pl version 1.0
# Copyright (C) Runmao Lin
# Contact (E-mail): linrunmao@caas.cn
# 
# This program is provided under the terms of a personal license to the recipient and 
# may only be used for the recipient's own research at an academic insititution.
# 
# For using this program in a company or for commercial purposes, a commercial license 
# is required.
# 
# The purpose of this program is used for transformation of data in a matrix.
################################################################################

use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;

my $matrix;
my $output;
my $help;
my $version=sprintf("%.1f",1.0);

GetOptions
(
	"matrix=s" => \$matrix,        # string
	"output=s" => \$output,        # string
	"help" => \$help                 # flag
);

################################################################################
# usage
############################# usage begin ######################################
my $usage= << "USAGE";

Example: perl $0 -matrix  matrix.txt  -output  change_matrix.txt 
version: 1.0
Options:
	-matrix <file>               one file contains scRNA-seq UMI values in matrix form, such as 'matrix.txt'
	-output <output file>        the output file in matrix form, such as 'change_matrix.txt'
	-help                        print help information

USAGE
############################## usage end #######################################

if ($help || !(defined $matrix) || !(defined $output))
{
	print $usage;
	exit;
}

my @months = qw(Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec);
my ($timesecond, $timeminute, $timehour, $timedaymonths, $timemonth, $timeyear, $timedayweek, $timedayYear, $timedaylightsavings) = localtime();
$timeyear+=1900;

print "[$timehour\:$timeminute\:$timesecond\, $months[$timemonth] $timedaymonths\, $timeyear] ";
print "Start to run the \'$0\' program ...\n";

# input matrix
# row_inforID	AAACCTGAGACGACGT	...
# AT1G01010	0	...

# output matrix
# row_inforID	AT1G01010	..
# AAACCTGAGACGACGT	0	..

################################################################################
# Reading 'matrix' file
################################################################################
print "Reading '$matrix' ...\n";

my $outID_infor="";
my @column_infor; # IDs for columns in the head line
my @row_infor; # IDs for rows
my %exp;
open T,"$matrix" || die "Cannot open the file '$matrix'.\n";
while(<T>)
{
	chomp;
	my @sp=split(/\s+/,$_);
	if($_=~/^row_inforID/ || $_=~/^column_inforID/ || $_=~/GeneID/ || $_=~/CellID/)
	{
		$outID_infor=$sp[0];
		for(my $i=1;$i<@sp;$i++)
		{
			push(@column_infor,$sp[$i]);
		}
	}
	elsif($_ ne "")
	{
		push(@row_infor,$sp[0]);
		for(my $i=1;$i<@sp;$i++)
		{
			$exp{$column_infor[$i-1]}{$sp[0]}=$sp[$i];
		}
	}
}
close T;

################################################################################
# Writing to the output file
################################################################################
print "Writing to the output file ...\n";
open OUT,">$output";
my $out=$outID_infor;
for(my $i=0;$i<@row_infor;$i++)
{
	$out.="\t$row_infor[$i]";
}
print OUT "$out\n";
for(my $i=0;$i<@column_infor;$i++)
{
	$out=$column_infor[$i];
	for(my $j=0;$j<@row_infor;$j++)
	{
		$out.="\t$exp{$column_infor[$i]}{$row_infor[$j]}";
	}
	print OUT "$out\n";
}
close OUT;

($timesecond, $timeminute, $timehour, $timedaymonths, $timemonth, $timeyear, $timedayweek, $timedayYear, $timedaylightsavings) = localtime();
$timeyear+=1900;
print "[$timehour\:$timeminute\:$timesecond\, $months[$timemonth] $timedaymonths\, $timeyear] ";
print "End of running the \'$0\' program.\n";

__END__
