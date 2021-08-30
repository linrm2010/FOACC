#! /usr/bin/perl -w
################################################################################
# comparison_Seurat-SC3-DIMMSC.type.select.step01.pl version 1.0
# Copyright (C) Runmao Lin
# Contact (E-mail): linrunmao@caas.cn
# 
# This program is provided under the terms of a personal license to the recipient and 
# may only be used for the recipient's own research at an academic insititution.
# 
# For using this program in a company or for commercial purposes, a commercial license 
# is required.
# 
# The purpose of this program is used to select cell IDs from 'comparison_Seurat-SC3-DIMMSC.type.txt', 
# which is generated after running 'comparison_Seurat-SC3-DIMMSC.pl'.
################################################################################

use strict;
use warnings;
use FindBin qw($Bin);
use Getopt::Long;

my $comparison_Seurat_SC3_DIMMSC_type; # for file 'comparison_Seurat-SC3-DIMMSC.type.txt'
my $comparison_Seurat_SC3_DIMMSC; # for file 'comparison_Seurat-SC3-DIMMSC.txt'
my $output; # the output file
my $help;
my $version=sprintf("%.1f",1.0);

GetOptions
(
	"comparison_Seurat_SC3_DIMMSC=s" => \$comparison_Seurat_SC3_DIMMSC,           # string
	"comparison_Seurat_SC3_DIMMSC_type=s" => \$comparison_Seurat_SC3_DIMMSC_type, # string
	"output=s" => \$output,                                                       # string
	"help" => \$help                                                              # flag
);

################################################################################
# usage
############################# usage begin ######################################
my $usage= << "USAGE";

Example: perl $0 -comparison_Seurat_SC3_DIMMSC  comparison_Seurat-SC3-DIMMSC.txt  -comparison_Seurat_SC3_DIMMSC_type  comparison_Seurat-SC3-DIMMSC.type.txt  -output  comparison_Seurat-SC3-DIMMSC.type.select.txt 
version: 1.0
Options:
	-comparison_Seurat_SC3_DIMMSC <file>       the 'comparison_Seurat-SC3-DIMMSC.txt' file 
	-comparison_Seurat_SC3_DIMMSC_type <file>  the 'comparison_Seurat-SC3-DIMMSC.type.txt' file 
	-output <file>                             such as 'comparison_Seurat-SC3-DIMMSC.type.select.txt' 
	-help                                      print help information

USAGE
############################## usage end #######################################

if ($help || !(defined $comparison_Seurat_SC3_DIMMSC) || !(defined $comparison_Seurat_SC3_DIMMSC_type) || !(defined $output))
{
	print $usage;
	exit;
}

# Seurat	SC3	DIMMSC	CellNumber
# 0	1	10	100

my $select_type="$Bin/comparison_Seurat-SC3-DIMMSC.type.select.step02.pl";

my @type_ID;
my @type_num;
open F,"$comparison_Seurat_SC3_DIMMSC_type" || die "Cannot open the file '$comparison_Seurat_SC3_DIMMSC_type'.\n";
while(<F>)
{
	chomp;
	my @sp=split(/\t/,$_);
	if($_ ne "" && !($_=~/^Seurat/))
	{
		if(!(defined $type_num[0]))
		{
			push(@type_ID,"$sp[0] $sp[1] $sp[2]");
			push(@type_num,$sp[3]);
		}
		else
		{
			my $tag=0;
			for(my $i=0;$i<@type_num;$i++)
			{
				if($type_num[$i]<$sp[3])
				{
					for(my $j=(scalar @type_num)-1;$j>=$i;$j--)
					{
						$type_ID[$j+1]=$type_ID[$j];
						$type_num[$j+1]=$type_num[$j];
					}
					$type_ID[$i]="$sp[0] $sp[1] $sp[2]";
					$type_num[$i]=$sp[3];
					$tag=1;
					last;
				}
			}
			if($tag==0)
			{
				push(@type_ID,"$sp[0] $sp[1] $sp[2]");
				push(@type_num,$sp[3]);
			}
		}
	}
}
close F;

my %a_ID;
my %b_ID;
my %c_ID;

open OUT,">$output";
for(my $i=0;$i<@type_ID;$i++)
{
	my @sp=split(/\s+/,$type_ID[$i]);
	if(!(defined $a_ID{$sp[0]}) && !(defined $b_ID{$sp[1]}) && !(defined $c_ID{$sp[2]}))
	{
		system(" perl $select_type -comparison_Seurat_SC3_DIMMSC $comparison_Seurat_SC3_DIMMSC -seurat_cluster_ID $sp[0] -sc3_cluster_ID $sp[1] -dimmsc_cluster_ID $sp[2] -output select_ID  ");
		open S,"select_ID" || die "Cannot open the file 'select_ID'.\n";
		while(<S>)
		{
			chomp;
			if($_ ne "")
			{
				print OUT "$_\n";
			}
		}
		close S;
		system(" rm select_ID ");
		$a_ID{$sp[0]}=1;
#		$b_ID{$sp[0]}=1;
#		$c_ID{$sp[0]}=1;
	}
}
close OUT;

__END__
