#! /usr/bin/perl -w
################################################################################
# cell_cluster_change.pl version 1.0
# Copyright (C) Runmao Lin
# Contact (E-mail): linrunmao@caas.cn
# 
# This program is provided under the terms of a personal license to the recipient and 
# may only be used for the recipient's own research at an academic insititution.
# 
# For using this program in a company or for commercial purposes, a commercial license 
# is required.
# 
# The purpose of this program is to detect the change of cell clusters for each cell.
################################################################################

use strict;
use warnings;
use Getopt::Long;

my $ori_cell_cluster_list; # list of cell clusters 
my $new_cell_cluster_list; # list of cell clusters 
my $output; # the output file
my $help;
my $version=sprintf("%.1f",1.0);

GetOptions
(
	"ori_cell_cluster_list=s" => \$ori_cell_cluster_list,         # string
	"new_cell_cluster_list=s" => \$new_cell_cluster_list,         # string
	"output=s" => \$output,                                       # string
	"help" => \$help                                                       # flag
);

################################################################################
# usage
############################# usage begin ######################################
my $usage= << "USAGE";

Example: perl $0 -ori_cell_cluster_list  ori_cell_cluster.list  -new_cell_cluster_list  new_cell_cluster.list  -output  cell_cluster_change.txt 
version: 1.0
Options:
	-ori_cell_cluster_list <file>         such as 'ori_cell_cluster.list'
	-new_cell_cluster_list <file>         such as 'new_cell_cluster.list'
	-output <file>                        the output file, such as 'cell_cluster_change.txt'
	-help                                 print help information

USAGE
############################## usage end #######################################

if ($help || !(defined $ori_cell_cluster_list) || !(defined $new_cell_cluster_list) || !(defined $output))
{
	print $usage;
	exit;
}

## 'seurat_cell.list'
# C0	/10t/linrm/single-cell/Arabidopsis_thaliana_sample3-4/200107U_S1_cellranger_v3/map_AthS4_AthS3_SeuratV4/test/cell.C0.txt

## cell.C0.txt
# # x
# # 1	AAACCCATCCAGTGTA-1

my %ori;
my %new;

open A,"$ori_cell_cluster_list" || die "Cannot open the file '$ori_cell_cluster_list'.\n";
while(<A>)
{
	chomp;
	my @sp=split(/\t/,$_);
	if($_ ne "")
	{
		open F,"$sp[1]" || die "Cannot open the file '$sp[1]'.\n";
		while(<F>)
		{
			chomp;
			my @csp=split(/\t/,$_);
			if($_ ne "" && !($_=~/^x/))
			{
				$ori{$csp[1]}=$sp[0];
			}
		}
		close F;
	}
}
close A;

open B,"$new_cell_cluster_list" || die "Cannot open the file '$new_cell_cluster_list'.\n";
while(<B>)
{
	chomp;
	my @sp=split(/\t/,$_);
	if($_ ne "")
	{
		open F,"$sp[1]" || die "Cannot open the file '$sp[1]'.\n";
		while(<F>)
		{
			chomp;
			my @csp=split(/\t/,$_);
			if($_ ne "" && !($_=~/^x/))
			{
				$new{$csp[1]}=$sp[0];
			}
		}
		close F;
	}
}
close B;

open OUT,">$output";
print OUT "CellID\tOriCluster\tNewCluster\n";
foreach my $id(sort keys %ori)
{
	if(!(defined $new{$id}))
	{
		$new{$id}=" --";
	}
	print OUT "$id\t$ori{$id}\t$new{$id}\n";
}
close OUT;

__END__
