#! /usr/bin/perl -w
################################################################################
# comparison_Seurat-SC3-DIMMSC.type.select.step02.pl version 1.0
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
use Getopt::Long;

my $comparison_Seurat_SC3_DIMMSC; # for file 'comparison_Seurat-SC3-DIMMSC.txt'
my $seurat_cluster_ID; # ID for Seurat cluster
my $sc3_cluster_ID; # ID for SC3 cluster
my $dimmsc_cluster_ID; # ID for DIMMSC cluster
my $output; # the output file
my $help;
my $version=sprintf("%.1f",1.0);

GetOptions
(
	"comparison_Seurat_SC3_DIMMSC=s" => \$comparison_Seurat_SC3_DIMMSC,    # string
	"seurat_cluster_ID=s" => \$seurat_cluster_ID,                          # string
	"sc3_cluster_ID=s" => \$sc3_cluster_ID,                                # string
	"dimmsc_cluster_ID=s" => \$dimmsc_cluster_ID,                          # string
	"output=s" => \$output,                                                # string
	"help" => \$help                                                       # flag
);

################################################################################
# usage
############################# usage begin ######################################
my $usage= << "USAGE";

Example: perl $0 -comparison_Seurat_SC3_DIMMSC  comparison_Seurat_SC3_DIMMSC.txt  -seurat_cluster_ID  A01  -sc3_cluster_ID  1  -dimmsc_cluster_ID  1  -output  comparison_Seurat_SC3_DIMMSC.select.cell.id 
version: 1.0
Options:
	-comparison_Seurat_SC3_DIMMSC <file>   the 'comparison_Seurat_SC3_DIMMSC.txt' file 
	-seurat_cluster_ID <string>            such as 'A01' 
	-sc3_cluster_ID <string>               such as '1'
	-dimmsc_cluster_ID <prefix>            such as '1'
	-output <file>                         such as 'comparison_Seurat_SC3_DIMMSC.select.cell.id' 
	-help                                 print help information

USAGE
############################## usage end #######################################

if ($help || !(defined $comparison_Seurat_SC3_DIMMSC) || !(defined $seurat_cluster_ID) || !(defined $sc3_cluster_ID) || !(defined $dimmsc_cluster_ID) || !(defined $output))
{
	print $usage;
	exit;
}

## the contents in file 'comparison_Seurat_SC3_DIMMSC.txt'
# Cell	Seurat	SC3	DIMMSC
# AAACCCAAGGGCAATC-1	C1	2	3

open OUT,">$output";
open F,"$comparison_Seurat_SC3_DIMMSC" || die "Cannot open the file '$comparison_Seurat_SC3_DIMMSC'.\n";
while(<F>)
{
	chomp;
	my @sp=split(/\t/,$_);
	if($_ ne "" && $sp[1] eq $seurat_cluster_ID && $sp[2] eq $sc3_cluster_ID && $sp[3] eq $dimmsc_cluster_ID)
	{
		print OUT "$sp[0]\n";
	}
}
close F;

__END__
