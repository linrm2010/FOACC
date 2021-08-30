#! /usr/bin/perl -w
################################################################################
# reformat_matrix.pl version 1.0
# Copyright (C) Runmao Lin
# Contact (E-mail): linrunmao@caas.cn
# 
# This program is provided under the terms of a personal license to the recipient and 
# may only be used for the recipient's own research at an academic insititution.
# 
# For using this program in a company or for commercial purposes, a commercial license 
# is required.
# 
# The purpose of this program is used to integrate scRNA-seq data from 'matrix.mtx', 
# 'features.tsv' and 'barcodes.tsv'.
################################################################################

use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;

my $matrix;
my $features;
my $barcodes;
my $min_expressed_genes_one_cell; # minimum of expressed genes in one cell, such as 200
my $min_cells_expressed_one_gene; # minimum of cells contain the same gene, such as 3
my $output_prefix; # the prefix for output files
my $help;
my $version=sprintf("%.1f",1.0);

GetOptions
(
	"matrix=s" => \$matrix,                                                # string
	"features=s" => \$features,                                            # string
	"barcodes=s" => \$barcodes,                                            # string
	"min_expressed_genes_one_cell=i" => \$min_expressed_genes_one_cell,    # numeric
	"min_cells_expressed_one_gene=i" => \$min_cells_expressed_one_gene,    # numeric
	"output_prefix=s" => \$output_prefix,                                  # string
	"help" => \$help                                                       # flag
);

################################################################################
# usage
############################# usage begin ######################################
my $usage= << "USAGE";

Example: perl $0 -matrix  matrix.txt  -features  features.tsv  -barcodes  barcodes.tsv  -min_expressed_genes_one_cell  200  -min_cells_expressed_one_gene  3  -output_prefix  reformat 
version: 1.0
Options:
	-matrix <file>                        the 'matrix.mtx' file for scRNA-seq data
	-features <file>                      the 'features.tsv' file for scRNA-seq data
	-barcodes <file>                      the 'barcodes.tsv' file for scRNA-seq data
	-output_prefix <prefix>               the prefix for output files, such as 'reformat'
	-min_expressed_genes_one_cell <int>   default: 200
	-min_cells_expressed_one_gene <int>   default: 3
	-help                                 print help information

USAGE
############################## usage end #######################################

if ($help || !(defined $matrix) || !(defined $features) || !(defined $barcodes) || !(defined $min_expressed_genes_one_cell) || !(defined $min_cells_expressed_one_gene) || !(defined $output_prefix))
{
	print $usage;
	exit;
}

if(!(defined $min_expressed_genes_one_cell))
{
	$min_expressed_genes_one_cell=200;
}

if(!(defined $min_cells_expressed_one_gene))
{
	$min_cells_expressed_one_gene=3;
}

my @months = qw(Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec);
my ($timesecond, $timeminute, $timehour, $timedaymonths, $timemonth, $timeyear, $timedayweek, $timedayYear, $timedaylightsavings) = localtime();
$timeyear+=1900;

print "[$timehour\:$timeminute\:$timesecond\, $months[$timemonth] $timedaymonths\, $timeyear] ";
print "Start to run the \'$0\' program ...\n";

## 'matrix.mtx' file format
# %%MatrixMarket matrix coordinate integer general
# %
# 35386 7310 10213173   # gene number,  cell population,  lines
# 35311 1 1
# 35257 1 1
# 35253 1 1

################################################################################
# Reading 'features.tsv' file
################################################################################
print "Reading 'features.tsv' ...\n";

my %genelist;
my $genelistnum;
if($features=~/\.gz$/)
{
	open GL,"gzip -dc $features | " || die "Cannot open the file '$features'.\n";
}
else
{
	open GL,"$features" || die "Cannot open the file '$features'.\n";
}
while(<GL>)
{
	chomp;
	my @sp=split(/\s+/,$_);
	if($_ ne "")
	{
		$genelistnum++;
		$genelist{$genelistnum}=$sp[0];
	}
}
close GL;

################################################################################
# Reading 'barcodes.tsv' file
################################################################################
print "Reading 'barcodes.tsv' ...\n";

my %barcodelist;
my $barcodelistnum;
if($barcodes=~/\.gz$/)
{
	open BC,"gzip -dc $barcodes | " || die "Cannot open the file '$barcodes'.\n";
}
else
{
	open BC,"$barcodes" || die "Cannot open the file '$barcodes'.\n";
}
while(<BC>)
{
	chomp;
	my @sp=split(/\s+/,$_);
	if($_ ne "")
	{
		$barcodelistnum++;
		$barcodelist{$barcodelistnum}=$sp[0];
	}
}
close BC;

################################################################################
# Reading 'matrix.mtx' file
################################################################################
print "Reading 'matrix.mtx' ...\n";

my %gene_infor;
my %cell_infor;
my %geneid;
my %cellid;
if($matrix=~/.gz$/)
{
	open F,"gzip -dc $matrix | " || die "Cannot open the file '$matrix'.\n";
}
else
{
	open F,"$matrix" || die "Cannot open the file '$matrix'.\n";
}
while(<F>)
{
	chomp;
	if($_=~/^\%\%/)
	{
		<F>;
		<F>;
	}
	elsif($_ ne "")
	{
		my @sp=split(/\s+/,$_);
		my $origeneid=$genelist{$sp[0]}; # origeneid, i.e., related ID in 'features.tsv'
		my $oricellid=$barcodelist{$sp[1]}; # oricellid, i.e., related ID in 'barcodes.tsv'
		$gene_infor{$origeneid}{$oricellid}=$sp[2];
		$cell_infor{$oricellid}{$origeneid}=$sp[2];
		$geneid{$origeneid}=1;
		$cellid{$oricellid}=1;
	}
}
close F;

my @geneidlist;
my @cellidlist;
my $geneidlistinfor;
my $cellidlistinfor;
foreach my $id(sort keys %geneid)
{
	push(@geneidlist,$id);
	$geneidlistinfor.="$id\t"
}
$geneidlistinfor=~s/\s+$//;
foreach my $id(sort keys %cellid)
{
	push(@cellidlist,$id);
	$cellidlistinfor.="$id\t";
}
$cellidlistinfor=~s/\s+$//;

################################################################################
# Writing results to '*_byGeneID.xls' file
################################################################################
print "Writing results to '$output_prefix\_MG$min_expressed_genes_one_cell\_MC$min_cells_expressed_one_gene\_byGeneID.xls' ...\n";

open OA,">$output_prefix\_MG$min_expressed_genes_one_cell\_MC$min_cells_expressed_one_gene\_byGeneID.xls";
print OA "GeneID\t$cellidlistinfor\n";
foreach my $gid(sort keys %gene_infor)
{
	my $min=0;
	my $outline=$gid;
	for(my $i=0;$i<@cellidlist;$i++)
	{
		if(!(defined $gene_infor{$gid}{$cellidlist[$i]}))
		{
			$gene_infor{$gid}{$cellidlist[$i]}=0;
		}
		$min++ if($gene_infor{$gid}{$cellidlist[$i]}>0);
		$gene_infor{$gid}{$cellidlist[$i]}=$gene_infor{$gid}{$cellidlist[$i]};
		$outline.="\t$gene_infor{$gid}{$cellidlist[$i]}";
	}
	print OA "$outline\n" if($min>=$min_cells_expressed_one_gene);
}
close OA;

################################################################################
# Writing results to '*_byCellID.xls' file
################################################################################
print "Writing results to '$output_prefix\_MG$min_expressed_genes_one_cell\_MC$min_cells_expressed_one_gene\_byCellID.xls' ...\n";

open OB,">$output_prefix\_MG$min_expressed_genes_one_cell\_MC$min_cells_expressed_one_gene\_byCellID.xls";
print OB "CellID\t$geneidlistinfor\n";
foreach my $cid(sort keys %cell_infor)
{
	my $min=0;
	my $outline=$cid;
	for(my $i=0;$i<@geneidlist;$i++)
	{
		if(!(defined $cell_infor{$cid}{$geneidlist[$i]}))
		{
			$cell_infor{$cid}{$geneidlist[$i]}=0;
		}
		$min++ if($cell_infor{$cid}{$geneidlist[$i]}>0);
		$cell_infor{$cid}{$geneidlist[$i]}=$cell_infor{$cid}{$geneidlist[$i]};
		$outline.="\t$cell_infor{$cid}{$geneidlist[$i]}";
	}
	print OB "$outline\n" if($min>=$min_expressed_genes_one_cell);
}
close OB;

($timesecond, $timeminute, $timehour, $timedaymonths, $timemonth, $timeyear, $timedayweek, $timedayYear, $timedaylightsavings) = localtime();
$timeyear+=1900;
print "[$timehour\:$timeminute\:$timesecond\, $months[$timemonth] $timedaymonths\, $timeyear] ";
print "End of running the \'$0\' program.\n";

__END__
