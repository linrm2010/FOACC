#! /usr/bin/perl -w
################################################################################
# ID_from_cell_list.pl version 1.0
# Copyright (C) Runmao Lin
# Contact (E-mail): linrunmao@caas.cn
# 
# This program is provided under the terms of a personal license to the recipient and 
# may only be used for the recipient's own research at an academic insititution.
# 
# For using this program in a company or for commercial purposes, a commercial license 
# is required.
# 
# The purpose of this program is obtain cell IDs from list file.
################################################################################

use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;

my $cell_list; # the file includes list of cell ID files
my $output; # the output file
my $help;
my $version=sprintf("%.1f",1.0);

GetOptions
(
	"cell_list=s" => \$cell_list,                            # string
	"output=s" => \$output,                                  # string
	"help" => \$help                                         # flag
);

################################################################################
# usage
############################# usage begin ######################################
my $usage= << "USAGE";

Example: perl $0 -cell_list  cell_ID_file.list  -output  cell.id 
version: 1.0
Options:
	-cell_list <file>                     such as 'cell_ID_file.list'
	-output <file>                        such as 'cell.id'
	-help                                 print help information

USAGE
############################## usage end #######################################

if ($help || !(defined $cell_list) || !(defined $output))
{
	print $usage;
	exit;
}

# the content of cell_ID_file.list
# A01	/Ath_roots_Ryu/step01/AthRoot.step01.Ath01_ID.txt
# A02	/Ath_roots_Ryu/step01/AthRoot.step01.Ath02_ID.txt

# the content of AthRoot.step01.Ath01_ID.txt
# x
# 1	AAACCTGAGACAGACC-1
# 2	AAACCTGAGATCCGAG-1

my %cell_ID;
open L,"$cell_list" || die "Cannot open the file '$cell_list'.\n";
while(<L>)
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
			if($_ ne "" && (scalar @csp)==2)
			{
				$cell_ID{$csp[1]}=1;
			}
		}
		close F;
	}
}
close L;

open OUT,">$output" || die "Cannot open the file '$output'.\n";
foreach my $id(sort keys %cell_ID)
{
	print OUT "$id\n";
}
close OUT;

__END__
