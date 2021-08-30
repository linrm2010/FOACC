#! /usr/bin/perl -w
################################################################################
# change_10x_chromium_lib_sequenced_data_format.pl version 1.0
# Copyright (C) Runmao Lin
# Contact (E-mail): linrunmao@caas.cn
# 
# This program is provided under the terms of a personal license to the recipient and 
# may only be used for the recipient's own research at an academic insititution.
# 
# For using this program in a company or for commercial purposes, a commercial license 
# is required.
# 
# The purpose of this program is to change the sequenced data for running Cell Ranger v2.1.1.
################################################################################

use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;

my $read1fq; # 1.fq 
my $read2fq; # 2.fq 
my $output_prefix; # the prefix for output files
my $help;
my $version=sprintf("%.1f",1.0);

GetOptions
(
	"read1fq=s" => \$read1fq,                  # string
	"read2fq=s" => \$read2fq,                  # string
	"output_prefix=s" => \$output_prefix,      # string
	"help" => \$help                           # flag
);

################################################################################
# usage
############################# usage begin ######################################
my $usage= << "USAGE";

Example: perl $0 -read1fq  1.fq.gz  -read2fq  2.fq.gz  -output_prefix  YSZT_1_L101 
version: 1.0
Options:
	-read1fq <file>         such as '1.fq.gz'
	-read2fq <file>         such as '2.fq.gz'
	-output_prefix <file>   the prefix for output files, such as 'YSZT_1_L101'
	-help                   print help information

USAGE
############################## usage end #######################################

if ($help || !(defined $read1fq) || !(defined $read2fq) || !(defined $output_prefix))
{
	print $usage;
	exit;
}

my @months = qw(Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec);
my ($timesecond, $timeminute, $timehour, $timedaymonths, $timemonth, $timeyear, $timedayweek, $timedayYear, $timedaylightsavings) = localtime();
$timeyear+=1900;

print "[$timehour\:$timeminute\:$timesecond\, $months[$timemonth] $timedaymonths\, $timeyear] ";
print "Start to run the \'$0\' program ...\n";

################################################################################
# Reading fastq files
################################################################################
print "Reading fastq files of '$read1fq' and '$read2fq' ...\n";

open OA,">$output_prefix\_R1_001.fastq";
open OB,">$output_prefix\_R2_001.fastq";
if($read1fq=~/\.gz$/)
{
	open A,"gzip -dc $read1fq | " || die "Cannot open the file '$read1fq'.\n";
}
else
{
	open A,"$read1fq" || die "Cannot open the file '$read1fq'.\n";
}
if($read2fq=~/\.gz$/)
{
	open B,"gzip -dc $read2fq | " || die "Cannot open the file '$read2fq'.\n";
}
else
{
	open B,"$read2fq" || die "Cannot open the file '$read2fq'.\n";
}
while(<A>)
{
	chomp;
	if($_=~/^@/)
	{
		print OA "$_\n";
		my $seq=<A>;
		chomp $seq;
		my $useseq=substr($seq,0,26);
		<A>;
		my $qual=<A>;
		chomp $qual;
		my $usequal=substr($qual,0,26);
		print OA "$useseq\n\+\n$usequal\n";
		my $Bid=<B>;
		chomp $Bid;
		$seq=<B>;
		chomp $seq;
		$useseq=substr($seq,0,98);
		<B>;
		$qual=<B>;
		chomp $qual;
		$usequal=substr($qual,0,98);
		print OB "$Bid\n$useseq\n\+\n$usequal\n";
	}
}
close A;
close B;
close OA;
close OB;

################################################################################
# compress output files
################################################################################
print "Compressing files of '$output_prefix\_R1_001.fastq' and '$output_prefix\_R2_001.fastq' ...\n";

system(" gzip $output_prefix\_R1_001.fastq ");
system(" gzip $output_prefix\_R2_001.fastq ");

($timesecond, $timeminute, $timehour, $timedaymonths, $timemonth, $timeyear, $timedayweek, $timedayYear, $timedaylightsavings) = localtime();
$timeyear+=1900;
print "[$timehour\:$timeminute\:$timesecond\, $months[$timemonth] $timedaymonths\, $timeyear] ";
print "End of running the \'$0\' program.\n";

__END__
