#!/usr/bin/perl -w
use strict;
use Getopt::Long;

# Copyright (c) 2025 Jinrong Huang <huangjinrong@genomics.cn>

my ($inBam,$outBam,$samtools,$help);

GetOptions(
	"inBam=s"    => \$inBam,
	"outBam=s"   => \$outBam,
	"samtools=s" => \$samtools,
	"help|h"     => \$help,
) or usage(1);

if ($help){
	usage(0);
}

unless ($inBam && $outBam && $samtools) {
	print "\nError: Input BAM, output BAM, and samtools path must be specified.\n\n";
	usage(1);
}

unless (-e $inBam) {
	die "\nError: Input BAM file '$inBam' does not exist.\n";
}

unless (-x $samtools) {
	die "\nError: Samtools executable '$samtools' not found or not executable.\n";
}

my %hash=();
my $PreviousSEQ;

open OT,"| $samtools view -b -S ->$outBam" or die $!;
open IN,"$samtools view -h $inBam|" or die $!;
while (<IN>){
	chomp;
	if (/^@/){
		print OT "$_\n";
		next;
	}
	next unless (/\s+NH:i:1\s+/ || /\s+NH:i:1$/);

	my ($FLAG,$SEQ)=(split /\t/)[1,9];

	next unless $FLAG < 256;
	# 0x100   256  SECONDARY      secondary alignment
	# 0x200   512  QCFAIL         not passing quality controls or other filters
	# 0x400  1024  DUP            PCR or optical duplicate
	# 0x800  2048  SUPPLEMENTARY  supplementary alignment

	my $UR=$1 if (/\s+(UR:Z:\w+)\s+/);
	my $Cx=$1 if (/\s+(Cx:i:\d+)\s+/);
	my $Cy=$1 if (/\s+(Cy:i:\d+)\s+/);

	my $seqInfo="$SEQ\t$UR\t$Cx\t$Cy";

	if (!defined $PreviousSEQ){
		$PreviousSEQ=$SEQ;
		$hash{$seqInfo}++;
		print OT "$_\n";
	}
	else {
		if ($SEQ eq $PreviousSEQ){
			if (exists $hash{$seqInfo}){
				$hash{$seqInfo}++;
				next;
			}
			else {
				$hash{$seqInfo}++;
				print OT "$_\n";
			}
		}
		else {
			%hash=();
			$PreviousSEQ=$SEQ;
			$hash{$seqInfo}++;
			print OT "$_\n";
		}
	}
}
close IN;
close OT;

print "Processing completed successfully.\n";

sub usage {
    my ($exit_code) = @_;
    
    print <<"USAGE";

Usage: $0 [options]

Options:
	--inBam FILE        Input BAM file (required)
	--outBam FILE       Output BAM file (required)
	--samtools PATH     Path to samtools executable (required)
	--help, -h          Display this help message

Description:
	This script filters Stereo-seq BAM files to retain only unique alignments (NH:i:1) while removing secondary, duplicate, and supplementary alignments.
	It also filters based on sequence uniqueness and preserve spatial coordinates (x,y) and UMI information.

Examples:
	Basic usage:
	$0 --inBam input.bam --outBam output.bam --samtools /usr/bin/samtools

	Display help:
	$0 --help

USAGE

	exit $exit_code;
}
