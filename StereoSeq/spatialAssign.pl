#!/usr/bin/perl -w
use strict;
use Getopt::Long;

# Copyright (c) 2025 Jinrong Huang <huangjinrong@genomics.cn>

my ($annotation,$input,$outdir,$suffix,$help);

GetOptions(
	"annotation=s"    => \$annotation,
	"input=s"	  => \$input,
	"outdir=s"        => \$outdir,
	"suffix=s"        => \$suffix,
	"help|h"          => \$help,
) or usage(1);

if ($help) {
	usage(0);
}

$suffix ||="REs.gz";

unless ($annotation && $input && $outdir) {
	print "\nError: Required parameters missing.\n\n";
	usage(1);
}

foreach my $file ($annotation, $input) {
	unless (-e $file) {
	die "Error: File '$file' does not exist.\n";
	}
}

`mkdir -p $outdir` unless (-e $outdir);

my %annotation;
if ($annotation=~/\.gz$/){open IN,"gunzip -cd <$annotation|" or die $!;}else{open IN,"<$annotation" or die $!;}
while (<IN>){
	chomp;
	# 26032,74621     14r-l2  Excitatory neurons
	my ($xy,$region,$celltype)=split /\t/;
	$annotation{$xy}="$region\t$celltype";
}
close IN;

my $name=(split /\//,$input)[-1];
$name=~s/\.$suffix$//;

my %REI;
my $output="$outdir/$name.$suffix";
if ($output=~/\.gz$/){open OT,"|gzip >$output" or die $!;}else{open OT,">$output" or die $!;}
if ($input=~/\.gz$/){open IN,"gunzip -cd <$input|" or die $!;}else{open IN,"<$input" or die $!;}
while (<IN>){
	chomp;
	# 20003,77362     chr19   55074976        A       1       0
	my ($xy,$chr,$pos,$ref,$cov,$edited)=split /\t/;
	if (exists $annotation{$xy}){
		print OT "$_\t$annotation{$xy}\n";
		$REI{$annotation{$xy}}{'cov'}+=$cov;
		$REI{$annotation{$xy}}{'edited'}+=$edited;
	}
	else {next;}
}
close IN;
close OT;

my $REI="$outdir/$name.REI.tsv";
if ($REI=~/\.gz$/){open OT,"|gzip >$REI" or die $!;}else{open OT,">$REI" or die $!;}
print OT "Region\tCelltype\tCoverage\tEdited\tREI\n";
foreach my $k(sort keys %REI){
	my $c=$REI{$k}{'cov'} || 0;
	my $e=$REI{$k}{'edited'} || 0;
	next unless $c >0;
	my $index=$e/$c;
	print OT "$k\t$c\t$e\t$index\n";
}
close OT;

print "spatialAssign completed successfully.\n";

sub usage {
    my ($exit_code) = @_;
    
    print <<"USAGE";

Usage: $0 [options]

Options:
	--annotation FILE      annotation
	--input FILE           Input file
	--outdir DIR           Output directory
	--suffix STR           Input file suffix [default: RE.gz]
        --help, -h             Display this help message

Description:
	Assign spatial annotations to RNA editing sites and calculate RNA editing index

Examples:
        Basic usage:
	$0 --annotation annotations.tsv --input RNA_editing_sites.REs.gz --outdir results --suffix REs.gz

        Display help:
        $0 --help

USAGE

        exit $exit_code;
}
