#!/usr/bin/perl -w
use strict;
use Getopt::Long;

# Copyright (c) 2025 Jinrong Huang <huangjinrong@genomics.cn>

my ($barcode2slide,$annotation,$input,$outdir,$suffix,$help);

GetOptions(
	"barcode2slide=s" => \$barcode2slide,
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

unless ($barcode2slide && $annotation && $input && $outdir) {
	print "\nError: Required parameters missing.\n\n";
	usage(1);
}

foreach my $file ($barcode2slide, $annotation, $input) {
	unless (-e $file) {
	die "Error: File '$file' does not exist.\n";
	}
}

`mkdir -p $outdir` unless (-e $outdir);

my %barcode2slide;
if($barcode2slide=~/\.gz$/){open IN,"gunzip -cd <$barcode2slide|" or die $!;}else{open IN,"<$barcode2slide" or die $!;}
while (<IN>){
	chomp;
	# AAACAACGAATAGTTC        17      1
	my ($barcode,$x,$y)=split /\t/;
	$barcode2slide{$barcode}="$x,$y";
}
close IN;

my %annotation;
if ($annotation=~/\.gz$/){open IN,"gunzip -cd <$annotation|" or die $!;}else{open IN,"<$annotation" or die $!;}
while (<IN>){
	chomp;
	# AAACAAGTATCTCCCA        Olfactory_area
	my ($barcode,$info)=split /\t/;
	$annotation{$barcode2slide{$barcode}}=$info if (exists $barcode2slide{$barcode});
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
	# 1,19    chr1    63142986        A       2       0
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
print OT "Annotation\tCoverage\tEdited\tREI\n";
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
	--barcode2slide FILE   Barcode to slide coordinates mapping file
	--annotation FILE      Barcode to annotation
	--input FILE           Input file
	--outdir DIR           Output directory
	--suffix STR           Input file suffix [default: RE.gz]
        --help, -h             Display this help message

Description:
	Assign spatial annotations to RNA editing sites and calculate RNA editing index

Examples:
        Basic usage:
	$0 --barcode2slide barcode_coordinates.tsv --annotation tissue_annotations.tsv --input RNA_editing_sites.REs.gz --outdir results --suffix REs.gz

        Display help:
        $0 --help

USAGE

        exit $exit_code;
}
