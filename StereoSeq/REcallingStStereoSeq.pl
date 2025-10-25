#!/usr/bin/perl -w
use strict;
use Getopt::Long;

# Copyright (c) 2025 Jinrong Huang <huangjinrong@genomics.cn>

my ($dataset,$bam,$suffix,$outdir,$samtools,$phred,$qual_cutoff,$help);

GetOptions(
	"dataset=s"       => \$dataset,
	"bam=s"           => \$bam,
	"suffix=s"        => \$suffix,
	"outdir=s"        => \$outdir,
	"samtools=s"      => \$samtools,
	"phred=i"         => \$phred,
	"qual_cutoff=i"   => \$qual_cutoff,
	"help|h"          => \$help,
) or usage(1);

if ($help) {
	usage(0);
}

$suffix ||= "bam";
$phred ||= 33;
$qual_cutoff ||= 20;

unless ($dataset && $bam && $outdir && $samtools) {
	print "\nError: Required parameters missing.\n\n";
	usage(1);
}

unless (-e $dataset) {
	die "\nError: Dataset file '$dataset' does not exist.\n";
}

unless (-e $bam) {
	die "\nError: BAM file '$bam' does not exist.\n";
}

unless (-x $samtools) {
	die "\nError: Samtools executable '$samtools' not found or not executable.\n";
}

`mkdir -p $outdir` unless (-e $outdir);

my %hash = (
	"A" => "G",
	"T" => "C",
);

my %dataset;
if($dataset=~/\.gz$/){open IN,"gunzip -cd <$dataset|" or die $!;}else{open IN,"<$dataset" or die $!;}
while (<IN>){
	chomp;
	next if $.==1;
	my ($chr,$pos,$ref)=(split /\t/)[0,1,2];
	$chr="chr$chr" unless ($chr =~/^chr/);
	$dataset{"$chr\t$pos"}=$ref;
}
close IN;

my %sites;
open IN,"$samtools view $bam|" or die $!;
while (<IN>){
	chomp;
	my ($FLAG,$CHR,$POS,$CIGAR,$SEQ,$QUAL)=(split /\t/)[1,2,3,5,9,10];
	next unless $FLAG < 256;
	$CHR="chr$CHR" unless ($CHR =~/^chr/);

	my $Cx=$1 if (/\s+Cx:i:(\d+)\s+/);
	my $Cy=$1 if (/\s+Cy:i:(\d+)\s+/);
	my $coordinate="$Cx,$Cy";

	my @SEQ =split //,$SEQ;
	my @QUAL =split //,$QUAL;
	my @CIGAR=$CIGAR=~/(\d+[A-Z])/g;

	my $refPos=$POS;
	my $seqPos =0;
	for(my $i=0;$i<@CIGAR;$i++){
		if($CIGAR[$i]=~/M/){
			$CIGAR[$i]=~s/M//;
			for (my $j=0;$j<$CIGAR[$i];$j++){
				my $currentPos=$refPos+$j;
				my $index="$CHR\t$currentPos";
				if (exists $dataset{$index}){
					my $I="$coordinate\t$index";
					$sites{$I}{'SEQ'}.= $SEQ[$seqPos];
					$sites{$I}{'QUAL'}.= $QUAL[$seqPos];
				}
				$seqPos++;
			}
			$refPos+=$CIGAR[$i];
		}
		elsif ($CIGAR[$i]=~/D/){
			$CIGAR[$i]=~s/D//;
			$refPos+=$CIGAR[$i];
		}
		elsif ($CIGAR[$i]=~/I/){
			$CIGAR[$i]=~s/I//;
			$seqPos+=$CIGAR[$i];
		}
		elsif ($CIGAR[$i]=~/N/){
			$CIGAR[$i]=~s/N//;
			$refPos+=$CIGAR[$i];
		}
		elsif ($CIGAR[$i]=~/S/){
			$CIGAR[$i]=~s/S//;
			$seqPos+=$CIGAR[$i];
		}
		else {
			die "Warning: Unhandled CIGAR $CIGAR[$i] operation.\n";
		}
	}
}
close IN;

my $name=(split /\//,$bam)[-1];
$name=~s/\.$suffix$//;

my $sam2base="$outdir/$name.sam2base.gz";
if($sam2base=~/\.gz$/){open OT,"|gzip >$sam2base" or die $!;}else{open OT,">$sam2base" or die $!;}
foreach my $k(sort keys %sites){
	my ($xy,$chr,$pos)=split /\t/,$k;
	my $index="$chr\t$pos";
	my $seq=$sites{$k}{'SEQ'};
	my $qual=$sites{$k}{'QUAL'};
	print OT "$k\t$dataset{$index}\t$seq\t$qual\n";
}
close OT;

my $output="$outdir/$name.REs.gz";
if($output=~/\.gz$/){open OT,"|gzip >$output" or die $!;}else{open OT,">$output" or die $!;}
if($sam2base=~/\.gz$/){open IN,"gunzip -cd <$sam2base|" or die $!;}else{open IN,"<$sam2base" or die $!;}
while (<IN>){
	chomp;
	my ($xy,$chr,$pos,$refbase,$seq,$qual)=split /\t/;
	my @base=split //,$seq;
	my @qual=split //,$qual;
	my ($cov,$alt)=("0","0");
	for (my $i=0;$i<@base;$i++){
		next if $base[$i] eq "N";
		my $score=ord($qual[$i])-$phred;
		next if $score < $qual_cutoff;
		if ($base[$i] =~/^$refbase$/i){
			$cov++;
		}
		elsif ($base[$i] =~/^$hash{$refbase}$/i){
			$alt++;
			$cov++;
		}
		else {
			next;
		}
	}
	print OT "$xy\t$chr\t$pos\t$refbase\t$cov\t$alt\n" if ($cov >0);
}
close IN;
close OT;

sub usage {
	my ($exit_code) = @_;
    
	print <<"USAGE";

Usage: $0 [options]

Required Options:
	--dataset FILE         Known RNA editing sites dataset file
	--bam FILE             Input BAM file
	--outdir DIR           Output directory
	--samtools PATH        Path to samtools executable (required)

Optional Options:
	--suffix STR           BAM file suffix [default: bam]
	--phred INT            PHRED quality offset [default: 33]
	--qual_cutoff INT      Quality score cutoff [default: 20]
	--help, -h             Display this help message

Description:
	A-to-I RNA editing calling based on Stereo-seq BAM files. Extracts spatial coordinates from Cx:i and Cy:i tags in the BAM file.

Examples:

	$0 --dataset editing_sites.txt --bam sample.bam --outdir results --samtools /usr/bin/samtools

	$0 --dataset editing_sites.gz --bam sample.bam --outdir results --samtools /path/to/samtools --phred 33 --qual_cutoff 20

USAGE

	exit $exit_code;
}

