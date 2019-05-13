#!/usr/bin/perl -w
use strict;
use warnings;
use Bio::DB::Sam;


#input reference genome
my $genome = $ARGV[0];

my $inputfile=$ARGV[1];

my $onebased = $ARGV[2]; # 0 for ) base dor 1 for 1 base dindexecs

my $genome_file;

if($genome eq "hg18"){
	$genome_file = '/databank/igenomes/Homo_sapiens/UCSC/hg18/Sequence/WholeGenomeFasta/genome.fa';
}elsif($genome eq "hg19"){
	$genome_file = '/databank/igenomes/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa';
}elsif($genome eq "mm9"){
	$genome_file = '/databank/igenomes/Mus_musculus/UCSC/mm9/Sequence/WholeGenomeFasta/genome.fa';
}else{
	$genome_file = $genome;
}

if($onebased != 0 && $onebased != 1){
	print "Pass 0 or 1 for specifying 0 or 1 based encoding";
	exit 2;
}


# load genome sequence
my $fai = Bio::DB::Sam::Fai->load("$genome_file");

open(IN, $inputfile) or die "Can't open input file $inputfile ! $!";

	while(<IN>){

		chomp;

		my ($chromosome, $start, $stop, @rest) = split(/\s+/, $_);

		# correct for bed coordinates
		if($onebased == 0){
			$start++;
		}
	
		# get sequence 
		my $fai_location = $chromosome.':'.$start.'-'.$stop;
		my $refseq = $fai->fetch($fai_location);

		$refseq=~s/a/A/g;
		$refseq=~s/g/G/g;
		$refseq=~s/t/T/g;
		$refseq=~s/c/C/g;

		# split sequence into bases
		my @seqsplit = split("", $refseq);

		# init base pair coordinates
		my $bstart = $start;
		my $bend = $start + 1;
		
		for(my $i = 0; $i < scalar(@seqsplit); $i++){

			# get base 
			my @bases = ("A", "C", "G", "T");
			my $current_base = $seqsplit[$i];
			# make variations and print
			foreach my $printbase (@bases){
				print "$chromosome\t$bstart\t$bend\t$printbase\n" if $printbase ne $current_base;
			}

			# count up positions
			$bstart++;
			$bend++;
		}
	

	}

close(IN);
