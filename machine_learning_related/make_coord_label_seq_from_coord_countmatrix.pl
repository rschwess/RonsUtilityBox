#!/usr/bin/perl -w
use strict;
use warnings;
use Bio::DB::Sam;

# Take a count martix of format chr\tstart\tend\toverlapLabelA\toverlapLabelB etc. with an optional header Line scane for "start" keyword
# take labels from header if present and replace by numbers --> print number to label table
# From the Selected genome(file) extract the Sequence specified
# Collapse the numeric lables and
# Print new file with chr start end label seq in traindata.txt

# usage: perl make_coord_label_seq_from_coord_countmatrix.pl <genome> <input_countmatrix.txt> <output_traindata> <output_label_number>

# Match and Get reference genome
my $genome = $ARGV[0];
my $genome_file;
if($genome eq "hg18"){
	$genome_file = '/databank/igenomes/Homo_sapiens/UCSC/hg18/Sequence/WholeGenomeFasta/genome.fa';
}elsif($genome eq "hg19"){
	$genome_file = '/databank/igenomes/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa';
}elsif($genome eq "mm9"){
	$genome_file = '/databank/igenomes/Mus_musculus/UCSC/mm9/Sequence/WholeGenomeFasta/genome.fa';
}else{
	print "Please select an implemented reference genome, hg19, hg18, mm9 !\n";
	exit 2;
}

# Access Reference Genome
my $fai = Bio::DB::Sam::Fai->load("$genome_file");
my $inputfile=$ARGV[1];
my $output_traindata=$ARGV[2];
my $output_label_number=$ARGV[3];

my ($chromosome, $start, $stop, @rest, @labels, $num_labels);

open(IN, $inputfile) or die "Can't open input file $inputfile ! $!";

  # get and check firstline
  my $firstline = <IN>;

  # if header line present --> get the label names from here
  if($firstline =~ /start/){
    ($chromosome, $start, $stop, @labels) = split(/\s+/, $firstline);
    $num_labels = @labels;
    $num_labels--;
    # open table output
    open(OUT1, ">$output_label_number") or die "Can't open Output File for writing Labeling $output_label_number $! !\n";
    for(my $i = 0; $i <= $num_labels; $i++){
      print OUT1 "$i\t$labels[$i]\n";
    }
    close(OUT1);
  }

  # open Output traindata file
  open(OUT2, ">$output_traindata") or die "Can't open Output File for writing Training Data $output_traindata $! !\n";

  while(<IN>){

    chomp;
    ($chromosome, $start, $stop, @rest) = split(/\s+/, $_);

    # correct for bed coordinates
    $start++;

    my $fai_location = $chromosome.':'.$start.'-'.$stop;
    my $refseq = $fai->fetch($fai_location); # fetch sequence
    $refseq=~s/a/A/g;  # clean up Sequence
    $refseq=~s/g/G/g;
    $refseq=~s/t/T/g;
    $refseq=~s/c/C/g;

    # reset to bed coordinates
    $start--;

    # go through labels overlaps and fetch active lables (1 or more) and concatinate them to a comma joined string
    my $conc_labels = "";
    my $j=0;
    foreach (@rest){
      $conc_labels .= ",$j" if $_ >= 1;
      $j++;
    }
    # remove first comma (cleanup string)
    $conc_labels =~s/^\,//g;
    # if no label is active --> skip sequence
    next if $conc_labels eq "";

		print OUT2 "$chromosome\t$start\t$stop\t$conc_labels\t$refseq\n";

    }

close(IN);
close(OUT2);
