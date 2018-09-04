#!/usr/bin/env perl

use strict;
use warnings;
$|=1;

############################################################
# Script to count multiple orf's inferred by Transdecoder
############################################################



my $seqdir = '/Users/bdorsey/Documents/Dioon/RNAseq/trinity_files_processed/';
# my $libdir = '/Volumes/HD2/RNAseq_cleaning/lib_dirs2.txt';
my $out = $seqdir .'multiple.orf.seqs_2.txt';
open (my $outfh, ">$out") or die;

my @libdir = qw/sonorense_both
purpusii2_RYBD010_index16
purpusii1_RYBD008_index14
merolae_RYBD007_index13
mejiae1_RYBD001_index2
mejiae2_both
edule2_RYBD011_index18
spinulosum1_RYBD004_index6
spinulosum2_RYBD012_index19
edule1_RYBD009_index15/;


my $count = 0;
my $filecount = 0;

INDEX:foreach my $index (@libdir) {
	my %seqhash;
	my $file = $seqdir . $index .'/transdecoder/best_candidates.eclipsed_orfs_removed.pep';
	open(my $filefh, "<$file") or die $!;
	
# 	>m.7 g.7  ORF g.7 m.7 type:5prime_partial len:160 (+) edule1_comp18156_c0_seq1_len=936:3-482(+)
# 	>m.5229 g.5229  ORF g.5229 m.5229 type:complete len:370 (-) edule2_comp48200_c0_seq1_len=2429:1319-2428(-)
	LINE: while (my $line = <$filefh>){
		chomp $line;
		next LINE unless ($line =~ /^>/  );
		if ($line =~ /^>/  ) {
			my @defline = (split /\s+/, $line);
			my $tgene = $defline[1];
			my $def = $defline[8];
			my $len = ':'.$defline[6];
			$def =~ s/:(.*)$//;
# 			print "$def\n";
# 			$seqhash{$def} = [] unless(exists $seqhash{$def});
# 			push @{$seqhash{$def}}, ($1, $len);
			$seqhash{$def} = {} unless (exists $seqhash{$def});
			$seqhash{$def}{$tgene} = $1;
		}
	}
	close $filefh or die $!;
# 	print "\r$count/$filecount";
	
# 	foreach (keys %seqhash){
# 		if (scalar(@{ $seqhash{$_} }) > 2 ) {
# 			print $outfh "$_\t";
# 			foreach my $pos (@{ $seqhash{$_} }) {
# 				
# 				print $outfh "$pos\t";
# 			}
# 			print $outfh "\n";
# 		}	
# 	}


	foreach my $trans (keys %seqhash)
		{
		my %hash = %{ $seqhash{$trans} };
		if ( scalar( keys %hash )  > 1 )
			{
			foreach my $orf ( keys %hash )
				{
				print $outfh "$trans\t$orf\t$hash{$orf}\n";
				}
			}
		}
	
}
close $outfh or die;
