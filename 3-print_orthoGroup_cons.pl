#!/Users/bdorsey/perl5/perlbrew/perls/perl-5.22.0/bin/perl


# SCRIPT TO PRINT CONSENSUS SEQUENCES OF ORTHOMCL ORTHOGROUPS

use strict;
use warnings;
use Bio::AlignIO;
use Bio::SimpleAlign; #qw(each_seq_with_id );


my $dir = '/Users/bdorsey/Documents/Dioon/RNAseq/trinity_files_processed/ortho_search/ortho/oneIndPerSpecies/I12Fasta/pep/';

opendir(my $dh, $dir) or die;
my @aln = grep {/\.afa/} readdir($dh);
closedir $dh or die $!;

my $out = $dir.'conseqs.fa';
open my $outfh, ">>$out" or die $!;

foreach my $f (@aln) {
# Create alignment object:
	my $file = $dir.$f;
	my $in = Bio::AlignIO->new(-file => $file);
	my $aln = $in->next_aln;

# Calculate consensus and print to common file	
	my $con = $aln->consensus_string();
	print $outfh '>'.$f."\n".$con."\n";

}

close $outfh or die $!;