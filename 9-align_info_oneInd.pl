

####################################################################
# Script to calculate several metrics for each orthogroup alignment
####################################################################

use strict;
use warnings;
use Bio::AlignIO;
use Bio::SimpleAlign; #qw(each_seq_with_id );
use List::MoreUtils;
use Scalar::Util qw(looks_like_number);
use Math::Round;

$| = 1;
my $dir = '/Users/bdorsey/Documents/Dioon/RNAseq/trinity_files_processed/ortho_search/ortho/oneIndPerSpecies/I12Fasta/nuc/';

# For all alignments
# opendir(my $dh, $dir) or die;
# my @aln = grep {/\.afa/} readdir($dh);
# closedir $dh or die $!;

# For only blast hits seqs mapping uniquely to P abies genome
my $idfile = $dir.'blastHit.pabiesUniq.IDs.txt';
open my $ifh, "<$idfile" or die $!;
my @aln = <$ifh>;
chomp @aln;
close $ifh or die $!;

my $numAlign = @aln;
print $numAlign."\n";

my $Iout = $dir.'blastHit.abiesMap.AlignInfo.txt';
open my $Ioutfh, ">>$Iout" or die $!;

my $Cout = $dir.'nuc.cons.all.fa';
# open my $Coutfh, ">>$Cout" or die $!;

print $Ioutfh join( "\t", qw|alignment num_seqs length no_gap_length mismatches pid opid percent_gap pis percent_pis nogap_pid nogap_opid nogap_pis nogap_percent_pis| ) ."\n";

my $counter = 0;
foreach my $f (@aln) {

# Create alignment object, check for gt 3 taxa:
	my $file = $dir.$f;
	my $in = Bio::AlignIO->new(-file => $file, -format => 'fasta');
	my $aln = $in->next_aln;
	my $taxa = $aln->num_sequences();
# 	next if $taxa < 4;
	
# Calculate consensus and print to common file	
# 	my $con = $aln->consensus_iupac();
# 	print $Coutfh '>'.$f."\n".$con."\n";
	
# Get info for alignment and print to file - prepended by '#alignment' line to delimit individual ortho groups
	&run_info($aln, $f);
		
	++$counter;
	print "\r$counter";

}

close $Ioutfh or die $!;
# close $Coutfh or die $!;
print "\n";

sub run_info
	{
	my ($aln, $alignment) = @_;

	my ($mismatches, $len, $len2, $pid, $opid, $percent_gap, $depth, $pis, $percent_pis, $nogap_pid, $nogap_opid, $nogap_pis, $nogap_percent_pis); #%pidhash, %opidhash, %gaphash, 
	
	if ( $aln -> is_flush() )
		{
		$len = $aln -> length();
		my $nogap = $aln -> remove_columns(['gaps']);
		$len2 = $nogap -> length();
		$depth = $aln -> num_sequences();
		if ( $len2 != 0)
			{
			my $matches = $nogap->match_line();
			$mismatches = $matches =~ s/\s/!/g;
			$mismatches = '0' unless ($mismatches);
			($nogap_pid, $nogap_pis) = &inform($nogap);
			$nogap_percent_pis = ($nogap_pis/$len2)*100;
			$nogap_opid = $nogap->overall_percentage_identity();
			$percent_gap = (1-($len2/$len))*100;
			}
		else
			{
			$nogap_pid = 'all gaps';
			$nogap_opid = 'all gaps';
			$percent_gap = 1;
			$nogap_pis = 'all gaps';
			$nogap_percent_pis = 'all gaps';
			}
		($pid, $pis) = &inform($aln);
		$percent_pis = $pis/$len;
		$opid = $aln->overall_percentage_identity();
		
# 		$pidhash{$alignment} = $pid;
# 		$opidhash{$alignment} = $opid;
# 		$gaphash{$alignment} = $percent_gap;
#
			my @list = ( $alignment, (nearest(0.001,( $depth, $len, $len2, $mismatches, $pid, $opid, $percent_gap, $pis, $percent_pis, $nogap_pid, $nogap_opid, $nogap_pis, $nogap_percent_pis) ) ) );
			print $Ioutfh join("\t", @list ) . "\n";
		}
	else
		{
		print "$alignment not flush\n";
		next;
		}
	}

sub inform { # Most of this code is taken from Bio::SimpleAlign::average_percentage_identity
	
	my ($self,@args) = @_;

	my @alphabet = ('A','B','C','D','E','F','G','H','I','J','K','L','M',
	'N','O','P','Q','R','S','T','U','V','W','X','Y','Z');

	my ($len, $total, $subtotal, $divisor, $subdivisor, @seqs, @countHashes, $pis);
	

	if (! $self->is_flush()){
		$self->throw("All sequences in the alignment must be the same length");
	}

	@seqs = $self->each_seq();
	$len = $self->length();

	# load the each hash with correct keys for existence checks

	for( my $index=0; $index < $len; $index++){
		foreach my $letter (@alphabet){
			$countHashes[$index]->{$letter} = 0;
		}
	}
	foreach my $seq (@seqs)
		{
		my @seqChars = split //, $seq->seq();
		for( my $column=0; $column < @seqChars; $column++ ){
			my $char = uc($seqChars[$column]);
			if (exists $countHashes[$column]->{$char}){
				$countHashes[$column]->{$char}++; #count number of each character in each column
			}
		}
	}

	$total = 0;
	$divisor = 0;
	$pis = 0;
	for(my $column =0; $column < $len; $column++){
		my %hash = %{$countHashes[$column]}; # a->1, c->0, t->5, etc
		my $vals = 0;
		my $present = 0;
		
		foreach (values %hash) {
			$present++ if ($_); # count bases present at position, i.e. if value for key (base) is gt 0
			$vals++ if ($_ >= 2); # count bases that are found 2 or more times
		}
		
		$pis++ if ($vals >= 2 ); # || ($vals == 1 && ($taxa == 3 || $present >= 2)PI if there are 2 or more multiples, or if there is one multiple and more than one base at a site
		
		$subdivisor = 0;
		foreach my $res (keys %hash){
			$total += $hash{$res}*($hash{$res} - 1); # count of a base present times count - 1; e.g. if a->3, add 3*2 to total; 
			$subdivisor += $hash{$res}; # sum num of bases
		}
		$divisor += $subdivisor * ($subdivisor - 1); # add 


	}
	return ($divisor > 0 ? ($total / $divisor )*100.0 : 0), $pis;
}


sub fasta_to_hash
	{
	my ($fastafh) = @_; 
	my $header;
	my $seq;
	my %sequence;

	while (<$fastafh>)
		{
		chomp($_);
 
		# If the line stars for ">" is a header, and we save the information in $header
		if ($_=~/^>(.*)/)
			{
 			# if we have something in $seq
			if ($seq)
				{
				$sequence{$header}=$seq;
				}
 
 
			$header = $1;
	# 		$counter++;
			# reset sequence!!!
			$seq    = '';
			}
		else
			{
			$seq.=$_;
			}
 
		}
	# Store last sequence
	$sequence{$header}=$seq;

	return \%sequence;
	}
