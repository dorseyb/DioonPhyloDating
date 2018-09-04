use strict;
use warnings;

my ($indir, $outdir, $ending, $tag) = @ARGV;

# die "in_directory, out_directory, file_ending, [tag]" unless (scalar @ARGV == 5);

$ending =~ s/^\.//;

opendir my $indh, $indir or die $!;
my @files = grep { /.*$ending$/ } readdir $indh;
$outdir =~ s|/?$|/|;
$indir =~ s|/?$|/|;

foreach my $file (@files) {
	my $outfile = $file;
	$outfile =~ s/afa/$tag.afa/;
	open my $ofh, ">", $outdir.$outfile or die $!;
	open my $fh, "<", $indir.$file or die $!;
	my $ref = fasta_to_hash($fh);
	my %hash = %$ref;
	my %newhash = ();
	my $check = 0;

# CHECK FOR ANY TAXON WITH MULTIPLE ALLELES
	foreach my $key (keys %hash) {
		my $query = $key;
		$query =~ s/\@.+$//;
		my $count = grep /$query/, keys %hash;
		$check++ if ($count > 1) 
	}
# IF MULTIPLE ALLELES PRESENT, DUPLICATE HOMOZYGOTES
	if ($check) {
		foreach my $key (keys %hash) {
			my $query = $key;
			$query =~ s/\@.+$//;
			my $count = grep /$query/, keys %hash;
		
			if ($count < 2) {
				my $hom = $key;
				$hom =~ s/_\d+$/_h/;
				$newhash{$hom} = $hash{$key};
				$newhash{$key} = $hash{$key};
			}
			else {
				$newhash{$key} = $hash{$key};
			}
		}
	}
	foreach my $new (keys %newhash) {
		print $ofh '>'.$new."\n".$newhash{$new}."\n";
	}
	
	
	close $fh or die $!;
	close $ofh or die $!;
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
		if ($_=~/^>(.+)/)
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



