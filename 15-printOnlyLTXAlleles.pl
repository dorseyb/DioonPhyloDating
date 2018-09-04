use strict;
use warnings;



# A SCRIPT TO FILTER OUT TAXA WITH >X ALLELES FROM ALIGNMENTS/FASTA FILES AND PRINT TO NEW FILES
# SET MAX # ALLELES BELOW
# CAN EXCLUDE OG TAXA AS WELL
# SET INDIR AND OUTDIR BELOW

my $dir = '/Volumes/HD2/FluidigmReads/160602_miseqrun/alignmentsWithOG/';
my $outdir = '/Volumes/HD2/FluidigmReads/160602_miseqrun/fourAlleleAlignmentsIG/';


my $maxAlleles = 4;

opendir(my $dh, $dir) or die $!;
my @files = grep { /afa$/ } readdir $dh;
closedir $dh or die $!;
my $total = scalar(@files);
my $count = 0;

my $Afile = '/Volumes/HD2/FluidigmReads/160602_miseqrun/alleleCounts.txt';
open my $afh, ">", $Afile or die $!;

foreach my $file (@files) {
	print "$file\n";
	my $outfile = $file;
	open my $ofh, ">", "$outdir$outfile" or die $!;
	
	my %allelecounts = ();
	open my $fh, "<", $dir.$file or die $!;
	my $fasref = fasta_to_hash($fh);
	my %fashash = %$fasref;
	
# 244_Ddotstevensonii:401_pre_1-634|cat|1|10|71|0.141	
# /.+\:([\d\w-]+)\|(\w+)\|(\d+)\|.+$/)

	foreach my $key (keys %fashash){
		my $keyN = $key;
		$keyN =~ s/\@.+$//;
		$allelecounts{$keyN}++;
	}
	foreach my $k (sort keys %fashash) {
		my $keyN = $k;
		$keyN =~ s/\@.+$//;
		
		print $afh "$file\t$k\t$allelecounts{$keyN}\n";
		next if ($k =~ /^\d+_[^D].+$/); #comment out to include OG
		my $count = $allelecounts{$keyN};
		next if ($count > $maxAlleles);
		
		print $ofh '>'.$k."\n".$fashash{$k}."\n";
		
	}
	close $ofh or die $!;

}

close $afh or die $!;


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
			$seq = '';
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
