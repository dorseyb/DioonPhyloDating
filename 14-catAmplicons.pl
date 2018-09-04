use strict;
use warnings;


my $dir = '/Volumes/HD2/FluidigmReads/160602_miseqrun/Dioon-phylo.reduced/occurrence.split_amplicon/';
# my $outdir = '/Volumes/HD2/FluidigmReads/160602_miseqrun/alignmentsNoOG/'; # change if including outgroups
my $outdir = '/Volumes/HD2/FluidigmReads/160602_miseqrun/fastaOG/';

my $stoutfile = $outdir.'stdoutFile.txt';

opendir (my $dh, $dir) or die $!;
my @files = grep { /^Amplicon/ } readdir $dh;
closedir $dh or die $!;

foreach my $file (@files) {
# 	PRINT MERGED FILES TO NEW DIRECTORY WITH OR WITHOUT OUTGROUPS DEPENDING
	if ($file =~ /merged/) {
		my $r1 = $file;
		$r1 =~ s/merged/read1/;
		next if (-s "$dir$r1"); # SKIP MERGED FILES IF F AND R FILES ALSO EXIST

		open my $mfh, "<", "$dir$file" or die $!;
		my $mref = fasta_to_hash($mfh);
		close $mfh or die $!;
		my %mhash = %$mref;		
		my $mout = "$outdir$file";
		open my $moutfh, ">", $mout or die $1;

		foreach my $mkey (keys %mhash) {
# 			next if ($mkey !~ /^\d+_D.+/); # comment out to include outgroups
			print "$file\n" if ($mkey !~ /^\d+_D.+/);
			print $moutfh '>'.$mkey."\n".$mhash{$mkey}."\n";
		}
		
		close $moutfh or die $!;
	}
# 	JOIN F AND R READS, ADD MERGED READS, AND PRINT TO FILE IN NEW DIRECTORY
	if ($file =~ /read1/) {
		my $file2 = $file;
		my $fileC = $file;
		$file2 =~ s/read1/read2/;
		$fileC =~ s/read1/cat/;
		open my $fh1, "<", "$dir$file" or die $!;
		open my $fh2, "<", "$dir$file2" or die $!;
		open my $cfh, ">", "$outdir$fileC" or die $!;
		
		my $r1ref = fasta_to_hash($fh1);
		my %r1hash = %$r1ref;
		
		my $r2ref = fasta_to_hash($fh2);
		my %r2hash = %$r2ref;
		
		foreach my $key1 (keys %r1hash) {
# 			next if ($key1 !~ /^\d+_D.+/); # comment out to include outgroups
			print "$file\n" if ($key1 !~ /^\d+_D.+/);
			my $key2 = $key1;
			$key2 =~ s/read1/read2/;
			my $catread = "$r1hash{$key1}$r2hash{$key2}";
			my $catdef = $key1;
			$catdef =~ s/read1/cat/;
			
			print $cfh '>'.$catdef."\n".$catread."\n";
		}
		my $mfile = $file;
		$mfile =~ s/read1/merged/;
		if (-s "$dir$mfile") {
			open my $mfh, "<", "$dir$mfile" or die $!;
			my $mref = fasta_to_hash($mfh);
			close $mfh or die $!;
			my %mhash = %$mref;
			foreach my $mkey (keys %mhash){
# 				next if ($mkey !~ /^\d+_D.+/); # comment out to include outgroups
				print "$mfile\n" if ($mkey !~ /^\d+_D.+/);
				print $cfh '>'.$mkey."\n".$mhash{$mkey}."\n";
			}
			
		}
		
		close $fh1 or die $!;
		close $fh2 or die $!;
		close $cfh or die $!;
	
	}
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
