# Script to extract 50 bp from facing ends of exons to use for primer design targeting introns
# Updating to concatenate extracted sequences in order to choose primer pairs around splice site
# rather than choosing single primers in up or down stream sequence individually

my $dir = '/Users/bdorsey/Documents/Dioon/RNAseq/trinity_files_processed/ortho_search/ortho/oneIndPerSpecies/I12Fasta/exons/';
my $fasta = $dir.'targetExons_noPos.fa';
my $upout = $dir.'upExons_secondRun.p3';
my $downout = $dir.'downExons_secondRun.p3';
my $upids = $dir.'upExonIDs.txt';
my $downids = $dir.'downExonIDs.txt';
my $catseq = $dir.'spliceSiteSeqs.p3';

open my $catfh, ">$catseq" or die $!;
open my $fh, "<$fasta" or die $!;
# open my $ufh, ">$upout" or die $!;
# open my $dfh, ">$downout" or die $!;
open my $uidfh, "<$upids" or die $!;
open my $didfh, "<$downids" or die $!;

my $fasref = fasta_to_hash($fh);
my %fashash = %{ $fasref };

# Go through up exon ids
# Foreach id capture id, intron number
# Get sequence from fasta hash
# Do same for down ids
# Concatenate up and down seqs
# Print to p3 file using intron number as id

my %cathash = ();

while (my $uid = <$uidfh>) {
	chomp $uid;
	my ($def,$Inum) = split '\|', $uid;
# 	print "$def\n";
	my $useq = substr($fashash{$def}, -50);
	$cathash{$Inum} = [ $useq ];
	print "$Inum\t$useq\n";
}

while (my $did = <$didfh>) {
	chomp $did;
	my ($def,$Inum) = split '\|', $did;
	my $dseq = substr($fashash{$def}, 0,50);
	push @{ $cathash{$Inum} }, $dseq;
	print "$Inum\t$dseq\n";
}

foreach my $I (keys %cathash) {
	my @a = @{ $cathash{$I} };
# 	my $seq = join('', @a);
	my $seq = $a[0].$a[1];
	print $catfh 'SEQUENCE_ID=Intron_'."$I\n".'SEQUENCE_TEMPLATE='."$seq\n=\n";
}



# while (my $def = <$uidfh>) {
# 	chomp $def;
# 	my $s = $def;
# 	my $num = $s =~ s/\|\d+//;
# # 	my @k = grep /$s/, (keys %fashash);
# # 	my $seq = $fashash{$k[0]};
# 	my $seq = $fashash{$s};
# 	my $pseq = substr($seq, -50);
# # 	print $ufh '>'."$def\n$pseq\n";
# 	print $ufh 'SEQUENCE_ID='."$def\n".'SEQUENCE_TEMPLATE='."$pseq\n=\n";
# }
# while (my $def = <$didfh>) {
# 	chomp $def;
# 	my $s = $def;
# 	$s =~ s/\|\d+//;
# # 	print "$s\n";
# 	my @k = grep /$s/, (keys %fashash);
# 	my $seq = $fashash{$k[0]};
# 	my $pseq = substr($seq, 0,50);
# # 	print $dfh '>'."$def\n$pseq\n";
# 	print $dfh 'SEQUENCE_ID='."$def\n".'SEQUENCE_TEMPLATE='."$pseq\n=\n";
# }

close $fh or die;
# close $ufh or die;
# close $dfh or die;
close $uidfh or die;
close $didfh or die;
close $catfh or die;


sub fasta_to_hash {
my ($fastafh) = @_; 
my $header;
my $seq;
my %sequence;

while (<$fastafh>) {
	chomp($_);
 
	# If the line starts with ">" its a header, and we save the information in $header
	if ($_=~/^>(.*)/) {
 
		# if we have something in $seq
		if ($seq) {
			$sequence{$header}=$seq;
		}
 
 
		$header = $1;
		# reset sequence!!!
		$seq    = '';
	} else {
		$seq.=$_;
	}
 
}
# Store last sequence
$sequence{$header}=$seq;

return \%sequence;
}
