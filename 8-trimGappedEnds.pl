# A script to trim gap postitions from beginning and end of an alignment
# and calculate a consensus

use strict;
use warnings;
$|=1;

my $dir = '/Users/bdorsey/Documents/Dioon/RNAseq/trinity_files_processed/ortho_search/ortho/oneIndPerSpecies/I12Fasta/';
# my @targets = qw/oneInd.og12_1007 oneInd.og12_1038 oneInd.og12_1039 oneInd.og12_1041 oneInd.og12_1044 oneInd.og12_1046 oneInd.og12_1052 oneInd.og12_1055 oneInd.og12_1067 oneInd.og12_1073 oneInd.og12_1075 oneInd.og12_1082 oneInd.og12_1083 oneInd.og12_1084 oneInd.og12_1090 oneInd.og12_1095 oneInd.og12_1110 oneInd.og12_1112 oneInd.og12_1115 oneInd.og12_1123 oneInd.og12_1135 oneInd.og12_1147 oneInd.og12_1148 oneInd.og12_1150 oneInd.og12_1162 oneInd.og12_1164 oneInd.og12_1172 oneInd.og12_1193 oneInd.og12_1212 oneInd.og12_1214 oneInd.og12_1218 oneInd.og12_1257 oneInd.og12_1272 oneInd.og12_1289 oneInd.og12_1294 oneInd.og12_1295 oneInd.og12_1301 oneInd.og12_1313 oneInd.og12_1334 oneInd.og12_1356 oneInd.og12_1373 oneInd.og12_1378 oneInd.og12_1448 oneInd.og12_1578 oneInd.og12_1629 oneInd.og12_1697 oneInd.og12_1855 oneInd.og12_1905 oneInd.og12_1908 oneInd.og12_1920 oneInd.og12_1938 oneInd.og12_1956 oneInd.og12_1998 oneInd.og12_2018 oneInd.og12_2089 oneInd.og12_2095 oneInd.og12_2103 oneInd.og12_2120 oneInd.og12_2135 oneInd.og12_2139 oneInd.og12_2154 oneInd.og12_2166 oneInd.og12_2174 oneInd.og12_2227 oneInd.og12_2235 oneInd.og12_2393 oneInd.og12_2442 oneInd.og12_2553 oneInd.og12_2844 oneInd.og12_2871 oneInd.og12_2887 oneInd.og12_2960 oneInd.og12_3046 oneInd.og12_3084 oneInd.og12_3089 oneInd.og12_3097 oneInd.og12_3158 oneInd.og12_3159 oneInd.og12_3236 oneInd.og12_3297 oneInd.og12_3333 oneInd.og12_3338 oneInd.og12_3354 oneInd.og12_3356 oneInd.og12_3366 oneInd.og12_3373 oneInd.og12_3411 oneInd.og12_3412 oneInd.og12_3444 oneInd.og12_3451 oneInd.og12_404 oneInd.og12_408 oneInd.og12_426 oneInd.og12_441 oneInd.og12_444 oneInd.og12_458 oneInd.og12_461 oneInd.og12_475 oneInd.og12_478 oneInd.og12_491 oneInd.og12_509 oneInd.og12_526 oneInd.og12_532 oneInd.og12_547 oneInd.og12_550 oneInd.og12_551 oneInd.og12_554 oneInd.og12_570 oneInd.og12_585 oneInd.og12_594 oneInd.og12_614 oneInd.og12_617 oneInd.og12_639 oneInd.og12_642 oneInd.og12_655 oneInd.og12_665 oneInd.og12_677 oneInd.og12_678 oneInd.og12_701 oneInd.og12_716 oneInd.og12_731 oneInd.og12_749 oneInd.og12_751 oneInd.og12_772 oneInd.og12_775 oneInd.og12_800 oneInd.og12_804 oneInd.og12_808 oneInd.og12_809 oneInd.og12_814 oneInd.og12_829 oneInd.og12_839 oneInd.og12_847 oneInd.og12_858 oneInd.og12_860 oneInd.og12_862 oneInd.og12_873 oneInd.og12_876 oneInd.og12_888 oneInd.og12_889 oneInd.og12_906 oneInd.og12_919 oneInd.og12_922 oneInd.og12_927 oneInd.og12_939 oneInd.og12_942 oneInd.og12_947 oneInd.og12_948 oneInd.og12_960 oneInd.og12_964 oneInd.og12_966 oneInd.og12_970 oneInd.og12_971 oneInd.og12_981 oneInd.og12_983 oneInd.og12_993 oneInd.og12_995/;
# my @targets = ('oneInd.og12_1108');

my $tfile = $dir.'blasthit.conIDs.txt';
my $confile = $dir.'/nuc/nuc.cons.all.fa';
my $newcon = $dir.'/nuc/blastHit.conseqs.thresholdTrim.fa';
my $lengthfile = $dir.'blastHit.thresholdTrimLengths.txt';

open my $tfh, "<$tfile" or die $!;
my @targets = <$tfh>;
chomp @targets;
close $tfh or die $!;

open my $lfh, ">$lengthfile" or die $!;
print $lfh "ID\talignLength\ttrimLength\n";

my $count = 0;

open my $newcfh, ">$newcon" or die $!;

open my $cfh, "<$confile" or die $!;
my $conref = &fasta_to_hash($cfh);
my %conhash = %{ $conref };
close $cfh or die $!;

# Open each file
foreach my $target (@targets) {
	
	my $file = $dir.'nuc/'.$target.'.afa';
	open my $fh, "<$file" or die $!;
	
	my $fasref = &fasta_to_hash($fh);
	my %fashash = %{ $fasref };
	
# 	my $out = $dir.'/nuc/'.$target.'_trim.afa';
# 	open my $outfh, ">$out" or die $!;
	
	my $longFront = 0;
	my $longEnd = 0;


# get pad lengths, store in hash, rank, choose length of sequence that adds 1/2 number of total seqs
	my %frontHash = ();
	my %endHash = ();
	my $numSeqs = keys %fashash;
	my $threshold = $numSeqs/2;
	
	foreach my $key (keys %fashash) {
		if ($fashash{$key} =~ /^(-+)/){
			$frontHash{$key} = length($1);
		}
		if 	($fashash{$key} =~ /(-+)$/){
			$endHash{$key} = length($1);
		}
	}

########################################
# THIS GOES ONE VALUE OF V TOO FAR!!!!

	my $numF = 1;
	my $numE = 1;
	foreach my $v ( sort { $a <=> $b } values %frontHash) {
# 		print "$numF\n$v\n$longFront\n";
		last if $numF >= $threshold;
		$longFront = $v;
		$numF++;
	}
	foreach my $v ( sort { $a <=> $b } values %endHash) {
		last if $numE >= $threshold;
		$longEnd = $v;
# 		$longEnd = $v unless $numE >= $threshold;
		$numE++;
	}

# Get length of padded ends, choose longest
# 	foreach my $key (keys %fashash) {
# 		if ($fashash{$key} =~ /^(-+)/){
# 			my $frontGap = length($1);
# 			$longFront = $frontGap if ($frontGap > $longFront);
# 		}
# 
# 		if 	($fashash{$key} =~ /(-+)$/){
# 			my $endGap = length($1);
# 			$longEnd = $endGap if ($endGap > $longEnd)
# 		}
# 
# 	}
# 	print $longFront."\n".$longEnd."\n";

	my $off = $longFront;
	my $l = 0-$longEnd;

# Extract seq between padded ends from consensus seq and print
	my $cid = $target.'.afa';
	my $conseq = $conhash{$cid};
	my $alength = length($conseq);
	my $cseq = '';
	if ($l < 0){
		$cseq = uc(substr($conseq, $off, $l));
	}
	elsif ($l == 0){
		$cseq = uc(substr($conseq, $off));
	}
	else { 
		die "$cid end gap length greater than 0\n";
	}
# 	print $off."\n".$l."\n";
	my $tlength = length($cseq);
	
	print $newcfh '>'.$cid."\n".$cseq."\n";
	print $lfh "$cid\t$alength\t$tlength\n";
	
	$count++;
	print "\r$count";

#######################################################
# Not writing new alignments, just conseqs above.
#  Extract seqs between padded ends and print
# 	my $longest = 0;
# 	foreach my $key2 (keys %fashash) {
# 		my $seq = '';
# 		if ($l < 0){
# 			$seq = substr($fashash{$key2}, $off, $l);
# 		}
# 		elsif ($l == 0){
# 			$seq = substr($fashash{$key2}, $off);
# 		}
# 		else { 
# 			die "end gap length greater than 0\n";
# 		}
# 		print $off."\n".$l."\n";
# 		print '>'.$key2."\n".$seq."\n";
# 	}
#######################################################

	close $fh or die $!;
# 	close $outfh or die $!;
}
print "\n";
close $newcfh or die $!;
close $lfh or die $!;



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
