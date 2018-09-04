# oneInd.og12_2555.afa	Pabies.nocontam.gmapdb	exon	168	250	86	-	.	ID=oneInd.og12_2555.afa.mrna1.exon1;Name=oneInd.og12_2555.afa;Parent=oneInd.og12_2555.afa.mrna1;Target=oneInd.og12_2555.afa 168 250 +;Length=83
# oneInd.og12_2555.afa	Pabies.nocontam.gmapdb	exon	251	469	87	-	.	ID=oneInd.og12_2555.afa.mrna1.exon2;Name=oneInd.og12_2555.afa;Parent=oneInd.og12_2555.afa.mrna1;Target=oneInd.og12_2555.afa 251 469 +;Length=219
# oneInd.og12_661.afa	Pabies.nocontam.gmapdb	exon	343	389	86	+	.	ID=oneInd.og12_661.afa.mrna1.exon1;Name=oneInd.og12_661.afa;Parent=oneInd.og12_661.afa.mrna1;Target=oneInd.og12_661.afa 343 389 +;Length=47
# oneInd.og12_661.afa	Pabies.nocontam.gmapdb	exon	390	504	87	+	.	ID=oneInd.og12_661.afa.mrna1.exon2;Name=oneInd.og12_661.afa;Parent=oneInd.og12_661.afa.mrna1;Target=oneInd.og12_661.afa 390 504 +;Length=115
# oneInd.og12_661.afa	Pabies.nocontam.gmapdb	exon	505	549	93	+	.	ID=oneInd.og12_661.afa.mrna1.exon3;Name=oneInd.og12_661.afa;Parent=oneInd.og12_661.afa.mrna1;Target=oneInd.og12_661.afa 505 549 +;Length=45
# oneInd.og12_661.afa	Pabies.nocontam.gmapdb	exon	550	609	87	+	.	ID=oneInd.og12_661.afa.mrna1.exon4;Name=oneInd.og12_661.afa;Parent=oneInd.og12_661.afa.mrna1;Target=oneInd.og12_661.afa 550 609 +;Length=60
# oneInd.og12_661.afa	Pabies.nocontam.gmapdb	exon	610	705	93	+	.	ID=oneInd.og12_661.afa.mrna1.exon5;Name=oneInd.og12_661.afa;Parent=oneInd.og12_661.afa.mrna1;Target=oneInd.og12_661.afa 610 705 +;Length=96
# oneInd.og12_661.afa	Pabies.nocontam.gmapdb	exon	706	791	92	+	.	ID=oneInd.og12_661.afa.mrna1.exon6;Name=oneInd.og12_661.afa;Parent=oneInd.og12_661.afa.mrna1;Target=oneInd.og12_661.afa 706 791 +;Length=86
# oneInd.og12_1097.afa	Pabies.nocontam.gmapdb	exon	3	81	97	+	.	ID=oneInd.og12_1097.afa.mrna1.exon1;Name=oneInd.og12_1097.afa;Parent=oneInd.og12_1097.afa.mrna1;Target=oneInd.og12_1097.afa 3 81 +;Length=79
# oneInd.og12_1097.afa	Pabies.nocontam.gmapdb	exon	82	262	90	+	.	ID=oneInd.og12_1097.afa.mrna1.exon2;Name=oneInd.og12_1097.afa;Parent=oneInd.og12_1097.afa.mrna1;Target=oneInd.og12_1097.afa 82 262 +;Length=181
# oneInd.og12_1097.afa	Pabies.nocontam.gmapdb	exon	263	360	88	+	.	ID=oneInd.og12_1097.afa.mrna1.exon3;Name=oneInd.og12_1097.afa;Parent=oneInd.og12_1097.afa.mrna1;Target=oneInd.og12_1097.afa 263 360 +;Length=98
# oneInd.og12_1097.afa	Pabies.nocontam.gmapdb	exon	361	792	88	+	.	ID=oneInd.og12_1097.afa.mrna1.exon4;Name=oneInd.og12_1097.afa;Parent=oneInd.og12_1097.afa.mrna1;Target=oneInd.og12_1097.afa 361 792 +;Length=432


use strict;
use warnings;

my $dir = '/Users/bdorsey/Documents/Dioon/RNAseq/trinity_files_processed/ortho_search/ortho/oneIndPerSpecies/I12Fasta/gmap/';
my $gff = $dir.'blastHit.abiesMap.targets.gff';
my $outA = $dir.'UpExonsForIntronsGT50.gff';
my $outB = $dir.'DownExonsForIntronsGT50.gff';
my $error = $dir.'AdjacentExonError.txt';

open my $efh, ">$error" or die $!;
open my $gffh, "<$gff" or die $!;
open my $outAfh, ">$outA" or die $!;
open my $outBfh, ">$outB" or die $!;

my @block;
my $prev = '';
my $blockCount = 0;
my $exCount = 0;
my $keepers = 0;
my $pairNum = 0;

while (<$gffh>) {
	next if ($_=~/#/);
	chomp;
	my @fields = split /\t/;
# 	print "$_\n" foreach @fields;
# 	die;

# If at a new locus	
	if ($fields[0] ne $prev) {
# Process old set of exons
		process_block(\@block) if @block;
		$blockCount++;
		@block = \@fields;
	}
	
# If in same locus add exon to block
	else {
		push @block, \@fields;
	}
# Set prev variable to current locus id	
	$prev = $fields[0];
}

# Process the last block
process_block(\@block) if @block;


print "Blocks:$blockCount\nExons:$exCount\nIntrons:$keepers\nPairs:$pairNum\n";


sub process_block {
# set variable equal to argument: ref to array of refs to arrays
	my ($A) = @_;
# 	deref array ref, new array is an array of refs to arrays
	my @exrefs = @{$A};
# 	Declare a hash to hold upstream and downstream exons
	my %hashA = ();
	my %hashB = ();
	
# 	Go through each element in exons array, each is a ref to an array
	if ($#exrefs == 0) {
		my @exon = @{$exrefs[0]};
		my $lineA = join('	', @exon);
		print $efh "$lineA\n";
		$exCount++;
	}
	else {
		$exCount++;
		for my $i (1..$#exrefs) {
			$exCount++;
			my @exon2 = @{$exrefs[$i]};
			my @exon1 = @{$exrefs[$i-1]};

			my $l1 = $exon1[4] - $exon1[3];
			my $l2 = $exon2[4] - $exon2[3];

			my $lineA = join('	', @exon1);#@A
			my $lineB = join('	', @exon2);#@B

# If exons are adjacent
			if ( $exon2[3] == ($exon1[4] + 1)) {

# If both are > 50bp
				if ($l1 >= 50 && $l2 >= 50){
					$pairNum++;
					$lineA =~ s/(ID=oneInd.*?);/$1|$pairNum;/;
					$lineB =~ s/(ID=oneInd.*?);/$1|$pairNum;/;
					$hashA{ $lineA }++;
					$hashB{ $lineB }++;
				}
				elsif ($l1 < 50) {
					print $efh "$lineA\n"
				}
				elsif ($l2 < 50) {	
					print $efh "$lineB\n";
				}	
			}
			else {
				print $efh "$lineA\n";
			}
		}
	}
	foreach my $k (sort byStart keys %hashA) {
		print $outAfh "$k\n";
		$keepers++;
	}
	foreach my $k (sort byStart keys %hashB) {
		print $outBfh "$k\n";
	}

	
}


close $efh or die $!;
close $gffh or die $!;
close $outAfh or die $!;
close $outBfh or die $!;

sub byStart { 
	($a =~ /exon\t(\d+)/)[0] <=> ($b =~ /exon\t(\d+)/)[0] 
}


# 
# 
# hold elements of line in memeory
# holde elements of next line
# compare second id with first id
# if the same id and the start of 2 = end of 1 + 1 then print exon 1 and 2
# save line 2 as new variable
# get next line and compare to line 2 in new variable
# if the same id and the start of 2 = end of 1 + 1 then print exon 1 and 2