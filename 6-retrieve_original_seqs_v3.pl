#!/usr/bin/env perl


# oneInd.og12_5993: pur1|m616g616 spi2|m1873g1873

use strict;
use warnings;
$|=1;

########################################################################
# Script to retrieve original nucleotide sequences 
# corresponding to chosen protein sequences after
# translation, ortholog search, paralog search, multi-orf filter

# This is version 3 to be used with one individual
# per species orthogroups AND retrieving the entire
# nuc sequence - not using position info.
#########################################################################

# my $pepdir = '/Users/bdorsey/Documents/Dioon/RNAseq/trinity_files_processed/ortho_search/ortho/bySpeciesResults/ogfasta/I20/all/rax_pruner_files/';

##### For retrieving non_paralog seqs
# my $pepdir = '/Users/bdorsey/Documents/Dioon/RNAseq/trinity_files_processed/ortho_search/ortho/bySpeciesResults/ogfasta/I20/all/';
my $dir = '/Users/bdorsey/Documents/Dioon/RNAseq/trinity_files_processed/ortho_search/ortho/oneIndPerSpecies/';
my $outdir = '/Users/bdorsey/Documents/Dioon/Microsats/';

my $oglistfile = $dir.'I12Fasta/no_paralog_groups.txt';
my $ogfile = $dir.'one.groups.I12.txt';

my $indexdir = '/Users/bdorsey/Documents/Dioon/RNAseq/trinity_files_processed/';
my $testfile = '/Users/bdorsey/Documents/Dioon/RNAseq/trinity_files_processed/ortho_search/ortho/bySpeciesResults/ogfasta/I20/all/rax_pruner_files/retrieve.checkfile.txt';
my (%spechash, %cutspec, %cuthash, %masterhash, %position_hash,	%orf_info_hash);


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

my @abbrev = qw/sono pur2 pur1 mero mej1 mej2 edu2 spi1 spi2 edu1/;
@spechash{@abbrev} = @libdir;

# Define hash to translate original library indices - in @libdir - to names used in fasta from chimera cutting
# cutspec{index} = abbreviated sp. name in defline
# e.g. cutspec{edule2_RYBD011_index18} = edule2

foreach (@libdir){
	if ($_ eq 'mejiae2_both' || $_ eq 'sonorense_both') {
		$cutspec{$_} = $_;
	}
	else{
		my $in = $_;
		$in =~ s/_.+//;
		$cutspec{$_} = $in;
	}
}

HASHES: # For each index, read through corresponding transdecoder fasta and parse deflines
	    # Create hashes that connect transgene ids and trinity ids and one that holds position and strand info for each tringene
	    # Then read the chimera cut fasta into a hash and assign that as a value to the key 'index' in %cuthash

foreach my $index (@libdir) {

# define species from hash created above if necessary
	my (@cds_id_list, $cutsp);
# 	if ($index eq 'sonorense_both') {
# 		$cutsp = $index;
# 	}
# 	elsif ($index eq 'mejiae2') {
# 		$cutsp = $index;
# 	}
# 	else{
		$cutsp = $cutspec{$index};
# 	}

# define the path to the orf seq file from transdecoder and open the file, then readline
	my $cds = $indexdir . $index . '/transdecoder/' . $index . '.orfs.cds';

	open(my $cdsfh, "<$cds") || die "can't open $cds $! ";
	while (<$cdsfh>){
		chomp;
# get list of deflines from cds file in form: >lcl|m1g1|ORF|g.1|m.1|type:5prime_partial|len:236|(+)|edule1_comp18136_c0_seq1_len=1117:3-710(+)
		push(@cds_id_list, $_) if ($_ =~ /^>/); 
	}
	close $cdsfh || die "can't close $cdsfh ";

# parse deflines	
	foreach my $cdsid (@cds_id_list) {
		my @array = split(/\|/, $cdsid); # split defline
		my $transgene = $array[1]; # get transdecoder gene id
		my $tringene = $array[8]; # get trinity gene id including transdecoder cds positions
		my $is_complete = $array[5];
# 				print "$is_complete\n";
		
############################################!!!!!!!!!!!!!!!!!!#############################################################
		# Originally, this next step causes the loss of different genes from the same transcript to be lost -
		# actually named the same thing - and not retrieved from cutfas
		# Need to include the position information in the retrieval from cutfas
		# Changed code to include position info and retrieve only orf sequences
##############################################################################################################################				
		
		if ($tringene =~ m/:(.+)/){; # do not remove positions from trinity gene id but save info in hash - this changed 5/13/14
			$position_hash{$tringene} = join('|', $1, $is_complete); # save position and strand
		}

# make hash of hashes with the first level the index and the second the pairs of trans genes and trinity genes
		$masterhash{$index} = {} unless (exists $masterhash{$index});
		$masterhash{$index}->{$transgene} = $tringene;
	}

# open cutfasta file and create for each index (key) a hash of defline=>seq	from cutfasta 
	my $cutfast = $indexdir . $index . '/chimera_cut_fasta/' . $cutsp . '.gymnos.cutfa';
	open my $cutfastfh, "<$cutfast" || die "can't open $cutfast";
	$cuthash{$index} = &fasta_to_hash($cutfastfh); # A reference to a hash of deflines->seqs for this index
	close $cutfastfh || die "can't close $cutfastfh ";
}


my $mult = '/Users/bdorsey/Documents/Dioon/RNAseq/trinity_files_processed/multiple.orf.seqs_2.txt';

my %multhash;
open my $multfh, "<$mult" or die $!;
while (<$multfh>)
	{
	my $id = (split /\t/)[0];
	$multhash{$id} = 1;
	}
close $multfh or die $!;

# Now we have:
# 1. Hash for each index with transgene keys and tringene values - each held as a value to key:index in masterhash
# 2. Hash of original fasta seqs for each index w/ defline keys ( >edule1_comp18146_c0_seq1_len=1124 ) and seq values - each held as a value to key:index in cuthash
# 3. Hash of abbreviated sp. and indexes
# 4. Hash of index->cutfile sp.
# 5. Hash of tringene => position and strand 
# 6. Hash of $tringene => 1 for all tringenes with multiple orfs - to be excluded


#################################################################################################
# NEW CODE HERE TO USE ONE INDIV PER SPECIES SEQS
#################################################################################################

# PSEUDOCODE

# oneInd.og12_1603: edu2|m2737g2737 mero|m6616g6616 pur1|m2232g2232 sono|m12125g12125 spi2|m4559g4559

# transgene - m1g1
# tringene - edule1_comp18136_c0_seq1_len=1124
# cutspecies - edule1
# og id - edu2|m8821g8821

# parse og id into abbrev. species and transgene
# translate transgene into tringene
# translate abbreviated sp to cutspecies

# get tringene from transgene
# get position and strand info from tringene
# use transgene and position/strand to extract sequence from cuthash

#  foreach og in group file
# 		seqid = elements of array built from ids in line
# 		foreach seqid
# 			split id into abbspecies and transgene
# 			index = spechash{abrev}; tringene = masterhash{index}{transgene}
# 			posinfo = poshash{tringene}; parse pos and strand;
# 			seq = cuthash{index}{tringene}; rc if strand = '-';
# 			print groupfasta seq

# Do we want position info or not?
# This only gives us cds regions from transgenes
# Good for primer design but potentially less informative
# Although, the sequenced loci likely include introns

# NOT USING POSITION INFO BECAUSE WE ARE MINING THESE SEQS FOR MICROSATS
#################################################################################################

# File of og, trin, trans, for multi-orf tringenes
# my $mlog = $dir.'multiOrf.log';
# open my $multilog, ">>$mlog" or die $!;

# Put non-paralog orthogroups into a hash
open my $oglistfh, "<$oglistfile" or die $!;
my @oglist = <$oglistfh>;
close $oglistfh;
chomp(@oglist);

my %oglisthash = ();
foreach (@oglist){
	$oglisthash{$_} = 1;
}

open my $ogfh, "<$ogfile";

# Create one fasta file to contain all seqs.
# We then map seqs to groups using mapSSRog.pl

my $outfile = $outdir. 'allFullSeqs_NoSono.fa';
open my $outfh, ">>$outfile" or die;

LINE: while (<$ogfh>) {
	my @grouparray = split;

# 	Store the orthogroup id
	my $gid = shift @grouparray;
	$gid =~ s/://;

	next LINE unless ($oglisthash{$gid});

# OPENNING FILES HERE CREATES SOME ZERO SIZE FILES FOR MULTI-ORF OGS
# PERHAPS MOVE OPEN TO WITHIN FOREACH LOOP OR USE SYSTEM CALL TO REMOVE ZERO SIZE FILE AT END

# CHANGED TO ONE COMMON FILE ABOVE
# 	Open file to print seqs
# 	my $gfile = $dir.'I12Fasta/nuc/'.$gid.'.fa';
# 	open my $outfh, ">>$gfile" or die $!;
	
# Go through each id in orthogroup and get sequence, print to file	
	SID: foreach my $sid (@grouparray) {
		
# 		Skip sono seqs b/c they are actually E. laevifolius!
		next if ($sid =~ /sono.*/);

# 	Store the abreviated species name and the transdecoder gene id
		my ($absp,$transgene) = split(/\|/, $sid);
		

# 	Store the appropriate index
		my $index = $spechash{$absp};
		my $tringene = $masterhash{$index}->{$transgene};

# 	Get pos/strand/type
		my $exon_info = $position_hash{$tringene};
		my ($start, $end, $strand) = ( $exon_info =~ m/(\d+)-(\d+)\(([\+-])\)/ );
		my $is_complete = ( split /\|/, $exon_info )[1];
		my $length = ($end - $start) + 1;
		$tringene =~ s/:.*//;

# 	Get tringene if not from multi-orf sequence
		if ($multhash{$tringene}) {
			next SID;
		}
		else {

 	# 	Get nuc seq
			my $seq = $cuthash{$index}->{$tringene};
# 			my $fullength = length $seq;
# 			my $seqtoprint = substr($seq, ($start-1), $length); # changed from start +1 to start -1

# 		RC if necessary
			if ($strand eq '-') {
				$seq = &rc($seq);	
			}
	# 	Print			
			my $defline = $sid;

			print $outfh '>'.$defline."\n".$seq."\n";
		}
	}
	
}

# close $multilog or die $!;
close $outfh or die $!;
close $ogfh or die $!;
# close $lengthfh or die $!;


# SUBS
#################################################################

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
		$seq = '';
	}
	else {
		$seq .= $_;
	}
 
}
# Store last sequence
$sequence{$header}=$seq;

return \%sequence;
}

sub rc {
	my ($seq) = @_;
	my $rc = $seq;
	$rc = reverse($rc);
	$rc =~ tr/ATGCatgc/TACGtacg/;
	return($rc);
	}


