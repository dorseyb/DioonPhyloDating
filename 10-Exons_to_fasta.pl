
use strict;
use warnings;


##########################################################################################################################################
# List of loci
# my %list = map {$_ => 1} ('bysp.og20_1739.all.pruned.cds.aln.fa','bysp.og20_169.all.pruned.cds.aln.fa','bysp.og20_741.all.pruned.cds.aln.fa','bysp.og20_1338.all.pruned.cds.aln.fa','bysp.og20_60b.all.pruned.cds.aln.fa','bysp.og20_1086.all.pruned.cds.aln.fa','bysp.og20_315.all.pruned.cds.aln.fa','bysp.og20_101.merged.pruned.cds.aln.fa','bysp.og20_640.all.pruned.cds.aln.fa','bysp.og20_1355.all.pruned.cds.aln.fa','bysp.og20_319.all.pruned.cds.aln.fa','bysp.og20_1132.all.pruned.cds.aln.fa','bysp.og20_100.all.pruned.cds.aln.fa','bysp.og20_1028.all.pruned.cds.aln.fa','bysp.og20_1346.all.pruned.cds.aln.fa','bysp.og20_200.all.pruned.cds.aln.fa','bysp.og20_1830.all.pruned.cds.aln.fa','bysp.og20_1372.all.pruned.cds.aln.fa','bysp.og20_1443.all.pruned.cds.aln.fa','bysp.og20_98a.all.pruned.cds.aln.fa','bysp.og20_1493.all.pruned.cds.aln.fa','bysp.og20_299.all.pruned.cds.aln.fa','bysp.og20_1894.all.pruned.cds.aln.fa','bysp.og20_1529.all.pruned.cds.aln.fa','bysp.og20_1632.all.pruned.cds.aln.fa','bysp.og20_62.all.pruned.cds.aln.fa','bysp.og20_389.all.pruned.cds.aln.fa','bysp.og20_606.all.pruned.cds.aln.fa','bysp.og20_352.all.pruned.cds.aln.fa','bysp.og20_1387.all.pruned.cds.aln.fa','bysp.og20_1000.all.pruned.cds.aln.fa','bysp.og20_1119.all.pruned.cds.aln.fa','bysp.og20_1503.all.pruned.cds.aln.fa','bysp.og20_401.all.pruned.cds.aln.fa','bysp.og20_1552.all.pruned.cds.aln.fa','bysp.og20_28.all.pruned.cds.aln.fa','bysp.og20_1418.all.pruned.cds.aln.fa','bysp.og20_402.all.pruned.cds.aln.fa','bysp.og20_455.all.pruned.cds.aln.fa','bysp.og20_376.all.pruned.cds.aln.fa','bysp.og20_954.all.pruned.cds.aln.fa','bysp.og20_217.all.pruned.cds.aln.fa','bysp.og20_278.all.pruned.cds.aln.fa','bysp.og20_58.all.pruned.cds.aln.fa','bysp.og20_772.all.pruned.cds.aln.fa','bysp.og20_346.all.pruned.cds.aln.fa','bysp.og20_37.all.pruned.cds.aln.fa','bysp.og20_1980.all.pruned.cds.aln.fa','bysp.og20_541.all.pruned.cds.aln.fa','bysp.og20_338.all.pruned.cds.aln.fa','bysp.og20_1566.all.pruned.cds.aln.fa','bysp.og20_996.all.pruned.cds.aln.fa','bysp.og20_1628.all.pruned.cds.aln.fa','bysp.og20_627.all.pruned.cds.aln.fa','bysp.og20_973.all.pruned.cds.aln.fa','bysp.og20_1470.all.pruned.cds.aln.fa','bysp.og20_117.all.pruned.cds.aln.fa','bysp.og20_292.all.pruned.cds.aln.fa','bysp.og20_1454.all.pruned.cds.aln.fa','bysp.og20_1154.all.pruned.cds.aln.fa','bysp.og20_1560.all.pruned.cds.aln.fa','bysp.og20_709.all.pruned.cds.aln.fa','bysp.og20_1553.all.pruned.cds.aln.fa','bysp.og20_80.all.pruned.cds.aln.fa','bysp.og20_415.all.pruned.cds.aln.fa','bysp.og20_109.all.pruned.cds.aln.fa','bysp.og20_807.all.pruned.cds.aln.fa','bysp.og20_1267.all.pruned.cds.aln.fa','bysp.og20_267.all.pruned.cds.aln.fa','bysp.og20_567.all.pruned.cds.aln.fa','bysp.og20_730.all.pruned.cds.aln.fa','bysp.og20_24.all.pruned.cds.aln.fa','bysp.og20_1109.all.pruned.cds.aln.fa','bysp.og20_1095.all.pruned.cds.aln.fa','bysp.og20_1961.all.pruned.cds.aln.fa','bysp.og20_1456.all.pruned.cds.aln.fa','bysp.og20_1588.all.pruned.cds.aln.fa','bysp.og20_1537.all.pruned.cds.aln.fa','bysp.og20_380.all.pruned.cds.aln.fa','bysp.og20_1481.all.pruned.cds.aln.fa','bysp.og20_1207.all.pruned.cds.aln.fa','bysp.og20_55.all.pruned.cds.aln.fa','bysp.og20_956.all.pruned.cds.aln.fa','bysp.og20_1557.all.pruned.cds.aln.fa','bysp.og20_279.all.pruned.cds.aln.fa','bysp.og20_570.all.pruned.cds.aln.fa','bysp.og20_2216.all.pruned.cds.aln.fa','bysp.og20_1156.all.pruned.cds.aln.fa','bysp.og20_1145.all.pruned.cds.aln.fa','bysp.og20_1101.all.pruned.cds.aln.fa','bysp.og20_274.all.pruned.cds.aln.fa','bysp.og20_737.all.pruned.cds.aln.fa','bysp.og20_680.all.pruned.cds.aln.fa','bysp.og20_1069.all.pruned.cds.aln.fa','bysp.og20_2517.all.pruned.cds.aln.fa','bysp.og20_1253.all.pruned.cds.aln.fa','bysp.og20_942.all.pruned.cds.aln.fa','bysp.og20_569.all.pruned.cds.aln.fa','bysp.og20_191.all.pruned.cds.aln.fa','bysp.og20_409.all.pruned.cds.aln.fa','bysp.og20_677.all.pruned.cds.aln.fa','bysp.og20_1153.all.pruned.cds.aln.fa','bysp.og20_509.all.pruned.cds.aln.fa','bysp.og20_81.all.pruned.cds.aln.fa','bysp.og20_1210.all.pruned.cds.aln.fa','bysp.og20_1683.all.pruned.cds.aln.fa','bysp.og20_351.all.pruned.cds.aln.fa','bysp.og20_940.all.pruned.cds.aln.fa','bysp.og20_84.all.pruned.cds.aln.fa','bysp.og20_1141.all.pruned.cds.aln.fa','bysp.og20_8.all.pruned.cds.aln.fa','bysp.og20_251.all.pruned.cds.aln.fa','bysp.og20_38.all.pruned.cds.aln.fa','bysp.og20_1489.all.pruned.cds.aln.fa','bysp.og20_379.all.pruned.cds.aln.fa','bysp.og20_767.all.pruned.cds.aln.fa','bysp.og20_1071.all.pruned.cds.aln.fa','bysp.og20_1703.all.pruned.cds.aln.fa','bysp.og20_476.all.pruned.cds.aln.fa','bysp.og20_1812.all.pruned.cds.aln.fa','bysp.og20_355.all.pruned.cds.aln.fa','bysp.og20_1477.all.pruned.cds.aln.fa','bysp.og20_1225.all.pruned.cds.aln.fa','bysp.og20_689.all.pruned.cds.aln.fa','bysp.og20_223.all.pruned.cds.aln.fa','bysp.og20_92.all.pruned.cds.aln.fa','bysp.og20_1079.all.pruned.cds.aln.fa','bysp.og20_88.all.pruned.cds.aln.fa','bysp.og20_147.all.pruned.cds.aln.fa','bysp.og20_120.all.pruned.cds.aln.fa','bysp.og20_233.all.pruned.cds.aln.fa','bysp.og20_435.all.pruned.cds.aln.fa','bysp.og20_920.all.pruned.cds.aln.fa','bysp.og20_41.all.pruned.cds.aln.fa','bysp.og20_26.all.pruned.cds.aln.fa','bysp.og20_769.all.pruned.cds.aln.fa','bysp.og20_305.all.pruned.cds.aln.fa','bysp.og20_959.all.pruned.cds.aln.fa','bysp.og20_819.all.pruned.cds.aln.fa','bysp.og20_48.all.pruned.cds.aln.fa','bysp.og20_246.all.pruned.cds.aln.fa','bysp.og20_809.all.pruned.cds.aln.fa','bysp.og20_19.all.pruned.cds.aln.fa','bysp.og20_3.all.pruned.cds.aln.fa','bysp.og20_1026.all.pruned.cds.aln.fa','bysp.og20_462.all.pruned.cds.aln.fa','bysp.og20_518.all.pruned.cds.aln.fa','bysp.og20_1309.all.pruned.cds.aln.fa','bysp.og20_997.all.pruned.cds.aln.fa','bysp.og20_14.all.pruned.cds.aln.fa','bysp.og20_1666.all.pruned.cds.aln.fa','bysp.og20_2146.all.pruned.cds.aln.fa','bysp.og20_133.all.pruned.cds.aln.fa','bysp.og20_885.all.pruned.cds.aln.fa','bysp.og20_203.all.pruned.cds.aln.fa','bysp.og20_682.all.pruned.cds.aln.fa','bysp.og20_755.all.pruned.cds.aln.fa','bysp.og20_263.all.pruned.cds.aln.fa','bysp.og20_342.all.pruned.cds.aln.fa','bysp.og20_304.all.pruned.cds.aln.fa','bysp.og20_15.all.pruned.cds.aln.fa','bysp.og20_1016.all.pruned.cds.aln.fa','bysp.og20_294.all.pruned.cds.aln.fa','bysp.og20_306.all.pruned.cds.aln.fa','bysp.og20_249.all.pruned.cds.aln.fa','bysp.og20_992.all.pruned.cds.aln.fa','bysp.og20_17.all.pruned.cds.aln.fa','bysp.og20_816.all.pruned.cds.aln.fa','bysp.og20_566.all.pruned.cds.aln.fa','bysp.og20_31.all.pruned.cds.aln.fa','bysp.og20_34.all.pruned.cds.aln.fa','bysp.og20_2272.all.pruned.cds.aln.fa','bysp.og20_316.all.pruned.cds.aln.fa','bysp.og20_780.all.pruned.cds.aln.fa','bysp.og20_793.all.pruned.cds.aln.fa','bysp.og20_272.all.pruned.cds.aln.fa','bysp.og20_679.all.pruned.cds.aln.fa','bysp.og20_124.edit.pruned.cds.aln.fa','bysp.og20_1930.all.pruned.cds.aln.fa','bysp.og20_584.all.pruned.cds.aln.fa','bysp.og20_1970.all.pruned.cds.aln.fa','bysp.og20_2186.all.pruned.cds.aln.fa','bysp.og20_908.all.pruned.cds.aln.fa','bysp.og20_1235.all.pruned.cds.aln.fa','bysp.og20_1031.all.pruned.cds.aln.fa','bysp.og20_1713.all.pruned.cds.aln.fa','bysp.og20_1297.all.pruned.cds.aln.fa','bysp.og20_2392.all.pruned.cds.aln.fa','bysp.og20_50.all.pruned.cds.aln.fa','bysp.og20_11.all.pruned.cds.aln.fa','bysp.og20_663.all.pruned.cds.aln.fa','bysp.og20_337.all.pruned.cds.aln.fa','bysp.og20_1068.all.pruned.cds.aln.fa','bysp.og20_665.all.pruned.cds.aln.fa','bysp.og20_207.all.pruned.cds.aln.fa','bysp.og20_461.all.pruned.cds.aln.fa','bysp.og20_221.all.pruned.cds.aln.fa','bysp.og20_322.all.pruned.cds.aln.fa','bysp.og20_426.all.pruned.cds.aln.fa','bysp.og20_944.all.pruned.cds.aln.fa','bysp.og20_27.all.pruned.cds.aln.fa','bysp.og20_672.all.pruned.cds.aln.fa','bysp.og20_947.all.pruned.cds.aln.fa','bysp.og20_348.all.pruned.cds.aln.fa','bysp.og20_660.all.pruned.cds.aln.fa','bysp.og20_552.all.pruned.cds.aln.fa','bysp.og20_932.all.pruned.cds.aln.fa','bysp.og20_4.all.pruned.cds.aln.fa','bysp.og20_116.all.pruned.cds.aln.fa','bysp.og20_1172.all.pruned.cds.aln.fa','bysp.og20_1394.all.pruned.cds.aln.fa','bysp.og20_33.all.pruned.cds.aln.fa','bysp.og20_496.all.pruned.cds.aln.fa','bysp.og20_190.all.pruned.cds.aln.fa','bysp.og20_1142.all.pruned.cds.aln.fa','bysp.og20_1295.all.pruned.cds.aln.fa','bysp.og20_1219.all.pruned.cds.aln.fa','bysp.og20_181.all.pruned.cds.aln.fa','bysp.og20_284.all.pruned.cds.aln.fa','bysp.og20_1474.all.pruned.cds.aln.fa','bysp.og20_1935.all.pruned.cds.aln.fa','bysp.og20_1201.all.pruned.cds.aln.fa','bysp.og20_73.all.pruned.cds.aln.fa','bysp.og20_617.all.pruned.cds.aln.fa','bysp.og20_1292.all.pruned.cds.aln.fa','bysp.og20_930.all.pruned.cds.aln.fa','bysp.og20_1332.all.pruned.cds.aln.fa','bysp.og20_801.all.pruned.cds.aln.fa','bysp.og20_452.all.pruned.cds.aln.fa','bysp.og20_630.all.pruned.cds.aln.fa','bysp.og20_336.all.pruned.cds.aln.fa','bysp.og20_703.all.pruned.cds.aln.fa','bysp.og20_77.all.pruned.cds.aln.fa','bysp.og20_2493.all.pruned.cds.aln.fa','bysp.og20_1480.all.pruned.cds.aln.fa','bysp.og20_13.all.pruned.cds.aln.fa','bysp.og20_244.all.pruned.cds.aln.fa','bysp.og20_9.all.pruned.cds.aln.fa','bysp.og20_914.all.pruned.cds.aln.fa','bysp.og20_583.all.pruned.cds.aln.fa','bysp.og20_618.all.pruned.cds.aln.fa','bysp.og20_1490.all.pruned.cds.aln.fa','bysp.og20_150.all.pruned.cds.aln.fa','bysp.og20_70.all.pruned.cds.aln.fa','bysp.og20_1244.all.pruned.cds.aln.fa','bysp.og20_276.all.pruned.cds.aln.fa','bysp.og20_69.all.pruned.cds.aln.fa','bysp.og20_817.all.pruned.cds.aln.fa','bysp.og20_2330.all.pruned.cds.aln.fa','bysp.og20_1508.all.pruned.cds.aln.fa','bysp.og20_2382.all.pruned.cds.aln.fa','bysp.og20_659.all.pruned.cds.aln.fa','bysp.og20_1308.all.pruned.cds.aln.fa','bysp.og20_1383.all.pruned.cds.aln.fa','bysp.og20_591.all.pruned.cds.aln.fa','bysp.og20_1774.all.pruned.cds.aln.fa','bysp.og20_107.all.pruned.cds.aln.fa','bysp.og20_2098.all.pruned.cds.aln.fa','bysp.og20_1679.all.pruned.cds.aln.fa','bysp.og20_1896.all.pruned.cds.aln.fa','bysp.og20_175.all.pruned.cds.aln.fa','bysp.og20_1080.all.pruned.cds.aln.fa','bysp.og20_544.all.pruned.cds.aln.fa','bysp.og20_2046.all.pruned.cds.aln.fa','bysp.og20_477.all.pruned.cds.aln.fa','bysp.og20_103.all.pruned.cds.aln.fa','bysp.og20_5.all.pruned.cds.aln.fa','bysp.og20_568.all.pruned.cds.aln.fa','bysp.og20_1890.all.pruned.cds.aln.fa','bysp.og20_924.all.pruned.cds.aln.fa','bysp.og20_57.all.pruned.cds.aln.fa','bysp.og20_177.all.pruned.cds.aln.fa','bysp.og20_377.all.pruned.cds.aln.fa','bysp.og20_216.all.pruned.cds.aln.fa','bysp.og20_794.all.pruned.cds.aln.fa','bysp.og20_777.all.pruned.cds.aln.fa','bysp.og20_1665.all.pruned.cds.aln.fa','bysp.og20_1732.all.pruned.cds.aln.fa','bysp.og20_622.all.pruned.cds.aln.fa','bysp.og20_51.all.pruned.cds.aln.fa','bysp.og20_972.all.pruned.cds.aln.fa','bysp.og20_259.all.pruned.cds.aln.fa','bysp.og20_937.all.pruned.cds.aln.fa');

# my $alldir = '/Users/bdorsey/Documents/Dioon/RNAseq/trinity_files_processed/ortho_search/ortho/bySpeciesResults/ogfasta/I20/all/';
# my $fasta = $alldir.'nuc_ortholog_fasta/orthogroup_cons_alignments/orthogroups_cons.fa';

# From iupac cons
# my $fasta = $alldir.'nuc_ortholog_fasta/orthogroup_cons_alignments/orthogroups_cons.fa';

# From findsplice.pl using majcon
# my $gff = '/Users/bdorsey/Documents/Dioon/RNAseq/trinity_files_processed/ortho_search/ortho/bySpeciesResults/ogfasta/I20/all/nuc_ortholog_fasta/cds.majcon.ss.exon_gt199.gff';

# From majcons
# my $fasta = $alldir.'nuc_ortholog_fasta/orthogroup_majcons_alignments/orthogroups_maj_cons.fa';
# my $gff = $alldir.'gmap_files/cds_ptaeda_majcon_exons_to_orthogroups.gff';
# my $newfasta = $alldir. 'nuc_ortholog_fasta/orthogroup_majcons_alignments/cds_og_majcons_ptaeda_exons_gt199.fa';

# my $newfasta = $alldir.'nuc_ortholog_fasta/orthogroup_cons_alignments/all_orthogroups_cons_exons_gt199.fa';
# my $gff = $alldir.'gmap_files/map_cds_exons_to_orthogroups_sorted.gff';

# my $gff = $alldir.'nuc_ortholog_fasta/cds.nomap.ss.exon.gff';
# my $newfasta = $alldir.'nuc_ortholog_fasta/orthogroup_cons_alignments/nomap_og_cons_ss.exons.fa';

# my $newgff = '/Users/bdorsey/Documents/Dioon/RNAseq/trinity_files_processed/probedesign/exons.utr..gff';

# From ptaeda/pabies intersection
# my $probedir = '/Users/bdorsey/Documents/Dioon/RNAseq/trinity_files_processed/probedesign/';
# my $gff = $probedir. 'majcon_exons_intersect_and_non.gff';

# 7-13-15 adding iupac b/c we want snps marked for aligning reads
# my $newfasta = $probedir. 'majcons_exons_intersect_and_non_iupac.fa';

# printing 'exons' from majcon splice site search using iupac con seqs
# my $newfasta = $probedir . 'og_cons_nomap_exons_from_majcon_ss_search.fa';
##########################################################################################################################################


# A script to read a gff file of inferred exon positions and print the exon sequences to file. Takes a gff file and a fasta file of
# original sequences as input. 
# Can output fasta file and/or a primer3 input file.
############################################################################



# For pcr-based experiment

my $dir = '/Users/bdorsey/Documents/Dioon/RNAseq/trinity_files_processed/ortho_search/ortho/oneIndPerSpecies/I12Fasta/';

# my $gff = $dir.'gmap/target_pabies_exons_to_orthogroups.gff';
# my $gff = $dir.'gmap/exonsForIntronsGT50.gff';

my $gff = $dir.'gmap/blastHit.abiesMap.targets.gff';

my $fasta = $dir.'nuc/blastHit.conseqs.thresholdTrim.fa';

my $newfasta = $dir.'exons/targetExons300-800.fa';

my $p3file = $dir.'exons/targetExons300-800.p3';

open my $gffh, "<$gff" or die $!;
open my $fastafh, "<$fasta" or die $!;
open my $newfh, ">$newfasta" or die $!;

my $fastref = &fasta_to_hash($fastafh);
my %fasthash = %$fastref;

# go through each line in gff
# capture seq name and exon coordinates
# if it is the first exon cut positions 0-start of exon and start to end
# if its not the first or last cut from start to end positions
# if its the last cut start and end and then end to last base
# print each of these exons to file with the defline <original_def>_exon_number/pre/post


# Changed defline to use gmap exon id 10-30-15


my ($last_end, $seq, $seqlen, $total, $passed);
my $current_id = '';
my $exon_num = 0;
my $first_gene = 1;
my $cutoff = 299;
my $max = 801;
my $p3 =1;

open my $p3fh, ">$p3file" or die $! if $p3;

LINE: while (my $line = <$gffh>) {
	next LINE if ($line =~ /#/);
	
	my @splt = split /\t/, $line;
	my $id = $splt[0];
	my ($start, $end) = @splt[3,4];
	my $length = ($end - $start) + 1;
	my $gmapexon = (split(/;/, $splt[8]))[0];
	$gmapexon =~ s/ID=//;

# if the og is in the set of loci
# 	if ($list{$id}){

# 	if we are in the middle of a gene
		if ($id eq $current_id) {
# 		increase exon number
			$exon_num++;

# 		get next exon seq
			my $cut = uc(substr($seq, 0, $length, ''));
			$total++ if ($cut);
			if ( (length($cut) > $cutoff) && (length($cut) < $max) ){
# 				print $newfh '>'.$id.'_exon'.$exon_num.'_'.$start.'-'."$end\n$cut\n";
				print $newfh '>'.$gmapexon.'_'.$start.'-'.$end."\n".$cut."\n";
				if ($p3) {
					print $p3fh 'Sequence_ID='.$gmapexon.'_'.$start.'-'.$end."\nSEQUENCE_TEMPLATE=".$cut."\n=\n";
				}	
				$passed++;
			}
# Testing
# 			else {
# 				print "$gmapexon.Firstblock\n";
# 			}
		}
		
# 	if we are at a new gene	
		else {

# 		new gene, not first gene
			if ($first_gene == 0){

# 			print final seq from last gene
				if ($seq){
					$total++;
					
					if ( (length($seq) > $cutoff) && (length($seq) < $max) ){
					my $n = $last_end+1;
					print $newfh '>'.$current_id.'_end_'.$n.'-'."$seqlen\n$seq\n";
					if ($p3){
						print $p3fh 'SEQUENCE_ID='.$current_id.'_end_'.$n.'-'."$seqlen\n".'SEQUENCE_TEMPLATE='."$seq\n=\n";
					}
					$passed++;
					}
# Testing
# 					else {
# 						my $n = $last_end+1;
# 						print $current_id.'_end_'.$n.'-'."$seqlen\n$seq\n";
# 					}

				}

				$seq = uc($fasthash{$id});
				$seqlen = length($seq);
				$exon_num = 0;

# 			get first seq from new gene (before first exon)	
				unless ($start == 1){
					my $pre = uc(substr($seq, 0, ($start-1), ''));
					$total++ if ($pre);
					if ( (length($pre) > $cutoff) && (length($pre) < $max) ){
						my $m = $start-1;
						print $newfh '>'."$id".'_pre_1-'."$m\n$pre\n";
						if ($p3) {
							print $p3fh 'SEQUENCE_ID='."$id".'_pre_1-'."$m\n".'SEQUENCE_TEMPLATE='."$pre\n=\n";
						}
						$passed++ ;
					}
# Testing
# 					else {
# 						my $m = $start-1;
# 						print "$id".'_pre_1-'."$m\n$pre\n";
# 					}
					
				}
# 			get first exon
				my $exon = uc( substr($seq, 0, $length, '') );
				$total++ if ($exon);
				$exon_num++;
				if ( (length($exon) > $cutoff) && (length($exon) < $max) ){
# 					print $newfh '>'.$id.'_exon'.$exon_num.'_'.$start.'-'."$end\n$exon\n";
					print $newfh '>'.$gmapexon.'_'.$start.'-'.$end."\n".$exon."\n";
					if ($p3) {
						print $p3fh 'SEQUENCE_ID='.$gmapexon.'_'.$start.'-'.$end."\n".'SEQUENCE_TEMPLATE='.$exon."\n=\n";
					}
					$passed++;
				}
# Testing
# 				else {
# 					print "$gmapexon.Fourthblock\n";
# 				}

			}
			
# 		first time through
			elsif ($first_gene == 1) {

				my $m = $start-1;
				$seq = uc( $fasthash{$id} );
				$seqlen = length($seq);
				$exon_num = 0;
				unless ($start == 1){
					my $pre = uc( substr($seq, 0, ($start-1), '') );
					$total++ if ($pre);
					if ( (length($pre) > $cutoff) && (length($pre) < $max) ){
						print $newfh '>'."$id".'_pre_1-'."$m\n$pre\n";
						if ($p3) {
							print $p3fh 'SEQUENCE_ID='."$id".'_pre_1-'."$m\n".'SEQUENCE_TEMPLATE='."$pre\n=\n";
						}
						$passed++;
					}
# Testing
# 					else {
# 						print "$gmapexon.Fifthblock\n";
# 					}
	
				}
				my $exon = substr($seq, 0, $length, '');
				$total++ if ($exon);
				$exon_num++;
				if ( (length($exon) > $cutoff) && (length($exon) < $max) ){
# 					print $newfh '>'.$id.'_exon'.$exon_num.'_'.$start.'-'."$end\n$exon\n";
					print $newfh '>'.$gmapexon.'_'.$start.'-'.$end."\n".$exon."\n";
					if ($p3) {
						print $p3fh 'SEQUENCE_ID='.$gmapexon.'_'.$start.'-'.$end."\n".'SEQUENCE_TEMPLATE='.$exon."\n=\n";
					}
					$passed++ ;
				}
# Testing			
# 				else {
# 					print "$gmapexon.Sixthblock\n";
# 				}

				$first_gene = 0;
			}
		}	
	
		$current_id = $id;
		$last_end = $end;
		next LINE;
}

close $gffh or die $!;
close $fastafh or die $!;
close $newfh or die $!;
close $p3fh or die $! if $p3;

print "Total = $total\n";
print "Passed cutoff = $passed\n";



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
