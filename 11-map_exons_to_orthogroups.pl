use strict;
use warnings;

# Convert gff files from gmap of orthogroups to several genomes into gff file mapping exons to 'genome' of all orthogroup cons
# Updated 10-22-15 to process targets seqs from pcr-based capture experiment

# my $gffdir = '/Users/bdorsey/Documents/Dioon/RNAseq/trinity_files_processed/ortho_search/ortho/bySpeciesResults/ogfasta/I20/all/gmap_files/';
my $gffdir = '/Users/bdorsey/Documents/Dioon/RNAseq/trinity_files_processed/ortho_search/ortho/oneIndPerSpecies/I12Fasta/gmap/';

# my $gff = $gffdir.'targetsToPabiesAll.gff';
# my $gff = $gffdir.'targetsToPtaedaAll.gff';
# my $gff = $gffdir.'no_overlap_pabiesToptaeda.gff';
my $gff = $gffdir.'gmap.11-3-15_blastpHits_thresh_to_pabies.uniq.gff';

# my $gff = $gffdir.'10-6-14.og_majcons_to_pabies_cat.sorted.gff';
# my $gff = $gffdir.'10-6-14.og_majcons_to_ptaeda_cat.sorted.gff';
# my $gff = $gffdir.'9-17-14.og_cons_to_ptaeda_cat_sorted.gff3';
# my $gff = $gffdir.'8-25-14.og_cons_to_Pabies_nocontam_cat.gff3';
# my $gff = $gffdir.'8-26-14.ws_cons_to_Pabies_cat.sorted.gff3';

# my $ngff = $gffdir.'target_pabies_exons_to_orthogroups.gff';
# my $ngff = $gffdir.'target_ptaeda_exons_to_orthogroups.gff';
# my $ngff = $gffdir.'no_overlap_pabiesToptaeda_targets.gff';
my $ngff = $gffdir.'blastHit.abiesMap.targets.gff';

# my $ngff = $gffdir. 'cds_pabies_majcon_exons_to_orthogroups.gff';
# my $ngff = $gffdir. 'cds_ptaeda_majcon_exons_to_orthogroups.gff';
# my $ngff = $gffdir . 'cds_ptaeda_exons_map_to_orthogroups.gff3';
# my $ngff = $gffdir . 'cds_exons_map_to_orthogroups_gt199.gff';
# my $ngff = $gffdir . 'map_ws_exons_to_orthogroups.gff';
# my $ngff = $gffdir. 'cds_pmume_exons_map_to_orthogroups.gff';

my $tfile =  '/Users/bdorsey/Documents/Dioon/RNAseq/trinity_files_processed/ortho_search/ortho/oneIndPerSpecies/I12Fasta/blastHitAbiesMapTargetIds.txt';
open my $tfh, "<$tfile" or die $!;
my @targets = <$tfh>;
chomp @targets;
close $tfh or die;

open my $newgff, ">$ngff" or die $!;
print $newgff "##gff-version   3\n#Conversion of exons mapped to P. abies genome to position on trimmed orthogroup cons\n";

open my $gffh, "<$gff" or die $!;

my $cutoff = 1; #199;
my $counter = 0;
LINE: while (my $line = <$gffh>) {
	if ($line =~ /exon/) {
		chomp $line;
		my @splt = split /\t/, $line;
		my $genome = $splt[1];
		my $strand = $splt[6];
		my $info = $splt[8];
# 		$line =~ /Target=(bysp\.og20_.+?)\s(\d+)\s(\d+)\s([+-])$/;
		$line =~ /Target=(oneInd.og12_.+?)\s(\d+)\s(\d+)\s([+-])$/;
# 		my ($target, $start, $end, $strand) = ($1, $2, $3, $4);
		my ($target, $start, $end) = ($1, $2, $3);		
		my $length = ($end-$start) + 1;
# 		my $gtarget = $target;
# 		$gtarget =~ s/\.afa//;
		next LINE unless (grep /$target/, @targets);
		$counter++;
		print "\r$counter";
		if ($length > $cutoff){
			print $newgff "$target\t$genome\texon\t$start\t$end\t$splt[5]\t$strand\t\.\t$info;Length=$length\n";
		}
	}
# 		elseif ($line =~ /mRNA/) {
# 			chomp $line;
# 			my @splt = split /\t/, $line;
# 			my $genome = $splt[1];
# 			my $info = $splt[8];
# 			$line =~ /Target=(bysp\.og20_.+?)\s(\d+)\s(\d+)\s([+-])$/;
# 			my ($target, $start, $end, $strand) = ($1, $2, $3, $4);
# 			print $newgff "$target\t$genome\tmRNA\t$start\t$end\t$splt[5]\t$strand\t\.\t$info\n";
# 			
# 		}
}
print"\n";
close $gffh or die $!;
close $newgff or die $!;


__DATA__
MA_5	Pabies.nocontam.gmapdb	gene	32797	33304	.	+	.	ID=bysp.og20_29.all.pruned.cds.aln.fa.path1;Name=bysp.og20_29.all.pruned.cds.aln.fa
MA_5	Pabies.nocontam.gmapdb	mRNA	32797	33304	.	+	.	ID=bysp.og20_29.all.pruned.cds.aln.fa.mrna1;Name=bysp.og20_29.all.pruned.cds.aln.fa;Parent=bysp.og20_29.all.pruned.cds.aln.fa.path1;coverage=31.1;identity=80.9;matches=342;mismatches=81;indels=0;unknowns=85
MA_5	Pabies.nocontam.gmapdb	exon	32797	33304	80	+	.	ID=bysp.og20_29.all.pruned.cds.aln.fa.mrna1.exon1;Name=bysp.og20_29.all.pruned.cds.aln.fa;Parent=bysp.og20_29.all.pruned.cds.aln.fa.mrna1;Target=bysp.og20_29.all.pruned.cds.aln.fa 503 1010 +
MA_5	Pabies.nocontam.gmapdb	CDS	32797	33303	80	+	0	ID=bysp.og20_29.all.pruned.cds.aln.fa.mrna1.cds1;Name=bysp.og20_29.all.pruned.cds.aln.fa;Parent=bysp.og20_29.all.pruned.cds.aln.fa.mrna1;Target=bysp.og20_29.all.pruned.cds.aln.fa 503 1009 +
MA_38	Pabies.nocontam.gmapdb	gene	5379	7545	.	+	.	ID=bysp.og20_817.all.pruned.cds.aln.fa.path1;Name=bysp.og20_817.all.pruned.cds.aln.fa
MA_38	Pabies.nocontam.gmapdb	mRNA	5379	7545	.	+	.	ID=bysp.og20_817.all.pruned.cds.aln.fa.mrna1;Name=bysp.og20_817.all.pruned.cds.aln.fa;Parent=bysp.og20_817.all.pruned.cds.aln.fa.path1;coverage=84.3;identity=83.0;matches=1058;mismatches=193;indels=23;unknowns=59
MA_38	Pabies.nocontam.gmapdb	exon	5379	6224	79	+	.	ID=bysp.og20_817.all.pruned.cds.aln.fa.mrna1.exon1;Name=bysp.og20_817.all.pruned.cds.aln.fa;Parent=bysp.og20_817.all.pruned.cds.aln.fa.mrna1;Target=bysp.og20_817.all.pruned.cds.aln.fa 245 1105 +
MA_38	Pabies.nocontam.gmapdb	CDS	5380	5700	76	+	0	ID=bysp.og20_817.all.pruned.cds.aln.fa.mrna1.cds1;Name=bysp.og20_817.all.pruned.cds.aln.fa;Parent=bysp.og20_817.all.pruned.cds.aln.fa.mrna1;Target=bysp.og20_817.all.pruned.cds.aln.fa 246 578 +
MA_38	Pabies.nocontam.gmapdb	CDS	5702	5711	75	+	0	ID=bysp.og20_817.all.pruned.cds.aln.fa.mrna1.cds2;Name=bysp.og20_817.all.pruned.cds.aln.fa;Parent=bysp.og20_817.all.pruned.cds.aln.fa.mrna1;Target=bysp.og20_817.all.pruned.cds.aln.fa 579 588 +
MA_38	Pabies.nocontam.gmapdb	CDS	5714	5867	78	+	1	ID=bysp.og20_817.all.pruned.cds.aln.fa.mrna1.cds3;Name=bysp.og20_817.all.pruned.cds.aln.fa;Parent=bysp.og20_817.all.pruned.cds.aln.fa.mrna1;Target=bysp.og20_817.all.pruned.cds.aln.fa 589 743 +
MA_38	Pabies.nocontam.gmapdb	CDS	5869	6224	82	+	0	ID=bysp.og20_817.all.pruned.cds.aln.fa.mrna1.cds4;Name=bysp.og20_817.all.pruned.cds.aln.fa;Parent=bysp.og20_817.all.pruned.cds.aln.fa.mrna1;Target=bysp.og20_817.all.pruned.cds.aln.fa 744 1105 +
MA_38	Pabies.nocontam.gmapdb	exon	6494	6720	86	+	.	ID=bysp.og20_817.all.pruned.cds.aln.fa.mrna1.exon2;Name=bysp.og20_817.all.pruned.cds.aln.fa;Parent=bysp.og20_817.all.pruned.cds.aln.fa.mrna1;Target=bysp.og20_817.all.pruned.cds.aln.fa 1106 1332 +
MA_38	Pabies.nocontam.gmapdb	CDS	6494	6720	86	+	2	ID=bysp.og20_817.all.pruned.cds.aln.fa.mrna1.cds5;Name=bysp.og20_817.all.pruned.cds.aln.fa;Parent=bysp.og20_817.all.pruned.cds.aln.fa.mrna1;Target=bysp.og20_817.all.pruned.cds.aln.fa 1106 1332 +
MA_38	Pabies.nocontam.gmapdb	exon	6821	6962	93	+	.	ID=bysp.og20_817.all.pruned.cds.aln.fa.mrna1.exon3;Name=bysp.og20_817.all.pruned.cds.aln.fa;Parent=bysp.og20_817.all.pruned.cds.aln.fa.mrna1;Target=bysp.og20_817.all.pruned.cds.aln.fa 1333 1474 +
MA_38	Pabies.nocontam.gmapdb	CDS	6821	6962	93	+	1	ID=bysp.og20_817.all.pruned.cds.aln.fa.mrna1.cds6;Name=bysp.og20_817.all.pruned.cds.aln.fa;Parent=bysp.og20_817.all.pruned.cds.aln.fa.mrna1;Target=bysp.og20_817.all.pruned.cds.aln.fa 1333 1474 +
MA_38	Pabies.nocontam.gmapdb	exon	7261	7322	90	+	.	ID=bysp.og20_817.all.pruned.cds.aln.fa.mrna1.exon4;Name=bysp.og20_817.all.pruned.cds.aln.fa;Parent=bysp.og20_817.all.pruned.cds.aln.fa.mrna1;Target=bysp.og20_817.all.pruned.cds.aln.fa 1475 1536 +
MA_38	Pabies.nocontam.gmapdb	CDS	7261	7322	90	+	2	ID=bysp.og20_817.all.pruned.cds.aln.fa.mrna1.cds7;Name=bysp.og20_817.all.pruned.cds.aln.fa;Parent=bysp.og20_817.all.pruned.cds.aln.fa.mrna1;Target=bysp.og20_817.all.pruned.cds.aln.fa 1475 1536 +
MA_38	Pabies.nocontam.gmapdb	exon	7509	7545	94	+	.	ID=bysp.og20_817.all.pruned.cds.aln.fa.mrna1.exon5;Name=bysp.og20_817.all.pruned.cds.aln.fa;Parent=bysp.og20_817.all.pruned.cds.aln.fa.mrna1;Target=bysp.og20_817.all.pruned.cds.aln.fa 1537 1573 +
MA_38	Pabies.nocontam.gmapdb	CDS	7509	7543	94	+	1	ID=bysp.og20_817.all.pruned.cds.aln.fa.mrna1.cds8;Name=bysp.og20_817.all.pruned.cds.aln.fa;Parent=bysp.og20_817.all.pruned.cds.aln.fa.mrna1;Target=bysp.og20_817.all.pruned.cds.aln.fa 1537 1571 +
MA_359	Pabies.nocontam.gmapdb	gene	2648	5516	.	-	.	ID=bysp.og20_574.all.pruned.cds.aln.fa.path1;Name=bysp.og20_574.all.pruned.cds.aln.fa
MA_359	Pabies.nocontam.gmapdb	mRNA	2648	5516	.	-	.	ID=bysp.og20_574.all.pruned.cds.aln.fa.mrna1;Name=bysp.og20_574.all.pruned.cds.aln.fa;Parent=bysp.og20_574.all.pruned.cds.aln.fa.path1;coverage=29.0;identity=86.8;matches=387;mismatches=59;indels=0;unknowns=3
MA_359	Pabies.nocontam.gmapdb	exon	2648	2726	87	-	.	ID=bysp.og20_574.all.pruned.cds.aln.fa.mrna1.exon6;Name=bysp.og20_574.all.pruned.cds.aln.fa;Parent=bysp.og20_574.all.pruned.cds.aln.fa.mrna1;Target=bysp.og20_574.all.pruned.cds.aln.fa 695 773 +
MA_359	Pabies.nocontam.gmapdb	CDS	2650	2726	86	-	1	ID=bysp.og20_574.all.pruned.cds.aln.fa.mrna1.cds6;Name=bysp.og20_574.all.pruned.cds.aln.fa;Parent=bysp.og20_574.all.pruned.cds.aln.fa.mrna1;Target=bysp.og20_574.all.pruned.cds.aln.fa 695 771 +
MA_359	Pabies.nocontam.gmapdb	exon	4005	4080	86	-	.	ID=bysp.og20_574.all.pruned.cds.aln.fa.mrna1.exon5;Name=bysp.og20_574.all.pruned.cds.aln.fa;Parent=bysp.og20_574.all.pruned.cds.aln.fa.mrna1;Target=bysp.og20_574.all.pruned.cds.aln.fa 619 694 +
MA_359	Pabies.nocontam.gmapdb	CDS	4005	4080	86	-	0	ID=bysp.og20_574.all.pruned.cds.aln.fa.mrna1.cds5;Name=bysp.og20_574.all.pruned.cds.aln.fa;Parent=bysp.og20_574.all.pruned.cds.aln.fa.mrna1;Target=bysp.og20_574.all.pruned.cds.aln.fa 619 694 +
MA_359	Pabies.nocontam.gmapdb	exon	4402	4476	82	-	.	ID=bysp.og20_574.all.pruned.cds.aln.fa.mrna1.exon4;Name=bysp.og20_574.all.pruned.cds.aln.fa;Parent=bysp.og20_574.all.pruned.cds.aln.fa.mrna1;Target=bysp.og20_574.all.pruned.cds.aln.fa 544 618 +
MA_359	Pabies.nocontam.gmapdb	CDS	4402	4476	82	-	0	ID=bysp.og20_574.all.pruned.cds.aln.fa.mrna1.cds4;Name=bysp.og20_574.all.pruned.cds.aln.fa;Parent=bysp.og20_574.all.pruned.cds.aln.fa.mrna1;Target=bysp.og20_574.all.pruned.cds.aln.fa 544 618 +
MA_359	Pabies.nocontam.gmapdb	exon	4646	4711	80	-	.	ID=bysp.og20_574.all.pruned.cds.aln.fa.mrna1.exon3;Name=bysp.og20_574.all.pruned.cds.aln.fa;Parent=bysp.og20_574.all.pruned.cds.aln.fa.mrna1;Target=bysp.og20_574.all.pruned.cds.aln.fa 478 543 +
MA_359	Pabies.nocontam.gmapdb	CDS	4646	4711	80	-	0	ID=bysp.og20_574.all.pruned.cds.aln.fa.mrna1.cds3;Name=bysp.og20_574.all.pruned.cds.aln.fa;Parent=bysp.og20_574.all.pruned.cds.aln.fa.mrna1;Target=bysp.og20_574.all.pruned.cds.aln.fa 478 543 +
MA_359	Pabies.nocontam.gmapdb	exon	5008	5070	90	-	.	ID=bysp.og20_574.all.pruned.cds.aln.fa.mrna1.exon2;Name=bysp.og20_574.all.pruned.cds.aln.fa;Parent=bysp.og20_574.all.pruned.cds.aln.fa.mrna1;Target=bysp.og20_574.all.pruned.cds.aln.fa 415 477 +
MA_359	Pabies.nocontam.gmapdb	CDS	5008	5070	90	-	0	ID=bysp.og20_574.all.pruned.cds.aln.fa.mrna1.cds2;Name=bysp.og20_574.all.pruned.cds.aln.fa;Parent=bysp.og20_574.all.pruned.cds.aln.fa.mrna1;Target=bysp.og20_574.all.pruned.cds.aln.fa 415 477 +
MA_359	Pabies.nocontam.gmapdb	exon	5427	5516	92	-	.	ID=bysp.og20_574.all.pruned.cds.aln.fa.mrna1.exon1;Name=bysp.og20_574.all.pruned.cds.aln.fa;Parent=bysp.og20_574.all.pruned.cds.aln.fa.mrna1;Target=bysp.og20_574.all.pruned.cds.aln.fa 325 414 +
MA_359	Pabies.nocontam.gmapdb	CDS	5427	5516	91	-	0	ID=bysp.og20_574.all.pruned.cds.aln.fa.mrna1.cds1;Name=bysp.og20_574.all.pruned.cds.aln.fa;Parent=bysp.og20_574.all.pruned.cds.aln.fa.mrna1;Target=bysp.og20_574.all.pruned.cds.aln.fa 325 414 +
MA_1249	Pabies.nocontam.gmapdb	gene	34189	35185	.	+	.	ID=bysp.og20_3056.all.pruned.cds.aln.fa.path1;Name=bysp.og20_3056.all.pruned.cds.aln.fa
MA_1249	Pabies.nocontam.gmapdb	mRNA	34189	35185	.	+	.	ID=bysp.og20_3056.all.pruned.cds.aln.fa.mrna1;Name=bysp.og20_3056.all.pruned.cds.aln.fa;Parent=bysp.og20_3056.all.pruned.cds.aln.fa.path1;coverage=40.6;identity=84.3;matches=697;mismatches=121;indels=9;unknowns=26
MA_1249	Pabies.nocontam.gmapdb	exon	34189	34932	83	+	.	ID=bysp.og20_3056.all.pruned.cds.aln.fa.mrna1.exon1;Name=bysp.og20_3056.all.pruned.cds.aln.fa;Parent=bysp.og20_3056.all.pruned.cds.aln.fa.mrna1;Target=bysp.og20_3056.all.pruned.cds.aln.fa 1241 1981 -
MA_1249	Pabies.nocontam.gmapdb	CDS	34189	34389	85	+	0	ID=bysp.og20_3056.all.pruned.cds.aln.fa.mrna1.cds1;Name=bysp.og20_3056.all.pruned.cds.aln.fa;Parent=bysp.og20_3056.all.pruned.cds.aln.fa.mrna1;Target=bysp.og20_3056.all.pruned.cds.aln.fa 1781 1981 -
MA_1249	Pabies.nocontam.gmapdb	exon	35080	35185	87	+	.	ID=bysp.og20_3056.all.pruned.cds.aln.fa.mrna1.exon2;Name=bysp.og20_3056.all.pruned.cds.aln.fa;Parent=bysp.og20_3056.all.pruned.cds.aln.fa.mrna1;Target=bysp.og20_3056.all.pruned.cds.aln.fa 1135 1240 -
MA_1692	Pabies.nocontam.gmapdb	gene	24487	26724	.	+	.	ID=bysp.og20_1034.all.pruned.cds.aln.fa.path1;Name=bysp.og20_1034.all.pruned.cds.aln.fa
MA_1692	Pabies.nocontam.gmapdb	mRNA	24487	26724	.	+	.	ID=bysp.og20_1034.all.pruned.cds.aln.fa.mrna1;Name=bysp.og20_1034.all.pruned.cds.aln.fa;Parent=bysp.og20_1034.all.pruned.cds.aln.fa.path1;coverage=95.9;identity=83.6;matches=1712;mismatches=323;indels=12;unknowns=200
MA_1692	Pabies.nocontam.gmapdb	exon	24487	26724	83	+	.	ID=bysp.og20_1034.all.pruned.cds.aln.fa.mrna1.exon1;Name=bysp.og20_1034.all.pruned.cds.aln.fa;Parent=bysp.og20_1034.all.pruned.cds.aln.fa.mrna1;Target=bysp.og20_1034.all.pruned.cds.aln.fa 89 2332 +
MA_1692	Pabies.nocontam.gmapdb	CDS	24487	24563	78	+	0	ID=bysp.og20_1034.all.pruned.cds.aln.fa.mrna1.cds1;Name=bysp.og20_1034.all.pruned.cds.aln.fa;Parent=bysp.og20_1034.all.pruned.cds.aln.fa.mrna1;Target=bysp.og20_1034.all.pruned.cds.aln.fa 89 165 +
MA_1692	Pabies.nocontam.gmapdb	CDS	24567	26724	83	+	2	ID=bysp.og20_1034.all.pruned.cds.aln.fa.mrna1.cds2;Name=bysp.og20_1034.all.pruned.cds.aln.fa;Parent=bysp.og20_1034.all.pruned.cds.aln.fa.mrna1;Target=bysp.og20_1034.all.pruned.cds.aln.fa 166 2332 +

