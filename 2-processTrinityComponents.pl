#!/usr/bin/env perl

use strict;
use warnings;
use Cwd;
use List::Util qw(sum);


########################################################################################################################################################################
# SCRIPT TO CALCULATE COVERAGE DEPTH FOR TRINITY COMPONENTS (GENES) AND SELECT THE ISOFORMS WITH THE HIGHEST PERCENT REPRESENTATION AMONG ALL ISOFORMS FOR A COMPONENT
# USES BOWTIE2, SAMTOOLS, GATK, RSEM (WHICH CALLS BOWTIE1)
# SCRIPT FOR CHOOSING HIGEST PERCENT ISO (HISO) IS MODIFIED FROM YANG AND SMITH 2012 PYTHON SCRIPT
########################################################################################################################################################################

# TAKES A TEXT FILE OF LIBRARY NAMES THAT CORRESPOND TO FINAL DIRECTORIES IN PATHS TO TRINITY.FASTA FILES AND IN MY FILE STRUCURE TO THE DIRECTORY CONATAINING /cleaned/<READS>
# RUNS THE PIPELINE FOR EACH TRINITY ASSEMBLY
# WHERE <READS> ARE FORWARD/REVERSE PAIRED READ FILES AND AN UNPAIRED READ FILE:
# MY FILE STRUCTURE FOR "CLEANED AND FILTERED" READ FASTQ FILES: /Volumes/HD2/RNAseq_working/<LIBRARY_NAME_HERE>/cleaned/
# IF FILES ARE NOT IN THIS STRUCTURE, SCRIPT SHOULD BE CHANGED ACCORDINGLY SO IT CAN FIND YOUR FILES


########################################################################################################################################################################
# SET UP DIRECTORIES
########################################################################################################################################################################


# DIRECTORY CONTAINING DIRECTORIES FOR EACH LIBRARY, WHICH CONTAIN THE READS - THE 'RNAseq_working' DIRECTORY IN THE PATH EXAMPLE ABOVE
my $seqdir = '/Volumes/HD2/RNAseq_working/';

# DIRECTORY CONTAINING DIRECTORIES FOR EACH LIBRARY, WHICH CONTAIN TRINITY.FASTA ASSEMBLIES
my $trindir = '/Users/bdorsey/Documents/Dioon/RNAseq/trinity_files/';

# SYSTEM CALL FOR GATK JAR FILE
my $gatk = 'java -jar /opt/local/NGS/GenomeAnalysisTK.jar';

# TEXT FILE CONTAINING ONE LIBRARY_NAME PER LINE - CORRESPONDS TO FINAL DIRECTORY NAMES CONTAINING 'CLEANED' FOLDER IN PATH ABOVE AND DIRECTORIES WITH TRINITY.FASTA FILES 
# my $libdir = $seqdir . 'lib_dirs2.txt';
# 
# my @libdir;
# open my $libdirFH, '<', $libdir or die "can't open $libdir";
# while( <$libdirFH> ) {
# 	chomp ( my $line = $_ );
# 	push ( @libdir, $line );
# 	}
# close( $libdirFH );

# FOR ONE INDEX ONLY
my @libdir = qw/purpusii2_RYBD010_index16/;

########################################################################################################################################################################
# START PIPELINE
########################################################################################################################################################################

foreach my $index (@libdir){
	
	my $indexdir = $trindir . $index;
	my $readdir = $seqdir . $index;
	my $logfile = $indexdir . '/' . $index . '_pipeline.log';
	
	open my $logFH, '>>', $logfile;
	
	select $logFH;
	$| = 1;
	*STDOUT = $logFH;
	*STDERR = $logFH;
	
	print "\n\n##########################################\nLOG FILE APPENDED HERE ON:\n";
	system("date >> $logfile");
	print "\n##########################################\n";
	
	
	#####################################
	# SET UP FILE NAMES
	#####################################

	my $inTime = time;
	print "\n\nStarting read mapping, coverage counts, and isoform selection pipeline for $index on:\n";
	my $inDate = system("date >> $logfile");
	print "\n\n";

	my $file1 = $readdir . '/cleaned/' . "$index" . '_cleaned_P1.fastq';
	my $compfile1 = $file1 . '.gz';
	if (!-e $file1 && -e $compfile1) {
		$file1 = &unzip($compfile1);
		}
	unless (-e  $file1) {
		die "$file1 not found";
		}
	my $file2 = $readdir . '/cleaned/' . "$index" . '_cleaned_P2.fastq';
	my $compfile2 = $file2 . '.gz';
	if (!-e $file2 && -e $compfile2) {
		$file2 = &unzip($compfile2);
		}
	unless (-e $file2){
		die "$file2 not found";
		}

	my $fileU = $readdir . '/cleaned/' . "$index" . '_cleaned_U.fastq';
	my $compfileU = $fileU . '.gz';
	if (!-e $fileU && -e $compfileU) {
		$fileU = &unzip($compfileU);
		}
	unless (-e $fileU){
		die "$fileU not found";
		}	
	
	my $file1U = $readdir . '/cleaned/' . "$index" . '_cleaned_P1U.fastq';
	my $compfile1U = $file1U . '.gz';
	if (!-e $file1U && -e $compfile1U) {
		$file1U = &unzip($compfile1U);
		}
	unless (-e $file1U){
		die "$file1U not found";
		}	

	my $assembly = $trindir . $index . '/Trinity.fasta';
	my $transcripts = $assembly;
	my $switch = $index . '_Trinity';
	$transcripts =~ s/Trinity/$switch/;
	my $trincomp = $assembly . '.bz2';

	if (!-e $assembly && -e $trincomp){
		$assembly = &bunzip($trincomp);
	}
	unless (-e $assembly or -e $transcripts){
		print "can't find Trinity file: $! \n";
		next;
	}

	########################
	# RANAME TRINITY.FASTA
	unless (-s $transcripts){
		my $rename = system ("mv $assembly $transcripts");
		if($rename){
			print "can't rename trinity file, $!\n";
			next;
		}
	}
	########################
	# FINAL OUTPUT FILE

	my $isos_picked_gt20 = $indexdir . '/' . $index . '_Hisos_gt20.fasta';
	
########################################################################################################################################################################
# BOWTIE ALIGNMENTS	AND PREPARE BAM FILES
########################################################################################################################################################################

	my $outdir = $indexdir . '/bowtie2_out/';
	
	unless(-e $outdir){
		my $mkcall=system("mkdir $outdir");
		if ($mkcall){
			print "can't make $outdir\n";
			next;
		}
	}
	
	chdir($outdir);
	
	my $unpairFail = $outdir . $index . '_unpaired_not_aligned.fastq';
	my $uncon = $outdir . $index . '_paired_notConcord_.fastq';
	my $ID = 'C18PTACXX_2';
	my $btpair = $outdir . $index . '_paired_bt.sam';
	my $btU = $outdir . $index . '_unpaired_bt.sam';

	#########################################################
	# CHECK FOR BOWTIE INDEX AND BUILD IF NOT ALREADY THERE #
	#########################################################

	# print $btfile;
	my $btlast = $outdir . $index . '.rev.1.bt2';

	if (-s $btlast){
		print "bowtie2 indices found including $btlast\n";
	
	}
	else{
		print "building bowtie index\n\n";
		my $build = system("bowtie2-build $transcripts $index > b2build.out 2>&1");
		if ($build){
			print "can't build index $!\n";
			next;
		}
	}
		
	#############################
	# START PAIRED READ MAPPING #
	#############################
	
	if (-e $btpair){
		print "\nPaired read bowtie alignment found at $btpair...\n";
	}
	else{
		print "starting paired read bowtie for $index\n";
		my $call1 = system("bowtie2 --very-sensitive -x $index -X 1000 --no-unal --un-gz $unpairFail --un-conc $uncon --phred33 --rg-id $ID --rg SM:$index --rg LB:$index -p 20 -1 $file1 -2 $file2 -S $btpair 2>> $logfile");
		if($call1){
			print "can't start bowtie paired run for $index: $!\n";
			next;
		}
	}
	
	#############################
	# START SINGLE READ MAPPING #
	#############################
	
	if (-e $btU){
		print "Single read bowtie alignment found at $btU...\n";
	}
	else{
		print "Starting single read bowtie for $index\n";
		my $call2 = system("bowtie2 --very-sensitive -x $index -U $fileU --no-unal --rg-id $ID --rg SM:$index --rg LB:$index -p 20 -S $btU 2>> $logfile");
		if($call2){
			print "can't start unpaired bowtie for $index: $!\n";
			next;
		}
	}

	#########################################################################
	# CONVERT SAM TO BAM, INDEX AND MERGE, CREATE FASTA DICTIONARY AND INDEX#
	#########################################################################
	
	########################
	# FILE NAMES
	
	my $pairbampre = $btpair;
	$pairbampre =~ s/.sam/_sort/;
	my $pairbam = $pairbampre . '.bam';
	
	my $singlebampre = $btU;
	$singlebampre =~ s/.sam/_sort/;
	my $singlebam = $singlebampre . '.bam';
	
	my $mergedbampre = $pairbampre;
	$mergedbampre =~ s/paired/merged/;
	my $mergedbam = $mergedbampre . '.bam';
	
	my $pairindex = $pairbam . '.bai';
	my $singleindex = $singlebam . '.bai';
	my $mergeindex = $mergedbam . '.bai';
	
	my $dict = $transcripts . '.dict';
	$dict =~ s/.fasta//;
	
	my $samind = $transcripts . '.fai';
	
	my $pic = 'java -jar /opt/local/NGS/picard-tools-1.99/CreateSequenceDictionary.jar';
	
	########################
	# RUN SAMTOOLS

	print "\nconverting sam files to bam\n";
	system("date >> $logfile");
	print "\n";

	if (-s ($pairbam)){
		print "$pairbam found, not recreating...\n";
	}
	else{
		print "converting paired sam to bam\n";
		&sam_to_bam_sort($btpair, $pairbampre, $index, $logfile);
	}
	
	if (-s ($singlebam)){
		print "$singlebam found, not recreating...\n";
	}
	else{
		print "converting single sam to bam\n";
		&sam_to_bam_sort($btU, $singlebampre, $index, $logfile);
	}
	
	if (-s $pairindex){
		print "$pairindex found, not recreating...\n";
	}
	else{
		print "indexing paired bam\n";
		&index_bam($pairbam, $logfile);
	}
	
	if (-s $singleindex){
		print "$singleindex found, not recreating...\n";
	}
	else{
		print "indexing single bam\n";
		&index_bam($singlebam, $logfile);
	}
	
	if (-s ($mergedbam)){
		print "$mergedbam found, not recreating...\n";
	}
	else{
		print "merging single and paired bam files\n";
		&merge_bams($singlebam, $pairbam, $mergedbam, $logfile);
	}
	
	if (-s $mergeindex){
		print "$mergeindex found, not recreating...\n";
	}
	else{
		print "indexing merged bam\n";
		&index_bam($mergedbam, $logfile);
	}
		
	if (-s $dict){
		print "$dict found, not recreating...\n";
	}
	else{
		print "\nmaking dictionary of $transcripts\n";

		my $mkdict = system("$pic R= $transcripts O= $dict >> $logfile 2>&1");
		if ($mkdict){
			print "can't make $transcripts dictionary, $!\n\n";
			next;
		}
	}
	
	if (-e $samind){
		print "$samind found, not recreating...\n";
	}
	else{
		print "indexing $transcripts\n";
		
		my $indcall = system("samtools faidx $transcripts >> $logfile 2>&1");
		if ($indcall){
		die "indexing $transcripts failed, ";
		
		}
	}
	
	unless (-s $mergeindex and -s $samind and -s $dict and -s $mergedbam and -s $singleindex and -s $pairindex and -s $singlebam and -s $pairbam) {
		die "samtools files not complete ";
	}
	
	
	print "\nFinished mapping and converting $index at:\n";
	my $outdate = system("date >> $logfile");
	print "\n";

########################################################################################################################################################################
# GET COVERAGE FOR ALL COMPS
########################################################################################################################################################################

	######################################
	# RUN GATK AND SAMTOOLS FOR COVERAGE #
	######################################


	my $gatdir = $indexdir . '/coverage';
	
	unless (-e $gatdir){
		my $mkgat = system("mkdir $gatdir");
		if ($mkgat){
			die "can't make $gatdir";
		}
	}
	
	chdir($gatdir);
	
	my $coveragefile = $gatdir . '/' . $index . '_cover';
	
	if (-s $coveragefile){
		print "found $coveragefile\nNot running GATK depthofCoverage tool\n";
	}
	else{
		print "\nrunning GATK depthofCoverage tool for $index\n";
		
		my $gatcall = system ("$gatk -T DepthOfCoverage -R $transcripts -I $mergedbam -pt library -o $coveragefile >> $logfile 2>&1");
	
		if ($gatcall){
			print "GATK failed for $index, $!";
			next;
		}
	}	
	
	my $samout = $gatdir . '/coverage_stats.txt';
	
	if (-s $samout){
		print "$samout found\nNot running samtools idxstats tool for $transcripts\n";
	}
	else{
		print "\nrunning samtools idxstats tool on $transcripts\n";
		
		my $sam = system ("samtools idxstats $mergedbam > $samout 2>> $logfile");
		if ($sam) {
			die "samtools failed at ";
		}
	}
	
	unless (-s $coveragefile && -s $samout){
		die "Coverage files from GATK and samtools not found, quitting ";
	}
	
########################################################################################################################################################################	
# PICK ISOFORMS
########################################################################################################################################################################

	######################################
	# RUN RSEM FOR ISOFORM PERCENTAGE    #
	######################################
	# USING ONLY THE SINGLE AND LEFT READS AS THE PAIRS ARE NOT INDEPENDENT
	
	my $rsem_dir = $indexdir . '/rsem';
	unless (-e $rsem_dir){
		my $mkrsem = system("mkdir $rsem_dir");
		if ($mkrsem){
			die "can't make $rsem_dir";
		}
	}

	my $fasta = $transcripts;
	my $rsemfile = $rsem_dir . '/RSEM.isoforms.results';
	my $iso_outfile = $rsem_dir . '/' . $index . '_Hiso.fasta';
	my $iso_list = $rsem_dir . '/' . $index . '_isolist.txt';
	
	
		
	if (-s $rsemfile){
			print "$rsemfile file found, not running RSEM\n";
	}
	else{
		print "\nstarting RSEM isoform percentage estimation for $index\n";
	
		chdir($rsem_dir);
		my $rsem_call = system("/opt/local/NGS/trinityrnaseq_r2013_08_14/util/RSEM_util/run_RSEM_align_n_estimate.pl --transcripts $transcripts --seqType fq --single $file1U --thread_count 20 > rsem.log 2>&1");

		if ($rsem_call){
			die "RSEM failed at ";
		}
	}

		
	########################################################
	# PICK REPRESENTATIVE ISOFORMS WITH HIGHEST PERCENTAGE #
	########################################################
	
	my $numberOfIsos;
	my $number_listlines;
	my @Hiso_picked;

	
	if ( -s $iso_outfile && -s $iso_list ){
	
		print "Found $iso_outfile\nSkipping isoform filtering step\n\n";
		print "Isolist file found, creating Hiso_picked array from this file\n";
	
		open my $isolistFH, '<', $iso_list;
		while (<$isolistFH>){
			chomp;
			push (@Hiso_picked, $_);
		}
	}
	else {	
	
		print "Picking representative isoforms using RSEM output\n";
	
		open my $rsemFH, '<', $rsemfile or die "cannot open rsem file ";

		my $last_gene = "";
		my @spls;
		my $spls;
		my $gene;
		my @block;

		###########################################################################################
		# READ RSEM FILE, GROUP ISOFORMS, AND PICK THE ONE WITH HIGHEST PERCENTAGE REPRESENTATION #

		my $rsemhead = <$rsemFH>;									# SKIP HEADER LINE
		while (<$rsemFH>) { 										# READ IN LINES OF RSEM FILE
		# 	chomp(my $line = $_);
			@spls = split;										# CONVERT TABBED LINE INTO ARRAY(/\t/, $line)
		# 	$spls = \@spls;										# MAKE A REFERENCE TO THE ARRAY
			$gene = $spls[1];									# ASSIGN THE CURRENT COMPONENT TO VAR GENE

			if (($gene eq $last_gene) or ((scalar @block) == 0)){		# CHECK IF SEQ IS PART OF CURRNT COMPONENT GROUP OR IF FIRST ISO
				push @block, [@spls];							# IF ISO IS PART OF CURRENT COMPONENT ADD REF TO LINE TO BLOCK ARRAY
				$last_gene = $gene;

			}else {												# IF NOT RUN PICK_ISO SUB ON CURRENT BLOCK OF ISOS
		# 		print "picking iso\n";
				push (@Hiso_picked, &pick_iso(\@block));
				@block = [@spls];								# REPLACE BLOCK ARRAY WITH LINE FOR NEXT COMPONENT - THIS IS NOT WORKING CORRECTLY!
				$last_gene = $gene;
			}


			# JUST FOR DEBUGGING - CHECKING BLOCK
		# 	foreach my $row (@block){
		# # 	print "$gene\n";
		# 		foreach my $column (@{$row}){					
		# 			print "$column\t"; 
		# 		}
		# 		print "\n";
		# 	}
		}
		push (@Hiso_picked, &pick_iso(\@block));				# PROCESS THE LAST BLOCK
	
		$numberOfIsos = @Hiso_picked;
		print "number of isos in Hiso_picked array: $numberOfIsos\n";
	
		close ($rsemFH);
			
		##################################
		# PRINT LIST OF ISOFORM SELECTED #

		open my $iso_listFH, '>', $iso_list;

		foreach my $iso ( @Hiso_picked ){
			print $iso_listFH "$iso\n";
		}
		close $iso_listFH;
		
		open $iso_listFH, '<', $iso_list;

		my @listlines;
		while ( <$iso_listFH> ){
			push ( @listlines, $_ );
		}	
		$number_listlines = @listlines;

		close $iso_listFH;

		print "iso list contains $number_listlines lines\n";
	
		close $iso_listFH;
	
		
		#########################################
		# PRINT ISOFORMS SELECTED TO FASTA FILE #
		
		my $ret = &fetch_seqs( $fasta, \@Hiso_picked, $iso_outfile );
		
		my @seqs_fetched;
		
		open my $iso_outfileFH, '<', $iso_outfile;
	
		while ( <$iso_outfileFH> ){
			push ( @seqs_fetched, $_ );
		}
		close $iso_outfileFH;
			
		my $num_seqs_fetched = grep ( />/, @seqs_fetched );
		print "Sequences fetched: $num_seqs_fetched\n";
		
		unless ( ($numberOfIsos == $number_listlines) && ($numberOfIsos == $num_seqs_fetched) && ($number_listlines == $num_seqs_fetched) ) {
			die "isoform output files do not match ";
		}
				
		if ( -s $iso_outfile ){
			print "isoforms written to: $iso_outfile\n";
		}
	

	}

########################################################################################################################################################################
# GET COVERAGE FOR ISOFORMS PICKED AND WRITE ISOS WITH GT 20X TO FILE
########################################################################################################################################################################


	####################################################################################################
	# REDUCE TRANSCRIPTS TO ONLY THOSE ISOFORMS PICKED AND WHICH HAVE AT LEAST 20 READS MAPPED TO THEM #
	# CREATE A HASH OF TRANSCRIPTS WITH AT LEAST 20 READS MAPPED TO STORE COVERAGE DATA                #
	####################################################################################################
	
	
# 	NOT RUNNING THE REDUCTION STEP, JUST MAKING A HASH OF ALL ISOS PICKED WITH VALUES = EMPTY ARRAYS
# 	my %dephash = map { $_ => [] } @Hiso_picked;
	
# 	open my $statFH, '<', $samout or die "cannot open $samout ";
# 	open my $iso_listFH, '<', $iso_list or die "cannot open $iso_list ";
# 	my @isoarray;
# 	
# 	while (<$iso_listFH>){
# 	push (@isoarray, $_);
# 	}
# 	close $iso_listFH;
# 	
# 	my %isohash = map { $_ => 1 } @isoarray;
# 		
# 	while (<$statFH>){
# 		my ($comp, $reads) = ( split /\t/ )[0,2];
# 		
# 		if ( ($reads  >= 20) && ($isohash{$comp}) ){  #grep /$comp/, @isoarray
# 			$dephash{$comp}=[];				
# 		}
# 	}		
# 	close $statFH;
			
# 		chomp(my $line = $_);
# 		if ($line =~ /(\w+)\t(\d+)\t(\d+)\t\d+/){
# 			my $comp = $1;
# 			my $reads = $3;
# 		}
	
	

# FOR DEBUGGING ONLY	
# 	my @tarray = keys %dephash;
# 	foreach my $c (@tarray[0..5]){
# 		print "From reducing loop:\n$c\n";
# 	}

	
# FOR DEBUGGING ONLY
# 	print "\nfrom sum loop\n";
# 	
# 	print "First position of comp40793_c1_seq1: $dephash{'comp40793_c1_seq1'}[0]\n";
# 	foreach my $c (@tarray[0..5]){
# 		my $S = sum @{ $dephash{$c} };
# 		my @a1 = @{ $dephash{$c} };
# 		my $S2 = sum @a1;
# 		print "\n$c: $S, $S2\n";
# 		foreach my $pos (@a1[0..5]){
# 			print "$pos\n";
# 		}
# 		
# 	my @varray = values %dephash;
# 	foreach my $v (@varray[0..5]){
# 		print "$v\n";
# 	}
# 
# 	}
# 	die;

	####################################################################################################
	#CHECK FOR FILES OF COVERAGE BY POSITION AND MEAN COVERAGE										   #
	# IF ALREADY EXIST, MAKE A HASH OF ISOS AND THEIR MEAN COVERAGE FROM THE MEAN COVERAGE FILE																				                   #
	####################################################################################################


	print "Getting coverage for isoforms picked\n";


	my $posfile = $gatdir . '/Hiso_coverage_byPos.txt';
	my $meanfile = $gatdir . '/Hiso_meanCoverage.txt';
	my %meanhash;
	my @posfile_lines;
	my @meanfile_lines;
	my %dephash;
	my %meanlinehash;
	my %poslinehash;

	if ( -s $posfile && -s $meanfile ){
		print "$posfile and $meanfile found, not recreating\nCreating meanhash for later\n";
		
		open my $posfileFH, '<', $posfile or die "cannot open $posfile ";
		open my $meanfileFH, '<', $meanfile or die "cannot open $meanfile ";
		
		while ( <$posfileFH> ){
			$poslinehash{$_} =1;
		}
		my $num_pos_lines = ( keys %poslinehash );
		
		while ( <$meanfileFH> ) {
			my ( $iso, $mean ) = ( split /,/ )[0,1];
			$meanhash{$iso} = $mean;
			
			$meanlinehash{$_} = 1;
		}
		my $num_mean_lines = ( keys %meanlinehash );
		
		print "Number of lines in By-position file: $num_pos_lines\n";
		print "Number of lines in mean depth file: $num_mean_lines\n";
			
		close $posfileFH;
		close $meanfileFH;

		my $mhashlength = ( keys %meanhash );
		
		unless ( ($num_pos_lines == $num_mean_lines) && ($num_pos_lines == $mhashlength) ){
			die "Position file and mean file are not the same length ";
		}
		
		system ( "date >> $logfile" );
		print "\n";
	}
	else{
	
		##########################################################################################################
		# IF FILES NOT FOUND:
		# CREATE ARRAYS OF COVERAGE BY POSITION FOR EACH COMP AND ASSIGN THEM TO THE KEYS (COMPS) IN THE DEPHASH #
		##########################################################################################################
	

		open my $coveragefileFH, '<', $coveragefile or die "cannot open $coveragefile ";
		while ( <$coveragefileFH> ){
			chomp ( my $line = $_ );
			if ( $line =~ /^(\w+):(\d+)\t(\d+)/ ) {
		
				my $comp = $1;
				my $pos = $2;
				my $dep = $3;
				my $apos = $pos - 1;
				$dephash{$comp}[$apos] = $dep;
			
			}
	
		}
		close $coveragefileFH;


		##############################################################
		# PRINT OUT CONTENTS OF EACH BY-POSITION ARRAY IN A CSV FILE #
		# PRINT OUT MEAN COVERAGE FOR EACH COMP IN A CSV FILE #
		##############################################################

		open my $posfileFH, '>', $posfile or die "cannot open $posfile ";
		open my $meanfileFH, '>', $meanfile or die "cannot open $meanfile ";
		
		foreach my $locus ( keys %dephash ){
			my @ar = @{ $dephash{$locus} };
			my $l = @ar;
			my $sum = sum @ar;
			my $mean = $sum/$l;

			print $posfileFH "$locus";

			foreach my $pos ( @ar ){
				print $posfileFH ",$pos";
			}
			print $posfileFH "\n";

			$meanhash{$locus} = $mean;
			print $meanfileFH "$locus, $mean\n";
	
		}
		close $posfileFH;
		close $meanfileFH;
		
		open $posfileFH, '<', $posfile or die "cannot open $posfile ";
		open $meanfileFH, '<', $meanfile or die "cannot open $meanfile ";
				
		while ( <$posfileFH> ){
			$poslinehash{$_} = 1;
		}
		my $num_pos_lines = ( keys %poslinehash );
		
		while ( <$meanfileFH> ){
			$meanlinehash{$_} = 1;
		}
		my $num_mean_lines = ( keys %meanlinehash );
		
		print "Number of lines in By-position file: $num_pos_lines\n";
		print "Number of lines in mean depth file: $num_mean_lines\n";
			
		close $posfileFH;
		close $meanfileFH;
		
		if ( -s $posfile ){
			print "isoform coverage by position written to $posfile\n";
		}
		else{ die "$posfile not created or found ";}
		
		if ( -s $meanfile ){

			print "mean isoform coverage printed to $meanfile\n";
		}
		else {die "$meanfile not created or found "};
		
		unless ( $num_pos_lines == $num_mean_lines ){
			die "Position file and mean file are not the same length ";
		}
		
		system( "date >> $logfile" );
		print "\n\n";
	}
	

	##################################################################
	# READ MEANHASH AND SELECT ISO WITH >20X COVERAGE, PRINT OUT LIST
	##################################################################

	my $gt20list = $rsem_dir . '/iso_gt20_list.txt';
	my @isos20;
	my $length_isos20;
	my %gt20list_lines;
	
	if ( -s $gt20list ){
		print "$gt20list found, not replacing...\nCreating iso array for checking fasta file later\n";

		foreach my $iso (keys %meanhash){

			if ( $meanhash{$iso} >= 20 ){
				push ( @isos20, $iso );
			}
		}
		
		open my $gt20listFH, '<', $gt20list;
		while ( <$gt20listFH> ){
			$gt20list_lines{$_} = 1;
		}
		my $num_gt20list_lines = ( keys %gt20list_lines );
		
		$length_isos20 = @isos20;
		
		unless ( $num_gt20list_lines == $length_isos20 ){
			die "Isos20 list array and list file not the same length ";
		}	
		
		close $gt20listFH;
		
		system("date >> $logfile");
		print "\n\n";
	}
	else{
		open my $gt20listFH, '>', $gt20list or die "cannot open $gt20list ";

		foreach my $iso ( keys %meanhash ){
			if ( $meanhash{$iso} >= 20 ){
				my $m = $meanhash{$iso};
				push (@isos20, $iso);
				print $gt20listFH "$iso, $m\n";
			}
		}
		close $gt20listFH;
		open $gt20listFH, '<', $gt20list or die "cannot open $gt20list ";
		
		while (<$gt20listFH>){
			$gt20list_lines{$_} = 1;
		}
		my $num_gt20list_lines = ( keys %gt20list_lines );
		
		$length_isos20 = @isos20;
		
		unless ($num_gt20list_lines == $length_isos20){
			die "Isos20 list array and list file not the same length ";
		}
		
		close $gt20listFH;
		
		print "list of isoforms with gt 20x coverage written to $gt20list\n";
		system("date >> $logfile");
		print "\n\n";
	}

	#############################################
	# FILTER BY-POSITION FILE ACCORDING TO >20X
	#############################################
	
	my $pos_two = $posfile;
	$pos_two =~ s/.txt/_gt20.txt/;
	my %pos_two_lines;
	
	if ( -s $pos_two ){
		print "$pos_two file found, not replacing...\n";
		
		open my $pos_twoFH, '<', $pos_two;
		
		while (<$pos_twoFH>){
			$pos_two_lines{$_} = 1;
		}
		my $num_pos_two_lines = ( keys %pos_two_lines );
		
		print "Filtered position file contains $num_pos_two_lines lines\nIsos20 array has $length_isos20 elements\n";
		
		unless ($num_pos_two_lines == $length_isos20){
			die "Lenght of isos20 array and filtered by-position file don't match ";
		}
		
		close $pos_twoFH;
		
		
		system("date >> $logfile");
		print "\n\n";

	}
	else{
	
		print "Filtering by-position file to only seq with >20x coverage\n";
		open my $pos_twoFH, '>', $pos_two or die "cannot open $pos_two ";
		open my $posfileFH, '<', $posfile or die "cannot open $posfile ";


		my %isos20hash = map { $_ => 1 } @isos20;

		while (<$posfileFH>) {
			my $comp = ( split /,/ )[0];
			print $pos_twoFH $_ if $isos20hash{$comp};
			
		}
		close $posfileFH;
		close $pos_twoFH;
		
		open $pos_twoFH, '<', $pos_two or die "cannot open $pos_two ";
		
		while (<$pos_twoFH>){
			$pos_two_lines{$_} = 1;
		}
		my $num_pos_two_lines = ( keys %pos_two_lines );
		
		print "Filtered position file contains $num_pos_two_lines lines\nIsos20 array has $length_isos20 elements\n";
		
		unless ($num_pos_two_lines == $length_isos20){
			die "Lenght of isos20 array and filtered by-position file don't match ";
		}
		
		close $pos_twoFH;
		print "\nBy position coverage for isoforms with gt 20x coverage written to $pos_two\n";
		system("date >> $logfile");
		print "\n\n";

	}


	######################################################################
	# PRINT ISOS CHOSEN THAT HAVE GT 20X COVERAGE TO FASTA FILE #
	######################################################################


	################################################################
# 	ONLY RETRIEVING 7166 SEQS WHITH ISO20 ARRAY OF 10029 ELEMENTS WHEN FILTERING ISO_OUTFILE
# 	FILTERING ORIGINAL TRANSCRIPTS FILE GIVES CORRECT 10029, BUT LINE CHECKING CODE IS NOT WORKING CORRECTLY
# 	THE 7166 IS THE INTERSECTION OF ISOS WITH HIGHEST PERCENTAGE AND THOSE WITH GT 20X COVERAGE, THIS IS THE CORRECT RESULT
# 	################################################################
	
	my %finalfasta;
		
	if (-s $isos_picked_gt20){
	
		print "$isos_picked_gt20 file found, not replacing...\n";
		system("date >> $logfile");
		print "\n\n";
	
	}
	else {
		
		print "Writing chosen isoform sequences to fasta file\n\n";
		
		my $code = &fetch_seqs($iso_outfile, \@isos20, $isos_picked_gt20);
	}

	open my $finalFH, '<', $isos_picked_gt20;
	while (<$finalFH>){
		chomp;
		if ($_ =~ /^>comp/){
		$finalfasta{$_} = 1;
		}
	}
	close $finalFH;
	
	my $fastalength = ( keys %finalfasta );
	print "\nNumber of seqs written to final fasta file: $fastalength\n";
	
	my %Hisohash = map { $_ => 1 } @Hiso_picked;
	my %iso20hash = map { $_ => 1 } @isos20;
	my $intersection = 0;
	
	foreach my $Hiso ( keys %Hisohash ) {
		$intersection++ if $iso20hash{$Hiso};
	}
	
	print "\nintersection: $intersection\n";
	
	if ( -s $isos_picked_gt20 && ($intersection == $fastalength) ){
	
		print "Isoforms with coverage gt 20x written to $isos_picked_gt20\n";
		print "Sequence number in file matches intersection of Hiso list and isos gt 20X\n";
		print"\nFinished Isoform coverage pipeline for $index at\n";
		system("date >> $logfile");
		print"\n";
	
	}
	else{
		die "$isos_picked_gt20 not created or does not match length of iso array ";
	}

	print "\ncompressing files\n";
	
# 	&compress($file1);
# 	&compress($file2);
# 	&compress($fileU);
# 	&compress($file1U);
#  	&compress($transcripts);
#  	&compress($dict);
#  	&compress($samind);
 	&compress($indexdir);
 	&compress($readdir);
 	
 	print "Moving to next library\n\n";
 	
 	close $logFH;

}


########################################################################################################################################################################
# SUBS
########################################################################################################################################################################


sub compress {
    my ($file) = @_;
    my $call = system("gzip -r $file &");
}

sub unzip  {
    my ($file) = @_;
    my $call = system("gunzip $file");
    my $new = $file;
    $new =~ s/\.gz//;
    return($new);
}

sub bunzip  {
    my ($file) = @_;
    my $call = system("bunzip2 $file");
    my $new = $file;
    $new =~ s/\.bz2//;
    return($new);
}

sub sam_to_bam_sort {
	my ($samfile, $bamfilepre, $index, $logfile) = @_;
	my $vcall = "samtools view -bSuh $samfile";
	my $scall = "samtools sort - $bamfilepre";
	my $pcall = join('|', $vcall, $scall);
	my $call = system("$pcall >> $logfile 2>&1 ");
	if ($call){
		return("sam_to_bam failed for $index");
	}

}

sub index_bam {
	my ($bamfile, $logfile) = @_;
	my $call = system("samtools index $bamfile >> $logfile 2>&1");
	if ($call){
		return("failed indexing $bamfile");
	}
}


sub merge_bams{
	my ($pbam, $sbam, $out, $logfile) = @_;
# 	foreach ($pbam, $sbam, $out){$_ .= '.bam'};
	my $call = system("samtools merge $out $pbam $sbam >> $logfile 2>&1");
	if ($call){
		return("merging bams failed: $!\n");
	}
}


###########################################################################
# CHOOSE ISOFORM WITH HIGHEST PERCENT REPRESENTATION PER COMPONENT (GENE) #
###########################################################################

sub pick_iso {
	my $block = shift;
# 	my @block = @{$refblock};
# 	my $block = \@block;
	if (scalar(@{$block}) == 1) {return $block->[0][0]}
	else {
		my ($transcript, $length, $IsoPct) = ($block->[0][0], $block->[0][2], $block->[0][7]);
		foreach my $i (1 .. $#{$block}) {
			my $perc = $block->[$i][7];
			if ( ($perc > $IsoPct) or ($perc == $IsoPct && $block->[$i][2] > $length) ) {
				($transcript, $length, $IsoPct) = ($block->[$i][0], $block->[$i][2], $perc);
			}
		}
		return $transcript;
	}
}

###########################################
# SELECT SEQS FROM INTERLEAVED FASTA FILE #
###########################################
sub fetch_seqs {
	
	my ($fasta, $Hiso_picked, $outfile) = @_;
	my @Hiso_picked = @{$Hiso_picked};
	my %Hisohash = map { $_ => 1 } @Hiso_picked;
	my $HisohashLength = ( keys %Hisohash );
	print "Hisohash has $HisohashLength elements\n";
	
	open my $fastaFH, '<', $fasta;
	open my $outfileFH, '>', $outfile;
	
	{	
		local $/ = '>';
		my $first = <$fastaFH>;
		while (<$fastaFH>){
			chomp;
			my $iso = ( split /\s/)[0];
				print $outfileFH ">$_" if $Hisohash{$iso};
		}
	
	}
	
	close $fastaFH;
	close $outfileFH;
	return (1);
}

# CHANGED THE PATTERN MATCH TO HASH LOOKUP
# sub fetch_seqs {
# 	my ($fasta, $Hiso_picked, $outfile) = @_;
# 	my @Hiso_picked = @{$Hiso_picked};
# 	open my $fastaFH, '<', $fasta;
# 	open my $outfileFH, '>', $outfile;
# # 	my $pat = join ("|", @Hiso_picked);
# 	my %Hisohash = map { $_ => 1 } @Hiso_picked;
# 
# # 	{
# # 		local $/ = '>';
# # 		while (my $line = <$fastaFH>){
# # 			if ($line =~ /($pat)/){
# # 				$line =~ s/>//;
# # 				print $outfileFH '>'.$line;
# # 			}	
# # 		}
# # 	
# # 	}
# 	{	
# 		local $/ = '>';
# 		while (<$fastaFH>){
# 			chomp;
# 			my $iso = ( split )[0];
# 			print $outfileFH '>'.$_ if $Hisohash{$iso};
# 		}
# 	
# 	}
# 	
# 	close $fastaFH;
# 	close $outfileFH;
# }
