use strict;
use warnings;

######################################################################################
# This is a script to process raw RNAseq reads. It is based on a script Sonal Singal #
# used in Bi et al. 2012. Please attribute to original authors if you use this or    #
# any modification thereof as I have only rearranged the code to fit my needs and    #
# have not added any substantial functionality.										 #
######################################################################################


##################################################################################
# a script to clean up reads (adaptor/duplicate/contamination removal), trimming #
# external dependencies: flash, trimmomatic, bowtie2, cutadapt, cope             #         
# written by Sonal Singhal, sonal.singhal1 [at] gmail.com, 28 Dec 2012           #
# version 1.5 -- made a little faster, still                                     #
##################################################################################


# my pipeline:
# 1. trimmomatic PE and cutadapt recursively (3x?)
# hard trim if kmers or GC content are odd at ends - like in adapter 2
# 2. remove duplicates ? NO
# 3. fixMatepairs 
# 4. mergeReads yes
# 5. removeContaminants

#PARENT DIRECTORY WITH ALL LIB DIRECTORIES
my $homedir = '/Users/bdorsey/Documents/Dioon/RNAseq/Raw.read.data/mejiae_script_test/';

#THIS FILE HAS A LIST OF DIRECTORIES, ONE FOR EACH LIBRARY
my $libdir = $homedir . 'lib_dirs.txt';

#THIS FILE HAS INFORMATION ABOUT THE LIBRARIES, HERE THE FIRST COLUMN HAS THE LIBRARY NAMES THE SECOND COLUMN HAS THE ADAPTOR AND THIRD HAS THE SEQ
my $libAdapt = $homedir . 'lib_adaptors.txt';

#CONTAMINANT GENOME IN FASTA FORMAT; HERE I AM USING HUMAN
my $contam = '/opt/local/NGS/bowtie-1.0.0/genomes/Human.rna.fa' ; #need to make this an index

#THE UNIVERSAL ADAPTOR SEQUENCE
my $uniad = 'AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT';

#PATHS TO DIFFERENT EXECUTABLES
my $trimmomatic = '/opt/local/NGS/bin/trimmomatic-0.30.jar';
my $cutadapt = '/opt/local/NGS/cutadapt-1.2.1/bin/cutadapt';
my $flash = '/opt/local/NGS/FLASH-1.2.6/flash';
my $bowtie = '/opt/local/NGS/bin/bowtie2-2.1.0/bowtie2';
#my $flashslow =  'flash_slower';
#my $cope = 'cope';

# MINIMUN ACCEPTABLE LENGTH FOR READS
my $minLength = 72;


# I WANT TO OPEN EACH DIRECTORY TO GAIN ACCESS TO READ FILES
# NEED TO HAVE A TEXT FILE WITH THE NAMES OF THE DIRECTORIES
# THEN HAVE SCRIPT OPEN THE DIRECTORY AND PULL FILES FROM THERE.
# NEED TO GET THE ADAPTORS FOR EACH LIBRARY


# PARSE LIBRARY DIR FILE INTO AN ARRAY
my @libdir;
open(IN, "<$libdir");
while(<IN>) {
	chomp(my $line = $_);
	push(@libdir, $line);
	}
close(IN);	

# MAKE A HASH OF DIRECTORY/LIB AND THEIR ADAPTOR SEQS
my %adhash;
open(IN, "<$libAdapt");
	while (<IN>) {
	chomp (my $line = $_);
	#SPLIT TABBED LINE INTO AN ARRAY
	my @list = split(/\t/,$line);
	#CREATE A HASH WITH KEYS = LIBRARY AND VALUES = A HASH WITH KEYS = NAME AND SEQ AND VALUES = NAME AND SEQS
	$adhash{$list[0]} = {'name' => $list[1], 'seq' => $list[2]};
	}
close (IN);	


# OPEN EACH DIRECTORY AND CREATE A DIRECTORY FOR OUTPUT
foreach my $dir (@libdir){
	my $file1 = $homedir . $dir . '/filtered/' . $dir . '_R1_filt.fastq';# if (-e  $file1) else die "$file1 not found, stopped";
	my $file2 = $homedir . $dir . '/filtered/' . $dir . '_R2_filt.fastq';# if (-e $file2) else die "$file2 not found, stopped";
	my $outdir = $homedir .$dir . '/filtered/';
	
	# CREATE A REF TO A HASH INDEXING THE ADAPTERS FOR THIS LIBRARY
	my $ad = {"uni" => $uniad, "uni_rc" => rc($uniad), "index" => $adhash{$dir}{'seq'}, "index_rc" => rc($adhash{$dir}{'seq'})};


#####################
# FIRST TRIMMING LOOP

	print "Starting the first adaptor trimming round... \n\n";
	my $starttrim = time;
	
	# ARRAY TO HOLD FILE NAMES OF TRIMMOMATIC OUTPUT
	my @cleaned1;
	
#######################
# TRIMMOMATIC
# INPUT IS ORIGINAL READ FILES (FILTERED)
# OUTPUT IS FOUR FILES: 2 OF PAIRED READS ('FILE*_TRIM_PAIRED*_OUT) AND 2 OF UNPAIRED READS (FILE*_UNPAIRED*_OUT)
	
	my $PEout = trimmomaticPE($file1, $file2, $ad);
	
	# JUST A FILE TO CHECK THE ORDER OF OUTPUT FILES - CAN COMMENT OUT
	open(OUT, ">$outdir . clean1_out");
		
		#FOR EACH OUTPUT FILE FROM TRIMMOMATIC, RUN CUTADAPT LOOKING FOR EACH ADAPTOR
		for my $trimoutfile (sort keys %{$PEout}){

			# INPUT IS EACH FILE FROM TIMMOMATIC - OUTPUT IS FOUR FILES WITH '_trim2' APPENDED
			my $cut1 = cutadapt($ad, $PEout->{$trimoutfile}, $PEout->{$trimoutfile}, 'trim2');

			# INPUT IS EACH FILE FROM CUTADAPT - OUTPUT IS EACH FILE WITH 'trim3' APPENDED
			my $cut2 = cutadapt($ad,$cut1,$PEout->{$trimoutfile},'trim3');

			# CREATE FILE NAME FOR OUPUT OF TRIMMING ROUND = FILENAME WITH '_cleaned1' APPENDED
			my $final = $PEout->{$trimoutfile} . '_cleaned1';
			
			# RENAME OUPUT OF CUT2 AS FINAL
			my $call = system("mv $cut2 $final");
			
			#ADD FILE NAME OF CLEANED1 FILES TO ARRAY @cleaned1
			push(@cleaned1, $final); 
			}
		# PRINT OUT FILE WITH FILE NAMES - JUST FOR CHECKING ORDER	
		for my $el (@cleaned1){	
		print OUT "$el\n";
		}
	close OUT;
	my $endtrim = time;
	my $time1 = int((time - $starttrim)/60);
	print "First trim round done in $time1 minutes! \nOn to merging read pairs \n\n";



####################
# FIX READ PAIRS
	print "Fixing read pairs...\n\n";
	my $mergestart = time;
	my $newpairs = fixMatePair(\@cleaned1,'fixed1');

# ARRAY OF READ FILENAMES FROM TRIM THAT HAVE MATES
	my @trimpairs = ($newpairs->{'1'},$newpairs->{'2'});
	my $trimpairs = \@trimpairs;
	
	
	# CHECKING OUTPUT
	open (OUTP, ">$outdir . trimp");
	for my $element (@trimpairs){
		print OUTP "$element\n";
		}
	close (OUTP);	

###################
# MERGE READ PAIRS
	# MERGE READS FROM TWO FILES - input is directory for renaming files, paired read files from trimming and fixMatePair, and base for renaming
	# OUTPUT IS NEW P1 AND P2 FILES (not in testing:REPLACING FILES FROM FIXMATEPAIR) AND AN M FILE (MERGED READS), AND A REF TO A HASH %newreads HOLDING THE FILNAMES
	print "Merging read pairs...\n\n";
	my $mergeout1 = mergeReads($trimpairs,"merged1"); 
	
	my $time2 = int((time-$mergestart)/60);
	print "\n\nFixed and merged read pairs in $time2 minutes... lets trim again! \n\n";



######################
# SECOND TRIMMING LOOP
# 
	print "Starting the second adaptor trimming round... \n\n";
	my $starttrim2 = time;
	my @clean2;

	# PE TRIMMING ON UNMERGED PAIRED READS FROM ABOVE
	my $PEout2 = trimmomaticPE($mergeout1->{'1'}, $mergeout1->{'2'}, $ad);
	
	# JUST A FILE TO CHECK THE ORDER OF OUTPUT FILES - CAN COMMENT OUT
	open(OUT, ">$outdir . clean2_out");

		
	#CAT UNPAIRED AND MERGED READS FROM ROUND 1 AND PEOUT2
	my $un = $outdir . $dir . '_for_cut2_u';	
	my $call1 = system("cat $mergout1->{'m'} $newpairs->{'u'} $PEout2->{'u1'} $PEout2->{'u2'} > $un");	

	my $forcut2 = {'1' => $PEout2->{'p1'}, '2' => $PEout2->{'p2'}, '3' => $un};	
		
	#FOR EACH OUTPUT FILE IN $forcut2, RUN CUTADAPT LOOKING FOR EACH ADAPTOR
	for my $k (sort keys %{$forcut2}) {	
		my $cut1 = cutadapt($ad, $forcut2->{$k}, $forcut2->{$k}, 'trim2');			
		my $cut2 = cutadapt($ad, $cut1, $forcut2{$k}, 'trim3');	
		my $final = $forcut2->{$k} . '_cleaned2';	
		my $call = system("mv $cut2 $final"); 		#output trimmed reads as file*.cleaned2	
		unlink($cut1, $cut2);	
		push(@cleaned2, $final); 					#add file name of cleaned2 files to array @cleaned2
		}	
	
	# PRINT OUT FILE WITH FILE NAMES - JUST FOR CHECKING ORDER	
	for my $el (@cleaned2){	
		print OUT "$el\n";
		}
	close (OUT);	
	my $time3 = int((time - $starttrim2)/60);
	print "Second trim round done in $time3 minutes! \nOn to merging read pairs... \n\n";

	
####################
# FIX READ PAIRS
	print "Fixing read pairs...\n\n";
	my $mergestart = time;
	my $newpairs2 = fixMatePair(\@cleaned2,'fixed2');

# ARRAY OF READ FILENAMES FROM TRIM THAT HAVE MATES
	my @trimpairs2 = ($newpairs2->{'1'},$newpairs2->{'2'});
	my $trimpairs2 = \@trimpairs2;
	
	
	# CHECKING OUTPUT
	open (OUTP, ">$outdir . trimp2");
	for my $element (@trimpairs2){
		print OUTP "$element\n";
		}
	close (OUTP);	
	
###################
# MERGE READ PAIRS
	my $mergestart = time;
	my $mergeout2 = mergeReads($trimpairs2,"merged2"); # merge reads from two files
	my $time2 = int((time-$mergestart)/60);
	print "Merged read pairs in $time2 minutes... time to clean up the mess! \n";


###################
# MAKE FINAL FILES
	my $finalUn = $outdir . 'cleaned/' . $dir . '_final_U';
	my $finalP1 = $outdir . 'cleaned/' . $dir . '_final_P1';
	my $finalP2 = $outdir . 'cleaned/' . $dir . '_final_P2';
	my $callU = system("cat $newpairs2->{'u'} $mergeout2->{'m'} > $finalUn");
	my $callp1 = system("mv $mergeout2->{'1'} $finalP1");
	my $callp2 = system("mv $mergeout2->{'2'} $finalP2");
# 	compress($finalUn);
# 	compress($finalP1);
# 	compress($finalP2);
# 	compress($file1);
# 	compress($file2);
# 
# 
# 	my $call_rm1 = system("rm $outdir*clean*");
# 	my $call_rm2 = system("rm $outdir*trim*");
# 	my $call_rm2 = system("rm $outdir*merge*");




####################
# Subroutines      #
####################

sub compress {
    my ($file) = @_;
    my $call = system("gzip $file");
}


	
sub trimmomaticPE {
	# need input1 and input2, pairedout1, pairedout2, unpaired1, unpaired2
	my ($file1, $file2, $ad) = @_;
	#my @readfiles = ($file1, $file2); #not needed
	my %trimoutfiles = ('p1' => $file1 . '_trim_paired1_out', 'u1' => $file1 .'_trim_unpaired1_out', 'p2' => $file2 . '_trim_paired2_out', 'u2' => $file2 . '_trim_unpaired2_out');
	my $adfile = $outdir . 'adfile.fa';
	my %ad = %{$ad};

	open(OUT, ">$adfile");
	print OUT ">PrefixPE/1", "\n", "TACACTCTTTCCCTACACGACGCTCTTCCGATCT", "\n";
	print OUT ">PrefixPE/2", "\n", "GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT", "\n";
	foreach my $key (keys %ad) {
		print OUT ">" . $key, "\n", $ad->{$key}, "\n";
		}
	close(OUT);	

	my $call1 = system("java -jar $trimmomatic PE -threads 12 -phred33 $file1 $file2 $trimoutfiles{'p1'} $trimoutfiles{'u1'} $trimoutfiles{'p2'} $trimoutfiles{'u2'} ILLUMINACLIP:$adfile:2:40:15 SLIDINGWINDOW:4:20 MINLEN:$minLength LEADING:3 TRAILING:3");
# 	unlink($adfile);
	return(\%trimoutfiles)
	}
	
# input to cutadapt will be:
# 1 the hash %ad with the adapters for this library - excluding the PE adapters
# 2 each of the output files from trimmomatic, the names of which are stored in %trimoutfiles,
# 3 the the name of the reads file,
# 4 a suffix i.e. 'trim1' or 'cut1'

sub cutadapt {
	my ($ad,$in,$base,$suffix) = @_;
	my $out  = $base . '_' . $suffix;	
	my $curRead = $in;
	my $tracker = 1;
	my %ad = %{$ad};
		
	foreach my $key (keys %ad) {																		# for every adapter
		my $out = $curRead . $tracker; 																	#define the loop var $out to be the file name plus the tracker# 
		my $call1 = system("$cutadapt -b $ad{$key} -O 4 -n 5 -f fastq $curRead -o $out -m $minLength"); 	#do adaptor trim
		unlink($curRead) unless($curRead eq $in);														#if the file name is not the original (without a tracker appended), delete it - this is the output from cutadapt iterations not original reads or trimmomatic results
		$curRead = $out; 																				#redefine sub var $curRead as loop var $out (i.e. curRead + tracker) to enter the result of current loop in next iteration
		$tracker++;
		}
	
	my $call2 = system("mv $curRead $out"); # rename the file curRead with the sub level var $out
	return($out);
	}


sub mergeReads {
    my ($trimpairs,$base) = @_;
	my @trimpairs = @{$trimpairs};
	
	#create filenames for new read files
	my $newread1 = $file1 . '_' . $base . '_p1'; #unmerged from file1
    my $newread2 = $file2 . '_' .$base .'_p2'; #unmerged from file2
	my $newreadm = $outdir . $dir . '_' . $base . '_m'; #reads merged

	# MERGE READ PAIRS WITH FLASH
	my $call1 = system("$flash  $trimpairs[0] $trimpairs[1] -m 5 -x 0.01");
    
    # RENAME OUTPUT FILES AS NEWREAD FILES CREATED ABOVE
    my $move1 = $outdir . 'out.extendedFrags.fastq';
    my $move2 = $outdir . 'out.notCombined_1.fastq';
	my $move3 = $outdir . 'out.notCombined_2.fastq';
    my $call2 = system("mv $move1 $newreadm");
    my $call3 = system("mv $move2 $newread1");
    my $call4 = system("mv $move3 $newread2");
    
    #my $call2 = system("cat $reads{'u'} $outdir\.extendedFrags.fastq > $newreadu");
    #my $call5 = system("rm $lib\.extendedFrags.fastq $lib\.hist*");

	# COLLECT NEW FILE NAMES IN A HASH AND RETURN A REFERENCE TO IT
    my %newreads = ('1' => $newread1,'2' => $newread2, 'm' => $newreadm);
    return(\%newreads);

}
#############################################################################
# Add mate pair sub because cutadapt results in unequal files

sub fixMatePair {

	# INPUT IS REF TO FILENAMES FROM TRIMMING (@cleaned*) AND BASE FOR NAMING OUTFILES
	# OUTPUT IS A REF TO %newreads HAS WITH FILENAMES OF PAIRED AND UNPAIRED SEQS
	my ($readarray,$base) = @_;
	my @cleaned = @{$readarray};
	my %pair;	

	foreach my $reads (@cleaned) {
		open(IN, "<$reads");
		while(<IN>) {
			chomp(my $line = $_);
			# COUNT NUMBER OF SEQS WITH EACH ID
			if ($line =~ m/^(@\S+)/) {
			    $pair{$1}++;
			    chomp(my $seq = <IN>); chomp(my $qualid = <IN>); chomp(my $qual = <IN>);
				}
			}
		close(IN);	
		}

	# CREATE FILENAMES FOR OUTPUT FILES BASED ON ORIGINAL FILE NAMES WITH $BASE AND _p* OR _u APPENDED
	my $out1 = $file1 . '_' . $base . '_p1';
	my $out2 = $file2 . '_' . $base . '_p2';
	my $outu = $outdir . $dir . '_' . $base . '_u';
	open(OUT1, ">$out1");
	open(OUT2, ">$out2");
	open(OUTU, ">$outu");

	# HASH TO HOLD OUTPUT FILENAMES
	my %newpairs = ('1' => $out1, '2' => $out2, 'u' => $outu);

	# GO THROUGH EACH FILE IN @cleaned
	foreach my $reads (@cleaned) {
		open(IN, "<$reads");

		# NAME THE VAR $FILE AS 1 OR 2 BASED ON THE INPUT FILENAME FROM @cleaned ARRAY
		my $file = $1 if $reads =~ m/\w+paired(\d)_out_cleaned\d/;
		while(<IN>) {
			chomp(my $line = $_);
			
			if ($line =~ m/^(@\S+)/) {
				# CREATE ID FROM SEQ ID IN FILE FOR USE IN %pairs HASH
				my $id = $1;
				my $seq = <IN>;
				my $qualid = <IN>;
				my $qual = <IN>;
				
				# CHECK FOR DUPS IN %pairs HASH AND ASSIGN READS TO OUTFILES ACCORDINGLY
				if ($pair{$id} == 2) {
					if ($file == 1) {
						print OUT1 $id . "\n" . $seq . $qualid . $qual;
						}
					else {
						print OUT2 $id . "\n" . $seq . $qualid . $qual;
						}
					}
				else {
					print OUTU $id . "\n" . $seq . $qualid . $qual;
					}
				}	
			}
		close(IN);	
		}	
	close(OUT1); close(OUT2); close(OUTU);	
	return(\%newpairs);
	}


sub rc {
	my ($seq) = @_;
	my $rc = $seq;
	$rc = reverse($rc);
	$rc =~ tr/ATGCatgc/TACGtacg/;
	return($rc);
	}
	
}