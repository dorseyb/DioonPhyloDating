#############################################################################################################################
# A SCRIPT TO QUERY THE CDD FOR TRANSPOSONS IN PEP SEQUENCES
# SENDS A GI NUMBER FROM BLAST SEARCH TO CDD AND GREPS RESPONSE FOR 'TRANSPOS'
# PRINTS SEQ ID AND OUTPUT TO FILE UPON SUCCESSFUL GREP
# MODIFIED FROM EARLIER VERSION TO USE A SINGLE FILE OF BLAST RESULTS RATHER THAN INDIVIDUAL FILES FOR EACH ORTHOGROUP
#############################################################################################################################


use strict;
use warnings;
use IPC::System::Simple qw(run capture);
use LWP::UserAgent;

$|=1;
my $dir = '/Users/bdorsey/Documents/Dioon/RNAseq/trinity_files_processed/ortho_search/ortho/oneIndPerSpecies/I12Fasta/pep/';
# my $gilist = $dir.'transposons.txt';
my $trans = $dir.'transposon_hits.txt';
my $url = "http://www.ncbi.nlm.nih.gov/Structure/bwrpsb/bwrpsb.cgi";
my $blastfile = $dir.'conseqs.all.blastp';
my $error = $dir.'trans_check_error.txt';
my %qhash = ();

# GET GI FOR EACH BLAST HIT IN BLAST RESULTS FILE AND MAKE A HASH OF OG->GI
open my $bfh, "<$blastfile" or die $!;

foreach my $line (<$bfh>) {
# my $line = 'oneInd.og12_1000pep.afa	gi|116780190|gb|ABK21582.1|	unknown [Picea sitchensis]	329	394	350	73.25	3e-175';
	chomp $line;
	my @a = split /\s/, $line;
	my $q = $a[0];
	my $gi = ( split(/\|/, $a[1]) )[1];
	$qhash{$q} = $gi;
}
# print "$q\n$gi\n";

# FILE TO PRINT RESULTS TO
open my $transfh, ">$trans" or die $!;
print $transfh "Orthogroup:Query	Hit type	PSSM-ID	From	To	E-Value	Bitscore	Accession	Short name	Incomplete	Superfamily	Definition\n";

# FILE TO PRINT ERROR MESSAGES
open my $errorfh, ">>$error" or die $!;

# SET COUNTS OF QUERIES AND HITS
my $filecount = 0;
my $transcount = 0;

# SUBMIT QUERY TO CDD AND CHECK FOR SUCCESSFUL SUBMISSION
QUERY: foreach my $f (keys %qhash) {
	$filecount++;
	
	my $query = $qhash{$f};
	my $browser = LWP::UserAgent->new;
	my $response = $browser->post(
		$url,
		[
		'smode'	=> 'prec',
		'useid1' => '0',
		'cddefl' => 'true',
		'qdefl'  => 'false',
		'dmode'  => 'rep',
		'tdata'  => "hits",
# 			( map {; queries => $_ } @queries )
		'queries' => $query
		],
	  );
	
	print $errorfh "Error: ", $response->status_line
	unless $response->is_success;
	my $cdsid;
	if($response->content =~ /^#cdsid\s+([a-zA-Z0-9-]+)/m) {
		$cdsid =$1;
	} else {
		print $errorfh "$f: Submitting the search failed, can't make sense of response:\t$response->content\n";
		next QUERY;
	}	

# CHECK FOR RESPONSE AND GATHER OUTPUT
	my $done = 0;
	my $status = -1;
	while ($done == 0) {
		my $time = 0;
		sleep(2);
				
		my $browser = LWP::UserAgent->new;
		my $response = $browser->post(
			$url,
			[
				'tdata' => "hits",
				'cdsid' => $cdsid
			],
		);
		print $errorfh "$f\tError: ".$response->status_line."\n"
			unless $response->is_success;

		if ($response->content =~ /^#status\s+([\d])/m) {
			$status = $1;
			if ($status == 0) {
				$done = 1;
	# 				print "\nSearch has been completed, retrieving results ..\n";
			} elsif ($status == 3) {
				next;				
			} elsif ($status == 1) {
				print $errorfh $f."\tInvalid request ID\n";
				next FILE;
			} elsif ($status == 2) {
				print $errorfh "$f\tInvalid input - missing query information or search ID\n";
	# 				print "\n$f: Invalid input - missing query information or search ID\n";
				next FILE;
			} elsif ($status == 4) {
				print $errorfh "$f\tQueue Manager Service error\n";
	# 				print "\n$f: Queue Manager Service error\n";
				next FILE;
			} elsif ($status == 5) {
				print $errorfh "$f\tData corrupted or no longer available\n";
	# 				print "\n$f: Data corrupted or no longer available\n";
				next FILE;
			} else {
	# 				print "\n$f: Checking search status failed, can't make sense of response";
				print $errorfh "$f\nChecking search status failed, can't make sense of response:\n";
				next FILE;
			}

			if ($response->content =~ /(^Q#\d.*transpos.*$)/m){
				print $transfh "$f\t$1\n";
				$transcount++;
			}
		
		}
	}		
	print "\r$transcount/$filecount";
}

close $errorfh or die $!;
close $transfh or die $!;
print "\n";
