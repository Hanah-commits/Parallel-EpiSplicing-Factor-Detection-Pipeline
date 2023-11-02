#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;

my ($scripts_dir, $results_dir);
my $timestamp = time;

########### MAIN PATHs - UPDATE IF NEEDED ######################################
BEGIN {
	$scripts_dir = "~/RBPmap_1.2/"; ## The main directory of RBPmap scripts and modules
	$results_dir = "~/EpiSplicing-Factor-Detection-Pipeline/" . $ARGV[0] . "_RBPmap/results" . basename($ARGV[2]); ## The main directory in which the results will be created. Chmod this directory 777 to enable writing permissions.
}
################################################################################
## MySQL connection definitions - UPDATE IF NEEDED: 
## (conf_file is optional and can contain the user and password instead of being hard-coded in the script)
my %mysql = (
	"host" => "localhost",
	"user" => "",
	"pass" => "",
	"port" => "3306",
	"conf_file" => $scripts_dir.".my.cnf", ## This file holds the user and password (instead of being hardcoded)
);
#################################################################################

use lib "$scripts_dir";
use RBPmap;
use SQL_RBPMAP;

##### Pathes and file names definitions
## Directories and files
my $UCSC_path = $scripts_dir . "UCSC/"; ## The main directory of the genome tables and UCSC utilities
my $motifs_dir = $scripts_dir . "RBP_PSSMs/";
my $working_dir = $results_dir . $timestamp . "/"; # The directory of the current run
my $logfile = $working_dir . "logfile.txt";
my $all_pred_file = $working_dir . "All_Predictions.txt";

##### Variables definitions
my %arg = (); # To hold the run parameters
## Default values
$arg{genome} = "human";
$arg{db} = "hg38";
$arg{input_type} = "coor";
$arg{window} = 50;
$arg{stringency} = "default";
$arg{pval_hit_consensus} = 0.01;
$arg{pval_sub_consensus} = 0.02;
$arg{pval_hit_PSSM} = 0.01;
$arg{pval_sub_PSSM} = 0.02;
$arg{COS} = "FALSE"; # In this version, the user cannot change this parameter
$arg{is_conservation} = "FALSE";
$arg{delete_fasta} = "TRUE";
my $MIN_SEQ_LENGTH = 21;
my $MAX_SEQ_LENGTH = 10000;
my $MAX_NUMBER_OF_SEQ = 5000;
my $MIN_MOTIF_LENGTH = 4;
my $MAX_MOTIF_LENGTH = 10;

## Genome <-> optional DB assemblies
my %genome_DB = (
	human => ["hg38", "hg19", "hg18"],
	mouse => ["mm10", "mm9"],
	drosophila => ["dm3"],
	other => [""],
);

my %origins = (
	human => "Hs/Mm",
	mouse => "Hs/Mm",
	drosophila => "Dm",
);

## Global variables
my $seq_method;
my @sequence = (); # An array of hashes (in the length of the processed sequence, after expansion and alignment) to hold all the sequence information. For each position: seq, position, coordinate, conservation score, state (for PSBM). (If no information: value is "NA")
my @Selected_Motifs = (); ## Should hold the proteins and motifs selected by the user with all the information saved in the %SF hash. For example: $Selected_Motifs[0]{protein} = NOVA1(human), $Selected_Motifs[0]{found_motif} = 0, $Selected_Motifs[0]{motifs}[0] = { motif => ycay, sf => NOVA1, origin => human, is_DB => 1, is_PSSM => 0}
my $db_motifs_str = "all_human"; # By default, selsct all human/mouse motifs
my $is_user_consensus = 0;
my $user_consensus_motifs = "";
my $is_user_pssm = 0;
my $user_pssm_file = "";
my $seq_working_dir;
my $pred_file;
my ($bedgraph_file, $bedgraph_for_user_file);
my $sequences;
my ($RBPmap, $SQL);
my ($ori, $chr, $start, $end) = ("", "", "", "");
my @t = (); # to hold the current time
my $seq_length;
my $err_string = "";
my $err_SQL = "";
my $found_motif_in_sequence = 0; # Should indicate if there are results for the sequence
my $total_seq_number = 0; # The total number of input sequences
my $valid_seq_number = 0; # The number of valid sequences
my $num_seq_with_results = 0; # Should indicate how many sequences have results in the current job (in order to know whether to display the All_Predictions file)
my ($selected_motifs, $new_motifs) = ("", ""); # To hold the selected motifs - in order to print in the results page
my $protein_counter = -1;

#########################
## Create an SQL object
$SQL = SQL_RBPMAP->new();

########################################################
### Create the working directory and set its permissions
if (system("mkdir $working_dir")) {
	print "\n". "*" x 50 . "\n";
	print "ERROR: mkdir $working_dir failed.\n$!\n";
	print "*" x 50 . "\n";
	exit;
}
if (system("chmod 777 $working_dir")) {
	print "\n". "*" x 50 . "\n";
	print "ERROR: chmod 777 $working_dir failed.\n$!\n";
	print "*" x 50 . "\n";
	exit;
}
## Open the logfile for writing
open (LOG, ">$logfile") or die "\nERROR: Cannot open the file $logfile for writing\n";

print LOG "Running Parameters:\n";
print LOG "-------------------\n\n";

#################################################
### Getting the parameters from the command line
my $proc = shift;
my $parameters = join (' ', @ARGV);
&parse_parameters($parameters);
	
################################################
#### Open and prepare the all-predictions file
unless (open (ALLPRED, ">$all_pred_file")){
	print LOG "\nError: Cannot open the main predictions output file $all_pred_file\n";
	print "\n". "*" x 50 . "\n";
	print "ERROR: Cannot open the main predictions output file $all_pred_file\n";
	print "*" x 50 . "\n";
	exit;
}
print ALLPRED "Predictions for job: $arg{job_name}\n";
print ALLPRED "Calculation parameters:\n";
if ($arg{genome} eq "other") {
	print ALLPRED "Genome: ".ucfirst($arg{genome})."\n";
}
else {
	print ALLPRED "Genome: ".ucfirst($arg{genome})." (".$arg{db}.")\n";
}

print ALLPRED "Selected motifs: $selected_motifs\n";
print ALLPRED "User-defined motifs: $new_motifs\n" if ($new_motifs ne "");
print ALLPRED "Stringency level: ".ucfirst($arg{stringency})."\n";
if ($arg{is_conservation} eq "TRUE") {
	print ALLPRED "Conservation filter: On\n";
}
else {
	print ALLPRED "Conservation filter: Off\n" ;
}
print ALLPRED "**********************************************\n";

############# Create a new RBPmap object ###########################################
$RBPmap = RBPmap->new();

##### Sending the general parameters to the RBPmap module
$RBPmap->set_mysql($mysql{host}, $mysql{user}, $mysql{pass}, $mysql{port}, $mysql{conf_file});
$RBPmap->set_general_parameters($arg{window}, $arg{pval_hit_consensus}, $arg{pval_sub_consensus}, $arg{pval_hit_PSSM}, $arg{pval_sub_PSSM}, $arg{genome}, $arg{db}, $arg{COS}, $arg{is_conservation}, \*LOG, $scripts_dir, $UCSC_path);

############# Calculation of background model and cutoff for each motif ################
print "\nRetrieving scores of background model...\nNote: the background model of user-defined motifs is calculated from scratch and may take a couple of minutes.\n";
foreach my $sf (@Selected_Motifs) {
	
	foreach my $query (@{$$sf{motifs}}){	
		
		## Set the match scoring function for the current motif
		my $match_function = &set_match_scoring_function($query);

		$err_string = $RBPmap->query($query, $match_function);
		if ($err_string ne "") {
			print "\n". "*" x 50 . "\n";
			print $err_string;
			print "*" x 50 . "\n";
			exit;
		}
	}
}
print "...Done!\n";
print "\n". "-" x 70 . "\n";

################################## Calculation per sequence ##################################################################
my $seq_count = 1;

### A loop over the sequences
foreach my $seq (@{$arg{seq}}) {
	
	### Create a directory and perform the calculation only for valid sequences
	if ($$seq{error} eq "") {
	
		$seq_length = length($$seq{sequence}); ## Calculate sequence length

		### Perform calculation
		print LOG "\n\nCalculating motifs for sequence: $$seq{title}\nSequence length: $seq_length\n";
		print "\nCalculating motifs for sequence: $$seq{title}\nSequence length: $seq_length\n";
		&sequence_calculation($seq, $seq_count);

		### Print the Predictions output file for the current sequence
		&print_prediction_file($seq);
				
		### Create the BedGraph file for the UCSC genome browser (only for human/mouse sequences with coordinates and when the total number of sequences <= $MAX_SEQ_FOR_BROWSER_DISPLAY)
		if ($ori =~ /[+-]/) {
			&print_bedgraph($seq);
		}
	}

	## The sequence has an error
	else {
		print LOG "\n\nSequence \'$$seq{title}\' cannot be analyzed\n$$seq{error}\n";

		print ALLPRED "\n$$seq{title}\n";
		print ALLPRED "==============================================================================\n";
		print ALLPRED "\n$$seq{error}\n";
		print ALLPRED "\n************************************************************************************************\n";
	}

	### Print the current sequence results
	&print_sequence_results($seq);

	print "-" x 70 . "\n";
	print LOG "-" x 70 . "\n";

	$seq_count++;
}
	
###########################################
### Delete the uploaded fasta files (if this option is on)
if ($arg{delete_fasta}) {

	my $fasta_files = $working_dir . "*.fasta";
	print LOG "\nRemoving fasta files...\n";
	print "\nRemoving fasta files...\n";
	print "-" x 70 . "\n";
	if (system ("rm -f $fasta_files")) {
		print LOG "\nsystem(\"rm -f $fasta_files\") has failed.\n$!\n";
		print "\nsystem(\"rm -f $fasta_files\") has failed.\n$!\n";
	}
}

close LOG;
close ALLPRED;

print "\nRBPmap has finished successfully.\n";
print "\nThe results can be found under: $working_dir\n";
print "\nThe main binding sites predictions file of all the sequences together: $all_pred_file\n\n";

exit(0);

#############################################################################################################
############################################### Functions ###################################################
#############################################################################################################
## Run the motifs calculation for the current sequence
sub sequence_calculation {
		
		my $seq = shift;
		my $seq_count = shift;

		my $is_conservation;
		my $match_function;
		my ($sf, $motif);
		my $found_motif;
		my $partial = 0;
		
		@sequence = ();
		$found_motif_in_sequence = 0;

		## Setting the output files names
		$seq_working_dir = $working_dir . "sequence" . $seq_count . "/";
		$pred_file = $seq_working_dir . "Predictions.txt";
		$bedgraph_file = $seq_working_dir . "bedgraph.txt";
		$bedgraph_for_user_file = $seq_working_dir . "RBPmap_custom_tracks.txt";

		## Create the working directory for the current sequence
		if (system("mkdir $seq_working_dir")) {
			print LOG "\nError: Cannot create directory $seq_working_dir\n";
			print "\n". "*" x 50 . "\n";
			print "ERROR: Cannot create directory $seq_working_dir\n";
			print "*" x 50 . "\n";
			exit;
		}	
		if (system("chmod 777 $seq_working_dir")) {
			print LOG "\nError: \'chmod 777 $seq_working_dir\' has failed\n";
			print "\n". "*" x 50 . "\n";
			print "ERROR: \'chmod 777 $seq_working_dir\' has failed\n";
			print "*" x 50 . "\n";
			exit;
		}

		### Analyze sequence / coordinates according to the input type
		($ori, $chr, $start, $end, $partial, $is_conservation, $$seq{warning}) = $RBPmap->sequence($arg{input_type}, $$seq{sequence}, $$seq{coor}, $$seq{title}, \@sequence, $seq_working_dir);

		if ($$seq{warning} ne "") {
			print "\nWarning: $$seq{warning}\n";
		}

		@t = localtime(time);
		print LOG "\n(sequence_calculation): $t[2]:$t[1]:$t[0] Starting calculation of score...\n";

		## A loop over the motifs
		foreach $sf (@Selected_Motifs){

			## Initialize the found_motif flag for the current protein in the current sequence 
			## (this flag indicates if a certain protein has at least one hit that passed the threshold in the current sequence - important for the results printing)
			$$sf{found_motif} = 0;

			foreach $motif (@{$$sf{motifs}}) {

				## Set the match scoring function for the current motif
				$match_function = &set_match_scoring_function($motif);
			
				$found_motif = $RBPmap->RUN_SFF($motif, \@sequence, length($$seq{sequence}), $is_conservation, $match_function);

				## If the current motif has hits in the current sequence, change the global flag found_motif_in_sequence (which indicates if there is any result for the current sequence) 
				if ($found_motif == 1 and $$sf{found_motif} == 0) {
					$$sf{found_motif} = 1;
					$found_motif_in_sequence = 1;
				}
			}
		}

		@t = localtime(time);
		print LOG "\n(sequence_calculation): $t[2]:$t[1]:$t[0] ...Done!\n";
}

####################################################################################################
###
sub ParseInput {

	my $sequences = shift;	
	
	my $err = "";
	my $i;
	my $array_ref;

	### Verify that the sequences string is not empty
	if ($sequences =~ /^\s*$/) {
		print LOG "\n(ParseInput): Error: no sequences - exiting...\n";
		print "\n". "*" x 50 . "\n";
		print "ERROR: no sequences - exiting...\n";
		print "*" x 50 . "\n";
		exit;
	}
	
	my @lines = split(/\n/, $sequences);
	
	## Determine whether the file is in FASTA format or coordinates and set the input_type parameter
	for ($i=0; $i<@lines ; $i++) {
		
		# Ignore empty lines
		if ($lines[$i] =~ /^\s*$/) {
			next;
		}

		# Format is FASTA -> call ParseSeq
		elsif ($lines[$i] =~ /^\s*>/ and $lines[$i+1] =~ /^[actugn-]+\s*$/i) {
			$arg{input_type} = "sequence";
			print LOG "Input type: sequence\n";
			$array_ref = ParSeq($sequences);
			last;
		}

		# Format is coordinates -> call ParseCoor
		elsif ($lines[$i] =~ /^(chr\w+)/) {

			if ($arg{genome} eq "other") {
				print "\n". "*" x 50 . "\n";
				print "ERROR: RBPmap cannot analyse genomic positions for genomes other than human, mouse or drosophila. For other genomes, the only acceptable input is sequences in FASTA format.\n";
				print "*" x 50 . "\n";
				print LOG "\nError: Genome is 'other' and input is in coordinates.\nExiting...\n";
				exit;
			}

			$arg{input_type} = "coor";
			print LOG "Input type: coordinates\n";
			$array_ref = ParseCoor($sequences);
			last;
		}

		# Unrecognized format -> print error
		else {
			print LOG "\n(ParseInput): Error: the input format is not recognized...\n";
			print "\n". "*" x 50 . "\n";
			print "ERROR: RBPmap cannot recognize the input format as either FASTA format or genomic coordinates.\n";
			print "*" x 50 . "\n";
			exit;
		}
	}

	return $array_ref;
}

#############################################################################################
### Separate sequences in fasta format
sub ParSeq{
	my $sequences = shift;
	my @seq = (); # array of hashes
	my $count = -1;
	my $error_seq_number = 0;
	
	my @lines = split(/\n/, $sequences);
	foreach my $line (@lines) {
		if ($line =~ /^\s*$/) {
			next;
		}
		if ($line =~ /^\s*>/) {
			###########################################
			### Validate the former sequence
			if ($count > -1) {
				&SeqValidation($seq[$count], "sequence");
			}
			###########################################
			$count++;
			$line =~ s/\s*$//; ## Remove only the whitespace characters in the end of the line (if any)
			$seq[$count]{title} = $line;
			$seq[$count]{title} =~ s/^\s*>//;
			$seq[$count]{sequence} = "";
			$seq[$count]{error} = "";
			$seq[$count]{coor} = "";
		}
		else {
			$line =~ s/\s//g;
			$line =~ s/u/t/gi;
			if ($count > -1) {
				$seq[$count]{sequence} .= $line;
			}
		}
	}

	if ($count > -1) {

		### Validate the last sequence
		&SeqValidation($seq[$count], "sequence");

		### Validate that the total number of sequences has not exceeded the maximum allowed
		$total_seq_number = @seq;
		if ($total_seq_number > $MAX_NUMBER_OF_SEQ) {
			print "\n". "*" x 50 . "\n";
			print "ERROR: Number of entries has exceeded the maximum allowed. Please divide your input into several RBPmap jobs, up to $MAX_NUMBER_OF_SEQ entries per job.\n";
			print "*" x 50 . "\n";
			print LOG "\nError: The number of entries ($total_seq_number) has exceeded the maximum allowed ($MAX_NUMBER_OF_SEQ).\n";
			exit;
		}
	}

	### No sequences at all -> print error page
	else {
		print LOG "\n(ParSeq): Error: the input format is not recognized... Exiting...\n";
		print "\n". "*" x 50 . "\n";
		print "ERROR: RBPmap could not recognize any valid sequence in FASTA format.\n";
		print "*" x 50 . "\n";
		exit;
	}

	### Assign the number of valid sequences (to be used later)
	foreach my $seq (@seq) {
		if ($$seq{error} eq "") {
			$valid_seq_number++;
		}
	}

=head
	##############################################################################
	## Print sequences to logfile
	foreach my $seq (@seq) {
		print LOG "Title: $$seq{title}\n";
		print LOG "Sequence: $$seq{sequence}\n";
		print LOG "Sequence length: ". length($$seq{sequence}) ."\n";
		if ($$seq{error} ne "") {
			print LOG "$$seq{error}\n\n";
		}
	}
=cut

	return(\@seq);
}

###########################################################################################################
### Parse the list of coordinates and get the sequences
sub ParseCoor{
	my $sequences = shift;
	my @seq = (); # array of hashes
	my $count = -1;
	my ($chr, $start, $end, $ori);	
	
	my @lines = split(/\n/, $sequences);
	foreach my $line (@lines) {

		# The line is empty
		if ($line =~ /^\s*$/) {
			next;
		}

		# Increment sequence counter and extract seq title
		else {
			$count++;
			$line =~ s/\s//g;
			$seq[$count]{title} = $line;
			$seq[$count]{title} =~ s/^\s*>//;
			$seq[$count]{sequence} = "";

			if ($line =~ /^(chr\w+)\s*:\s*(\S+)\s*-\s*(\S+)\s*:\s*([+-])\s*$/i) {

				$chr = $1;
				$start = $2;
				$end = $3;
				$ori = $4;

				$start = "1" if ($start eq "0"); # In case the start coor is '0', change it to '1' (this is also done in the genome browser - probably that's what the user meant)

				## Clean the start and end from commas (if any) and verify they are numbers
				$line =~ s/,//g;
				$start =~ s/,//g;
				$end =~ s/,//g;
				unless ($start =~ /^\d+$/ and $end =~ /^\d+$/) {
					$seq[$count]{error} = "ERROR: the coordinates format is invalid.\nRe-enter a genomic position in the following format: chromosome:start-end:strand\n";
					next;
				}

				## Get the sequences from the DB using the coordinates (send start-1 since GetSeq returns a sequence not including the start coor and including the end coor).
				$seq[$count]{coor} = $line;
				($seq[$count]{sequence}, $seq[$count]{error}) = &GetSeq($chr, $start-1, $end, $ori, $arg{db}, $working_dir);

				### Sequence validation
				if ($seq[$count]{error} eq "") {
					&SeqValidation($seq[$count], "coor");
				}
					
			}

			## The sequence is invalid - print an error message
			else {
				if ($seq[$count]{error} eq "") {
					$seq[$count]{error} = "ERROR: the coordinates format is invalid.\nRe-enter a genomic position in the following format: chromosome:start-end:strand\n";
				}
			}
		}
	}

	###########################################
	## Validate that the total number of sequences has not exceeded the maximum
	$total_seq_number = @seq;
	if ($total_seq_number > $MAX_NUMBER_OF_SEQ) {
		print "\n". "*" x 50 . "\n";
		print "ERROR: Number of entries has exceeded the maximum allowed. Please divide your input into several RBPmap jobs, up to $MAX_NUMBER_OF_SEQ entries per job.\n";
		print "*" x 50 . "\n";
		print LOG "\nError: The number of entries ($total_seq_number) has exceeded the maximum allowed ($MAX_NUMBER_OF_SEQ).\n";
		exit;
	}

	### Assign the number of valid sequences (to be used later)
	foreach my $seq (@seq) {
		if ($$seq{error} eq "") {
			$valid_seq_number++;
		}
	}

	return(\@seq);
}

###########################################################################################################
sub SeqValidation {

	my $ref = shift;
	my $format = shift;
	
	my $seq_length = length($$ref{sequence});

	## Verify legal characters (in 'sequence' format)
	if ($format eq "sequence") {
		unless ($$ref{sequence} =~ /^[actgn-]+$/i) {
			$$ref{error} = "ERROR: the sequence contains illegal characters. RBPmap supports valid DNA/RNA sequences in FASTA format.\n";
		}
	}
			
	if ($$ref{error} eq "") {
	
		## Verify that the sequence is not empty
		if ($seq_length == 0) {
			if ($format eq "sequence") {
				$$ref{error} = "ERROR: the sequence is empty. Please re-enter a valid sequence in FASTA format.";
			}
			else {
				$$ref{error} = "ERROR: the segment length is 0. Please re-enter a valid genomic position in coordinates format (chromosome:start-end:strand).";
			}
		}
			
		## Verify minimal length
		elsif ($seq_length < $MIN_SEQ_LENGTH) {
			$$ref{error} = "ERROR: the sequence is too short. The minimal sequence length is $MIN_SEQ_LENGTH bp.";
		}
			
		## Verify maximal length
		elsif ($seq_length > $MAX_SEQ_LENGTH) {
			$$ref{error} = "ERROR: the sequence exceeds $MAX_SEQ_LENGTH bp and cannot be calculated, please divide it into smaller segments.";
		}
	}
}

###########################################################################################################
# Get A sequence from coordinates. 
# Note: returns a sequence not including the start coordinate and including the end coor
# The sequence is returned in the right orientation (according to strand)
sub GetSeq{
	my ($chr, $str, $end, $ori, $DB, $working_dir)=@_;
	my $out = $working_dir . $chr."_".$str."_".$end.".fasta";
	my @seq = ();
	my $seq = "" ;
	my $err = "";
	my $real_start = $str+1; # To hold the correct start coordinate for an error massege

	my $nibfrag_path = $UCSC_path . "nibFrag";
	my $nib_path = $UCSC_path . $DB . "/" . $chr . ".nib";
	my $command = "$nibfrag_path -name=$chr:$str-$end:$ori $nib_path $str $end $ori $out";
	
	## The command has failed (probably something is wrong with coordinates)
	if (system($command)) {
		$err = "ERROR: Cannot extract the sequence from the following coordinates: $chr:$real_start-$end:$ori . Verify that the coordinates are correct and fit to the selected genome and database assembly.";
		return ($seq, $err);
	}

	## No output file (probably something is wrong with coordinates)
	if (-z $out or !(-e $out)) {
		$err = "ERROR: Cannot extract the sequence from the following coordinates: $chr:$real_start-$end:$ori . Verify that the coordinates are correct and fit to the selected genome and database assembly.";
		return ($seq, $err);
	}
	## OK
	else {
		open (IN,$out);
		while(<IN>){chomp;push (@seq,$_)}
		close(IN);
		foreach(@seq){if ($_ !~/>/){$seq=$seq.$_}}
		return ($seq, $err);
	}
}

#############################################################################################
#transpose matrix
sub transpose{
	my(@a)=@_;
	my  @b=();

	for my $i(0..$#a){
		for my $j(0..$#{$a[$i]}){
			$b[$j][$i]=$a[$i][$j];
		}
	}
	return(@b);
}

###########################################################################################
## Set an RGB color (without black, white or grey => R != G != B)
sub set_color {
	
	my ($R, $G, $B);
	my $color;

	$R = int(rand(200));
	$G = int(rand(200));
	$B = int(rand(200));

	$color = $R . "," . $G . "," . $B;
	return $color;
}

##########################################################################################################################################
sub print_prediction_file {
	
	my $seq = shift;

	my $chr_coor = "";

	unless (open (PRE, ">$pred_file")) {
		print LOG "\n(print_prediction_file): Error: Cannot open Predictions output file $pred_file for writing\n";
		print "\n". "*" x 50 . "\n";
		print "(print_prediction_file): ERROR: Cannot open Predictions output file $pred_file for writing\n";
		print "*" x 50 . "\n";
		exit;
	}
	unless (open (ALLPRED, ">>$all_pred_file")){
		print LOG "\n(print_prediction_file): Error: Cannot open All Predictions output file $all_pred_file for appending\n";
		print "\n". "*" x 50 . "\n";
		print "(print_prediction_file): ERROR: Cannot open All Predictions output file $all_pred_file for appending\n";
		print "*" x 50 . "\n";
		exit;
	}
	
	my $i = 0; # Counter of the proteins
	my $j = 0; # Counter of the positions
	my ($sf, $motif, $kmer, $WR_score, $Z_score, $P_value);

	print PRE "\nPredictions for sequence: $$seq{title}\n";
	print PRE "Calculation parameters:\n";
	if ($arg{genome} eq "other") {
		print PRE "Genome: ".ucfirst($arg{genome})."\n";
	}
	else {
		print PRE "Genome: ".ucfirst($arg{genome})." (".$arg{db}.")\n";
	}
	print PRE "Selected motifs: $selected_motifs\n";
	print PRE "User-defined motifs: $new_motifs\n" if ($new_motifs ne "");
	print PRE "Stringency level: ".ucfirst($arg{stringency})."\n";
	if ($arg{is_conservation} eq "TRUE") {
		print PRE "Conservation filter: On\n";
	}
	else {
		print PRE "Conservation filter: Off\n" ;
	}
	print PRE "======================================================================================================================\n";
	
	print ALLPRED "\n$$seq{title}\n";
	print ALLPRED "==============================================================================\n";


	#### There are results for the sequence => 
	#### 1. Increment $num_seq_with_results (important to know whether to display the All_Predictions file)
	#### 2. Print the results
	if ($found_motif_in_sequence) {
		$num_seq_with_results++;

		### A loop over the proteins
		for $i(0..$#Selected_Motifs) {

			$sf = $Selected_Motifs[$i]{protein};

			### There are results for this protein:
			### Print title
			if ($Selected_Motifs[$i]{found_motif}) {
				print PRE "\nProtein: $sf\n";
				print ALLPRED "\nProtein: $sf\n";

				### There can be 3 formats of printing:
				### 1. The genome is not defined - no coordinates and only WR_score
				if ($arg{genome} eq "other") {
					print PRE sprintf("%-17s\t%-10s\t%-10s\t%-5s\t%-6s\n", "Sequence Position", "Motif", "K-mer", "Score", "Cutoff");
					print ALLPRED sprintf("%-17s\t%-10s\t%-10s\t%-5s\t%-6s\n", "Sequence Position", "Motif", "K-mer", "Score", "Cutoff");
				}
				else {
					
					## 2. Have coordinates
					if ($ori =~ /[+-]/) {
						print PRE sprintf("%-17s\t%-18s\t%-10s\t%-10s\t%-7s\t%-8s\n", "Sequence Position", "Genomic Coordinate", "Motif", "K-mer", "Z-score", "P-value");
						print ALLPRED sprintf("%-17s\t%-18s\t%-10s\t%-10s\t%-7s\t%-8s\n", "Sequence Position", "Genomic Coordinate", "Motif", "K-mer", "Z-score", "P-value");
					}
					
					## 3. No coordinates
					else {
						print PRE sprintf("%-17s\t%-10s\t%-10s\t%-7s\t%-8s\n", "Sequence Position", "Motif", "K-mer", "Z-score", "P-value");
						print ALLPRED sprintf("%-17s\t%-10s\t%-10s\t%-7s\t%-8s\n", "Sequence Position", "Motif", "K-mer", "Z-score", "P-value");
					}
				}
				
				### A loop over the positions
				for $j(0..$#sequence) {
			
					## Position is in the original sequence
					if ($sequence[$j]{position} ne "NA") {

						## motif doesn't exceed the sequence length
						if ($sequence[$j]{position} <= $seq_length-length($sequence[$j]{$sf}{motif})+1) {

							$motif = $sequence[$j]{$sf}{motif};
							$kmer = $sequence[$j]{$motif}{kmer};
							$WR_score = $sequence[$j]{$sf}{WR_score};
					
							## Print hit if WR_score > 0 
							if($WR_score ne "NA" and $WR_score > 0){
						
								### There can be 3 formats of printing:
								### 1. The genome is not defined - no coordinates and only WR_score
								if ($arg{genome} eq "other") {
									print PRE sprintf("%-17s\t%-10s\t%-10s\t%-5s\t%-6s\n", $sequence[$j]{position}, $motif, $kmer, $WR_score, $sequence[$j]{$motif}{cutoff});
									print ALLPRED sprintf("%-17s\t%-10s\t%-10s\t%-5s\t%-6s\n", $sequence[$j]{position}, $motif, $kmer, $WR_score, $sequence[$j]{$motif}{cutoff});
								}
								else {

									$Z_score = $sequence[$j]{$sf}{Z_score};
									$P_value = $sequence[$j]{$sf}{P_value};

									## 2. Have coordinates
									if ($ori =~ /[+-]/) {
										$chr_coor = $chr . ":" . $sequence[$j]{coordinate};
										print PRE sprintf("%-17s\t%-18s\t%-10s\t%-10s\t%-7s\t%-.2e\n", $sequence[$j]{position}, $chr_coor, $motif, $kmer, $Z_score, $P_value);
										print ALLPRED sprintf("%-17s\t%-18s\t%-10s\t%-10s\t%-7s\t%-.2e\n", $sequence[$j]{position}, $chr_coor, $motif, $kmer, $Z_score, $P_value);
									}
									## 3. No coordinates
									else {
										print PRE sprintf("%-17s\t%-10s\t%-10s\t%-7s\t%-.2e\n", $sequence[$j]{position}, $motif, $kmer, $Z_score, $P_value);
										print ALLPRED sprintf("%-17s\t%-10s\t%-10s\t%-7s\t%-.2e\n", $sequence[$j]{position}, $motif, $kmer, $Z_score, $P_value);
									}
								}
							}
						}
						
						# Zero the positions in the downstream edge where the motif exceeds the sequence for the bedgraph file
						else {
							$sequence[$j]{$sf}{WR_score} = 0;
						}
					}
				}
			}
		}
	}
		
	### No motifs for the sequence
	else {
		print PRE "\nNo motifs found.\n";
		print ALLPRED "\nNo motifs found.\n";
	}

	print ALLPRED "\n************************************************************************************************\n";

}

##########################################################################################################################################
sub print_bedgraph {

	my $seq = shift;

	my ($i, $j) = (0,0);
	my $prev_coor = 0;
	my $motifs_of_SF = "";
	my $num_motifs_in_SF = 0;
	my $sf;
	my $Z_score;
	my $color;

	my %tracks = (
					"hg38" => {"conservation" => "cons20way", "genes" => "knownGene"},
					"hg19" => {"conservation" => "cons46way", "genes" => "knownGene"},
					"hg18" => {"conservation" => "cons44way", "genes" => "knownGene"},
					"mm10"  => {"conservation" => "cons60way", "genes" => "knownGene"},
					"mm9"  => {"conservation" => "cons30way", "genes" => "knownGene"},
					"dm3"  => {"conservation" => "multiz15way", "genes" => "flyBaseGene"},
	);

	### Open and print to the bedgraph file that will be displayed via the URL
	unless (open (BED, ">$bedgraph_file")) {
			print LOG "\n(print_bedgraph): Cannot open BedGraph file $bedgraph_file for writing\n";
			print "RBPmap has encountered a problem and cannot display the results in the UCSC genome browser.";
			return;
	}
	print BED "browser position ".$chr.":".$start."-".$end."\n";
	print BED "browser hide all\n";
	print BED "browser dense $tracks{$arg{db}}{genes}\n";
	print BED "browser full $tracks{$arg{db}}{conservation}\n" if ($arg{is_conservation} eq "TRUE");
	print BED "\n";

	### Open and print the bedgraph file that will be given to the user for download and contains RBPmap tracks only
	unless (open (USER_BED, ">$bedgraph_for_user_file")) {
			print LOG "\n(print_bedgraph): Cannot open BedGraph file $bedgraph_for_user_file for writing\n";
			print "RBPmap has encountered a problem and cannot provide the custom tracks file for download.";
			return;
	}
	print USER_BED "browser position ".$chr.":".$start."-".$end."\n\n";

	### Print the proteins tracks
	
	## A loop over the proteins
	for $i(0..$#Selected_Motifs){

		$sf = $Selected_Motifs[$i]{protein};
		
		## Write the protein in the BedGraph file only if it has at least one hit
		if ($Selected_Motifs[$i]{found_motif}) {

			$color = &set_color;
			$num_motifs_in_SF = @{$Selected_Motifs[$i]{motifs}};
			
			# Create a list of the motifs of the current protein 
			$motifs_of_SF = ""; 
			foreach  (@{$Selected_Motifs[$i]{motifs}}) {
				$motifs_of_SF .= "$$_{motif},";
			}
			$motifs_of_SF =~ s/,$//;

			if ($num_motifs_in_SF == 1) {
				print BED "track type=bedGraph name=\"$sf\" description=\"RBPmap results for $sf, motif: $motifs_of_SF\" visibility=full color=$color maxHeightPixels=30 autoScale=off viewLimits=0:5 windowingFunction=maximum yLineOnOff=on yLineMark=1.6\n";
				print USER_BED "track type=bedGraph name=\"$sf\" description=\"RBPmap results for $sf, motif: $motifs_of_SF\" visibility=full color=$color maxHeightPixels=30 autoScale=off viewLimits=0:5 windowingFunction=maximum yLineOnOff=on yLineMark=1.6\n";
			}
			else {
				print BED "track type=bedGraph name=\"$sf\" description=\"RBPmap results for $sf, motifs: $motifs_of_SF\" visibility=full color=$color maxHeightPixels=30 autoScale=off viewLimits=0:5 windowingFunction=maximum yLineOnOff=on yLineMark=1.6\n";
				print USER_BED "track type=bedGraph name=\"$sf\" description=\"RBPmap results for $sf, motifs: $motifs_of_SF\" visibility=full color=$color maxHeightPixels=30 autoScale=off viewLimits=0:5 windowingFunction=maximum yLineOnOff=on yLineMark=1.6\n";
			}	

			## '+' strand
			if ($ori eq "+") {
						
				## A loop over the positions
				for $j(0..$#sequence){

					$Z_score = $sequence[$j]{$sf}{Z_score};

					## position is part of the original sequence
					if ($sequence[$j]{position} ne "NA") {

						# Zero the uncalculated edges (NA)
						if($Z_score eq "NA"){ 
							$Z_score = 0;
						}
					
						$prev_coor = $sequence[$j]{coordinate} - 1;
						print BED "$chr $prev_coor $sequence[$j]{coordinate} $Z_score\n";
						print USER_BED "$chr $prev_coor $sequence[$j]{coordinate} $Z_score\n";
					}
				}
			}

			## '-' strand
			else {
						
				## A loop over the positions in a reverse order for the minus strand
				for ($j=$#sequence ; $j>=0 ; $j--) {
					
					$Z_score = $sequence[$j]{$sf}{Z_score};

					## position is part of the original sequence
					if ($sequence[$j]{position} ne "NA") {

						# Zero the uncalculated edges (NA)
						if($Z_score eq "NA"){ 
							$Z_score = 0;
						}
					
						$prev_coor = $sequence[$j]{coordinate} - 1;
						print BED "$chr $prev_coor $sequence[$j]{coordinate} $Z_score\n";
						print USER_BED "$chr $prev_coor $sequence[$j]{coordinate} $Z_score\n";
					}
				}
			}
			
			print BED "\n";
			print USER_BED "\n";
		}
	}

	close BED;
	close USER_BED;
}

#######################################################################################################
## Sets the match scoring function for the current motif
sub set_match_scoring_function {

	my $query_ref = shift;
	my $match_function;
	
	if ($arg{genome} ne "other" and $$query_ref{is_PSSM}) {
		$match_function = "pssm";
	}
	else {
		$match_function = "consensus";
	}

	return $match_function;
}

#####################################################################################################
sub parse_parameters {

	my $parameters = shift;
	my $upload_seq_file;
	my $found = 0;
	my $supported_db = "";

	#################################
	## Print help
	if ($parameters =~ /-help/ or $parameters eq "" or !defined($parameters)) {
		&print_usage;
		exit;
	}

	print "\nRunning Parameters:\n";
	print "-------------------\n";
	print LOG "\nRunning Parameters:\n";
	print LOG "-------------------\n";

	#############################################
	### Get the input file and save the sequences
	if ($parameters =~ /-input\s+(\S+)/i) {
		$upload_seq_file = $1;
		open (FILE, $upload_seq_file) or die "Cannot open the input file $upload_seq_file\n";
		while (<FILE>) { 
			$sequences .= $_;
		}
		close FILE;
		print "\nInput file: $upload_seq_file\n";
		print LOG "\nInput file: $upload_seq_file\n";
	}
	else {
		print "\n". "*" x 50 . "\n";
		print "ERROR: missing input file.\n";
		print "*" x 50 . "\n";
		&print_usage;
		exit;
	}

	##################
	## Parse the input and save the sequences is an array of hashes of the following form:
	## $arg{seq}[0]{title} = <seq 1 title>, $arg{seq}[0]{sequence} = <original seq 1>, $arg{seq}[0]{proc_sequence} = <processed seq 1 (expanded, after alignment)>, $arg{seq}[0]{error} = <error for the sequence>
	## Returns a reference to the array.
	$arg{seq} = &ParseInput($sequences);

	####################################
	### Get the genome and db parameters
	if ($parameters =~ /-genome\s+(\S+)/i) {
		$arg{genome} = $1;
		unless (exists($genome_DB{$arg{genome}})) {
			print LOG "\nError: the genome \'$arg{genome}\' is not supported. Available genomes: \'human\', \'mouse\', \'drosophila\' or \'other\'.\n";
			print "\n". "*" x 50 . "\n";
			print "ERROR: the genome \'$arg{genome}\' is not supported. Available genomes: \'human\', \'mouse\', \'drosophila\' or \'other\'.\n";
			print "*" x 50 . "\n";
			exit;
		}
	}

	if ($parameters =~ /-db\s+(\S+)/i) {
		$arg{db} = $1;

		## Verify that the genome and assembly fit together
		$found = 0;
		foreach  (@{ $genome_DB{$arg{genome}} }) {
			if ($arg{db} eq $_) {
				$found = 1;
			}
			$supported_db .= $_ . ", ";
		}
		$supported_db =~ s/,\s+$//;
		if ($found == 0) {
			print LOG "\nError: the genome and DB assembly don't fit.\nExiting...\n";
			print "\n". "*" x 50 . "\n";
			print "ERROR: the database assembly '$arg{db}' do not fit the genome '$arg{genome}'.\nSupported assemblies for $arg{genome}: $supported_db.\n";
			print "*" x 50 . "\n";
			exit;
		}
	}
	# Set it according to the defaults 
	else {
		$arg{db} = $genome_DB{$arg{genome}}[0];
	}
	print "\nGenome: $arg{genome}\n";
	print "\nDatabase assembly: $arg{db}\n" unless ($arg{genome} eq "other");
	print LOG "\nGenome: $arg{genome}\n";
	print LOG "\nDatabase assembly: $arg{db}\n" unless ($arg{genome} eq "other");

	####################################
	### Set the stringency
	if ($parameters =~ /-stringency\s+(\S+)/i) {
		$arg{stringency} = $1;
		if ($arg{stringency} eq "high"){
			$arg{pval_hit_consensus} = 0.005;
			$arg{pval_sub_consensus} = 0.01;
			$arg{pval_hit_PSSM} = 0.005;
			$arg{pval_sub_PSSM} = 0.01;
		}
		elsif ($arg{stringency} eq "default") {
			$arg{pval_hit_consensus} = 0.01;
			$arg{pval_sub_consensus} = 0.02;
			$arg{pval_hit_PSSM} = 0.01;
			$arg{pval_sub_PSSM} = 0.02;
		}
		else {
			print "\n". "*" x 50 . "\n";
			print "ERROR: the 'stringency' parameter is wrong. Available stringency levels: high/default\n";
			print "*" x 50 . "\n";
			exit;
		}
	}
	print LOG "\nConsensus Hit P-value: $arg{pval_hit_consensus}, Suboptimal P-value: $arg{pval_sub_consensus}\n";
	print LOG "\nPSSM Hit P-value: $arg{pval_hit_PSSM}, Suboptimal P-value: $arg{pval_sub_PSSM}\n";
	print "\nStringency level: $arg{stringency}\n";

	####################################
	### Set the conservation filter
	if ($parameters =~ /-conservation\s+(\S+)/i) {
		my $conservation = $1;
		if ($conservation eq "on") {
			if ($arg{genome} eq "other") {
				$arg{is_conservation} = "FALSE";
				print "\nConservation filter: off (cannot apply conservation for genome 'other')\n";
				print LOG "\nConservation filter: off\n";
			}
			else {
				$arg{is_conservation} = "TRUE";
				print "\nConservation filter: on\n";
				print LOG "\nConservation filter: on\n";
			}
		}
		elsif ($conservation eq "off") {
			$arg{is_conservation} = "FALSE";
			print "\nConservation filter: off\n";
		}
		else {
			print "\n". "*" x 50 . "\n";
			print "ERROR: the 'conservation' flag should be followed by on/off\n";
			print "*" x 50 . "\n";
			exit;
		}
	}

	########################
	### Set the job name
	if ($parameters =~ /-job_name\s+(\S+)/i) {
		$arg{job_name} = $1;
	}
	else {
		$arg{job_name} = $timestamp;
	}
	print "\nJob name: $arg{job_name}\n";
	print LOG "\nJob name: $arg{job_name}\n";

	#########################
	### Set the 'delete fasta' option
	if ($parameters =~ /-delete_fasta\s+(\S+)/i) {
		my $delete_fasta = $1;
		if ($delete_fasta eq "on") {
			$arg{delete_fasta} = "TRUE";
		}
		elsif ($delete_fasta eq "off") {
			$arg{delete_fasta} = "FALSE";
		}
		else {
			print "\n". "*" x 50 . "\n";
			print "ERROR: the 'delete_fasta' flag should be followed by on/off\n";
			print "*" x 50 . "\n";
			exit;
		}
	}

	##############################
	### Add the selected motifs to the main motifs array
	if ($parameters =~ /-db_motifs\s+(\S+)/i) {
		$db_motifs_str = $1;
		print LOG "\nDB motifs string: \'$db_motifs_str\'\n";
	}
	if ($parameters =~ /-consensus\s+(\S+)/i) {
		$is_user_consensus = 1;
		$user_consensus_motifs = $1;
		print LOG "\nUser-defined consensus motifs string: \'$user_consensus_motifs\'\n";
	}
	if ($parameters =~ /-pssm\s+(\S+)/i) {
		$is_user_pssm = 1;
		$user_pssm_file = $1;
		print LOG "\nUser-defined PSSM file: \'$user_pssm_file\'\n";
	}
	&AddMotifs;

	print "\n". "-" x 70 . "\n";
	print LOG "\n". "-" x 70 . "\n";
}

########################################################################################
### Adds new motifs added by the user to the selected motifs hash
sub AddMotifs {

	my $DB_motifs_ref;
	my @new_motifs = (); 
	my @SFs = ();
	my ($new_motif, $sf, $current_sf) = ("", "", "");
	my $motif_in_sf_counter = 0;
	my $i;
	my @list = ();

	print LOG "\nMotifs selection: ";
	print "\nMotifs selection: ";

	## Connect to RBPmap DB
	$err_SQL = $SQL->db('RBPmap', $mysql{host}, $mysql{user}, $mysql{pass}, $mysql{port}, $mysql{conf_file});
	if ($err_SQL ne "") {
		print "$err_SQL\n";
		exit;
	}

	### Select all the motifs stored in the DB
	if ($db_motifs_str eq "all") {
		($DB_motifs_ref, $err_SQL) = $SQL->FromSQL("SELECT * FROM Motifs_updated order by protein_name asc, motif_origin desc, motif_type desc");
		if ($err_SQL ne "") {
			print "$err_SQL\n";
			exit;
		}
		&add_DB_motifs($DB_motifs_ref);		
		$selected_motifs = "All RBPmap motifs";
		print LOG "All RBPmap motifs\n";
		print "All RBPmap motifs\n";
	}

	### Select all the human/mouse motifs
	elsif ($db_motifs_str eq "all_human") {
		($DB_motifs_ref, $err_SQL) = $SQL->FromSQL("SELECT * FROM Motifs_updated where motif_origin in ('human','mouse') order by protein_name asc, motif_type desc");
		if ($err_SQL ne "") {
			print "$err_SQL\n";
			exit;
		}
		&add_DB_motifs($DB_motifs_ref);		
		$selected_motifs = "All Human/Mouse motifs";
		print LOG "All Human/Mouse motifs\n";
		print "All Human/Mouse motifs\n";
	}

	### Select all the drosophila motifs
	elsif ($db_motifs_str eq "all_drosophila") {
		($DB_motifs_ref, $err_SQL) = $SQL->FromSQL("SELECT * FROM Motifs_updated where motif_origin in ('drosophila') order by protein_name asc, motif_type desc", \*LOG);
		if ($err_SQL ne "") {
			print "$err_SQL\n";
			exit;
		}
		&add_DB_motifs($DB_motifs_ref);		
		$selected_motifs = "All drosophila motifs";
		print LOG "All Drosophila motifs\n";
		print "All Drosophila motifs\n";
	}

	elsif ($db_motifs_str eq "none") {
		$selected_motifs = "None";
		print LOG "None\n";
		print "None\n";
	}

	### Custom selection - add motifs from the list
	else {
		
		## Verify that the proteins are given in the correct format
		unless ($db_motifs_str =~ /^\w+(,\w+)*/) {
			print "\n". "*" x 50 . "\n";
			print "Error: Invalid motifs selection.\n";
			print "*" x 50 . "\n";
			print "The available options for the selection of motifs from RBPmap database (-db_motifs <options>):\n";
			print "\t'all' - All RBPmap stored motifs.\n";
			print "\t'all_human' - All RBPmap human/mouse stored motifs.\n";
			print "\t'all_drosophila' - All RBPmap drosophila stored motifs.\n";
			print "\t<protein1,protein2,...> - Select all the motifs of the mentioned proteins names (primary names as appear in RBPmap database, without spaces).\n";
			print "\t'none' - None of RBPmap stored motifs (requires providing user-defined motifs).\n";
			exit;
		}

		## Quote each protein name in the string (for SQL)
		$db_motifs_str =~ s/(\w+)/\'$1\'/g; 

		## Select the motifs of the requested proteins from the database
		if ($arg{genome} eq "human" or $arg{genome} eq "mouse") {
			($DB_motifs_ref, $err_SQL) = $SQL->FromSQL("SELECT * FROM Motifs_updated where protein_name in ($db_motifs_str) and motif_origin in ('human','mouse') order by protein_name asc, motif_origin desc, motif_type desc");
			if ($err_SQL ne "") {
				print "$err_SQL\n";
				exit;
			}
		}
		elsif ($arg{genome} eq "drosophila") {
			($DB_motifs_ref, $err_SQL) = $SQL->FromSQL("SELECT * FROM Motifs_updated where protein_name in ($db_motifs_str) and motif_origin in ('drosophila') order by protein_name asc, motif_origin desc, motif_type desc");
			if ($err_SQL ne "") {
				print "$err_SQL\n";
				exit;
			}
		}
		else {
			($DB_motifs_ref, $err_SQL) = $SQL->FromSQL("SELECT * FROM Motifs_updated where protein_name in ($db_motifs_str) order by protein_name asc, motif_origin desc, motif_type desc", \*LOG);
			if ($err_SQL ne "") {
				print "$err_SQL\n";
				exit;
			}
		}
		&add_DB_motifs($DB_motifs_ref);
		print LOG "$selected_motifs\n";
		print "$selected_motifs\n";
	}

	##############################
	# New motifs
	
	### 1. Consensus motifs
	if ($is_user_consensus) {

		### Add motifs from the list
		@list = split(",", $user_consensus_motifs);
		
		## No motifs - print error
		if (@list == 0) {
			print "\n". "*" x 50 . "\n";
			print "ERROR: the user-defined consensus motifs list should be of the following format: 'protein1:motif1,protein2:motif2,...' (separated by commas with no spaces).\nFor example: '-consensus PTBP1:cucucu,MBNL1:ygcuky'.\n";
			print "*" x 50 . "\n";
			exit;
		}

		foreach (@list) {
			if ($_ =~ /(\S+):(\w+)/i) {
				push(@SFs, $1);
				push(@new_motifs, $2);
			}
			else {
				print "\n". "*" x 50 . "\n";
				print "ERROR: the user-defined consensus motifs list should be of the following format: 'protein1:motif1,protein2:motif2,...' (separated by commas with no spaces).\nFor example: '-consensus PTBP1:cucucu,MBNL1:ygcuky'.\n";
				print "*" x 50 . "\n";
				exit;
			}
		}
		
		$current_sf = "";
		## A loop over the new motifs
		for ($i=0 ; $i<=$#new_motifs ; $i++) {
			if ($new_motifs[$i] ne "") {
				$new_motif = $new_motifs[$i];
				$sf = $SFs[$i];

				unless ($sf =~ /^User/) {
					$sf  = "User_".$sf;
				}
				
				print LOG "\nNew consensus motif: $new_motif, Protein: $sf\n";
				print "\nNew consensus motif: $new_motif, Protein: $sf\n";
				
				# The motif contains illegal characters - print error and exit
				unless ($new_motif =~ /^[acgturymkswbdhvn]+$/i) {
					print "\n". "*" x 50 . "\n";
					print "ERROR: The user-defined motif \'$new_motif\' contains illegal characters. Consensus motifs may contain nucleotide symbols in IUPAC notation only.\n";
					print "*" x 50 . "\n";
					print LOG "Error: The user-defined motif \'$new_motif\' contains illegal characters.\nExiting...\n";
					exit;
				}
				# The motif exceeds the permitted length - print error and exit
				if (length($new_motif) > $MAX_MOTIF_LENGTH) {
					print "\n". "*" x 50 . "\n";
					print "ERROR: The user-defined motif \'$new_motif\' exceeds the permitted length (4-10 bp).\n";
					print "*" x 50 . "\n";
					print LOG "Error: The user-defined motif \'$new_motif\' exceeds the permitted length (4-10 bp).\nExiting...\n";
					exit;
				}
				elsif (length($new_motif) < $MIN_MOTIF_LENGTH) {
					print "\n". "*" x 50 . "\n";
					print "ERROR: The user-defined motif \'$new_motif\' is shorter than the permitted length (4-10 bp).\n";
					print "*" x 50 . "\n";
					print LOG "Error: The user-defined motif \'$new_motif\' is shorter than the permitted length (4-10 bp).\nExiting...\n";
					exit;
				}

				# New protein - add it to the Selected_Motifs array
				if ($sf ne $current_sf) {
					$current_sf = $sf;
					$protein_counter++;

					$Selected_Motifs[$protein_counter]{protein} = $sf;
					$Selected_Motifs[$protein_counter]{found_motif} = 0;
				}
					
				# Add the motif hash to the current protein in the Selected_Motifs array
				$new_motif = lc($new_motif);
				$new_motif =~ tr/t/u/;
				#push (@{$Selected_Motifs[$protein_counter]{motifs}}, {motif => $new_motif, sf => $sf, origin => "NA", is_PSSM => 0, is_DB => 0});
				push (@{$Selected_Motifs[$protein_counter]{motifs}}, {motif => $new_motif, sf => $sf, origin => "human", is_PSSM => 0, is_DB => 0});
				$new_motifs .= "$sf:$new_motif, ";
			}
		}
	}

	### 2. PSSM motif
	if ($is_user_pssm) {

		my @pssm = ();
		my ($A, $C, $G, $T);
		my @frequencies = ();
		my @consensus = ();
		my $consensus = "";
		my $pssm_counter = 0;

		open(PSSM, $user_pssm_file) or die "Cannot open the probability matrix file $user_pssm_file. Verify that the path is correct and the file is not empty.\n";

		$current_sf = "";
		while (<PSSM>) { 
			if ($_ =~ /^MOTIF\s+(\S+)?/) {
				
				## Save the previous pssm (if it's not the first one)
				if ($pssm_counter > 0) {

					# The PSSM exceeds the permitted length - print error and exit
					if (length($consensus) > $MAX_MOTIF_LENGTH) {
						print "\n". "*" x 50 . "\n";
						print "ERROR: The user-defined PSSM motif \'$sf\' exceeds the permitted length (4-10 bp).\n";
						print "*" x 50 . "\n";
						print LOG "Error: The user-defined PSSM motif \'$sf\' exceeds the permitted length (4-10 bp).\nExiting...\n";
						exit;
					}
					elsif (length($consensus) < $MIN_MOTIF_LENGTH) {
						print "\n". "*" x 50 . "\n";
						print "ERROR: The user-defined PSSM motif \'$sf\' is shorter than the permitted length (4-10 bp).\n";
						print "*" x 50 . "\n";
						print LOG "Error: The user-defined PSSM motif \'$sf\' is shorter than the permitted length (4-10 bp).\nExiting...\n";
						exit;
					}

					# New protein - add the protein to the Selected_Motifs array
					if ($sf ne $current_sf) {
						$current_sf = $sf;
						$protein_counter++;
						$motif_in_sf_counter = 0;

						$Selected_Motifs[$protein_counter]{protein} = $sf;
						$Selected_Motifs[$protein_counter]{found_motif} = 0;
					}

					# Add the motif hash to the current protein in the Selected_Motifs array
					#push (@{$Selected_Motifs[$protein_counter]{motifs}}, {motif => $consensus, sf => $sf, origin => "NA", is_PSSM => 1, is_DB => 0});
					push (@{$Selected_Motifs[$protein_counter]{motifs}}, {motif => $consensus, sf => $sf, origin => "human", is_PSSM => 1, is_DB => 0});
					push (@{$Selected_Motifs[$protein_counter]{motifs}[$motif_in_sf_counter]{PSSM}}, @pssm);
					print LOG "\nNew PSSM. Consensus: $consensus, Protein name: $sf\n";
					print "\nNew PSSM. Consensus: $consensus, Protein name: $sf\n";
					$new_motifs .= "$sf:$consensus, ";

					$motif_in_sf_counter++;
				}

				## Assign a motif name
				$pssm_counter++;
				@pssm = ();
				$consensus = "";

				# The protein name is defined in the file
				if (defined $1) {
					$sf = "User_".$1;
				}
				else {
					$sf = "User_PSSM_" . $pssm_counter;
				}

				print LOG "\nNew PSSM: $_";
			}
			elsif ($_ =~ /(\d\.\d+)\s+(\d\.\d+)\s+(\d\.\d+)\s+(\d\.\d+)/) {
				$A = $1;
				$C = $2;
				$G = $3;
				$T = $4;

				push(@pssm, {a => $A, c => $C, g => $G, t => $T, u => $T});

				# Create the consensus
				@frequencies = (
									{letter => "a", freq => $A},
									{letter => "c", freq => $C},
									{letter => "g", freq => $G},
									{letter => "u", freq => $T},
				);

				@consensus = sort { $$b{freq} <=> $$a{freq} } @frequencies;
				$consensus .= $consensus[0]{letter};

				print LOG "$A $C $G $T\n";
			}
		}
		close PSSM;
		
		## Save the last pssm
		if ($pssm_counter > 0) {

			# The PSSM exceeds the permitted length - print error and exit
			if (length($consensus) > $MAX_MOTIF_LENGTH) {
				print "\n". "*" x 50 . "\n";
				print "ERROR: The user-defined PSSM motif \'$sf\' exceeds the permitted length (4-10 bp).\n";
				print "*" x 50 . "\n";
				print LOG "Error: The user-defined PSSM motif \'$sf\' exceeds the permitted length (4-10 bp).\nExiting...\n";
				exit;
			}
			elsif (length($consensus) < $MIN_MOTIF_LENGTH) {
				print "\n". "*" x 50 . "\n";
				print "ERROR: The user-defined PSSM motif \'$sf\' is shorter than the permitted length (4-10 bp).\n";
				print "*" x 50 . "\n";
				print LOG "Error: The user-defined PSSM motif \'$sf\' is shorter than the permitted length (4-10 bp).\nExiting...\n";
				exit;
			}

			# New protein - add the protein to the Selected_Motifs array
			if ($sf ne $current_sf) {
				$current_sf = $sf;
				$protein_counter++;
				$motif_in_sf_counter = 0;

				$Selected_Motifs[$protein_counter]{protein} = $sf;
				$Selected_Motifs[$protein_counter]{found_motif} = 0;
			}

			# Add the motif hash to the current protein in the Selected_Motifs array
			#push (@{$Selected_Motifs[$protein_counter]{motifs}}, {motif => $consensus, sf => $sf, origin => "NA", is_PSSM => 1, is_DB => 0});
			push (@{$Selected_Motifs[$protein_counter]{motifs}}, {motif => $consensus, sf => $sf, origin => "human", is_PSSM => 1, is_DB => 0});
			push (@{$Selected_Motifs[$protein_counter]{motifs}[$motif_in_sf_counter]{PSSM}}, @pssm);
			print LOG "New PSSM. Consensus: $consensus, Protein: $sf\n";
			print "\nNew PSSM. Consensus: $consensus, Protein: $sf\n";
			$new_motifs .= "$sf:$consensus, ";
			$new_motifs .= "$sf:$consensus, ";
			
			$motif_in_sf_counter++;
		}

		## No motif was found (probably the file format is not correct) => print error page and exit
		else {
			print LOG "No PSSM motif was found in the file\n";
			print "\n". "*" x 50 . "\n";
			print "ERROR: no motif was found in the uploaded probability matrix file. Verify that the motif format is correct (MEME format).\n";
			print "*" x 50 . "\n";
			exit;
		}
	}

	$new_motifs =~ s/, $//;

	### Verify that there is at least one motif selected
	if (@Selected_Motifs == 0) {
		print LOG "Error: no valid motif was selected/provided.\n";
		print "\n". "*" x 50 . "\n";
		print "ERROR: no valid motif was selected/provided.\n";
		print "*" x 50 . "\n";
		print "\nMotif selection options:\n";
		print "\n-db_motifs <options>: Select one or more motifs from RBPmap database (default is 'all_human' - all RBPmap motifs).\n";
		print "\tAvailable options:\n";
		print "\t'all' - All RBPmap stored motifs.\n";
		print "\t'all_human' - All RBPmap human/mouse stored motifs.\n";
		print "\t'all_drosophila' - All RBPmap drosophila stored motifs.\n";
		print "\t<protein1,protein2,...> - Select all the motifs of the mentioned proteins names (primary names as appear in RBPmap database, without spaces).\n";
		print "\t'none' - None of RBPmap stored motifs (requires providing user-defined motifs).\n";
		print "\n-consensus <protein1:motif1,protein2:motif2,...>: An option to provide one or more user-defined consensus motifs. Valid motifs are 4-10 characters long and contain IUPAC symbols only, no spaces between the motifs.\n";
		print "\n-pssm <probability matrix file path>: An option to provide one or more user-defined PSSM motifs. The matrix file should be written in MEME format.\n";

		exit;
	}
}

#############################################################################################
sub add_DB_motifs {

	my $DB_motifs_ref = shift;
	my $i;
	my $motif_in_sf_counter = 0;
	my ($sf, $current_sf, $full_name, $origin, $motif, $motif_type) = ("", "", "", "", "", "");
	my @pssm = ();
	my $PSSM_file;
	my ($A, $C, $G, $T);

	for ($i=0 ; $i<=$#{$DB_motifs_ref} ; $i++) {
		$sf = $$DB_motifs_ref[$i][1];
		$origin = $$DB_motifs_ref[$i][2];
		$motif = $$DB_motifs_ref[$i][3];
		$motif_type = $$DB_motifs_ref[$i][4];

		$full_name = $sf."(".$origins{$origin}.")";
			
		## New protein - add the protein to the Selected_Motifs array
		if ($full_name ne $current_sf) {
			$motif_in_sf_counter = 0;
			$protein_counter++;
			$current_sf = $full_name;

			$Selected_Motifs[$protein_counter]{protein} = $full_name;
			$Selected_Motifs[$protein_counter]{found_motif} = 0;
		}

		## Add the motif

		## 1. PSSM - read from file
		if ($motif_type eq "pssm") {
			push (@{$Selected_Motifs[$protein_counter]{motifs}}, {motif => $motif, sf => $sf, origin => $origin, is_PSSM => 1, is_DB => 1});

			@pssm = ();
				
			$PSSM_file = $motifs_dir . $sf . "_" . $motif . "_" . $origin . "_PSSM.txt";
			unless (open (PSSM_FILE, $PSSM_file)) {
				print LOG "\nError: Cannot open $PSSM_file for reading\n";
				print "\n". "*" x 50 . "\n";
				print "ERROR: Cannot open $PSSM_file for reading\n";
				print "*" x 50 . "\n";
				exit;
			}
			while (<PSSM_FILE>) {
				if ($_ =~ /^\d+\s+(\d\.\d+)\s+(\d\.\d+)\s+(\d\.\d+)\s+(\d\.\d+)/) {
					$A = $1;
					$C = $2;
					$G = $3;
					$T = $4;

					push(@pssm, {a => $A, c => $C, g => $G, t => $T, u => $T});
				}
			}
			close PSSM_FILE;
				
			push (@{$Selected_Motifs[$protein_counter]{motifs}[$motif_in_sf_counter]{PSSM}}, @pssm);
		}

		## 2. Consensus
		else {
			push (@{$Selected_Motifs[$protein_counter]{motifs}}, {motif => $motif, sf => $sf, origin => $origin, is_PSSM => 0, is_DB => 1});
		}

		$motif_in_sf_counter++;
		$selected_motifs .= "$full_name:$motif, ";
	}
	$selected_motifs =~ s/, $//;
}

##########################################################################################################################
sub print_usage {
	print "\nUSAGE:\n\tRBPmap.pl -input <input file path> [options]\n";
	
	print "\n\tNote: by default, RBPmap searches for all the human/mouse motifs stored in RBPmap database.\n";

	print "\nINPUT:\n";
	print "\tRBPmap mandatory input is a DNA/RNA sequence or a list of sequences, to which the selected motifs are mapped. The query sequences can be given in two formats:\n";
	print "\t1. Sequences in FASTA format. For example:\n";
	print "\t\t>seq1\n\t\ttagataatttgacttgtttttactatta\n\t\t>seq2\n\t\tgtgtcctcccaaaaggcaacactggcagt\n";
	print "\t2. A list of genomic coordinates in the following format: chromosome:start-end:strand.\n\t\tFor example:\n";
	print "\t\tchr6:82518025-82518125:-\n\t\tchr15:97489441-97489541:+\n";

	print "\nOPTIONAL FLAGS (order is not important):\n----------------------------------------\n";
	
	print "\n\t-help: Print RBPmap manual\n";
	
	print "\nGENOME SELECTION:\n";
	print "\n\t-genome <'human'/'mouse'/'drosophila'/'other'>\n\t\tThe genome of the query sequences (default is 'human').\n";
	print "\n\t-db <'hg38'/'hg19'/'hg18'/'mm10'/'mm9'/'dm3'>\n\t\tThe database assembly of the query sequences (default is 'hg38').\n";
	
	print "\nMOTIF SELECTION:\n";
	print "\n\t-db_motifs <options>\n\t\tSelect one or more motifs from RBPmap database (default is 'all_human' - all RBPmap motifs). Available options:\n";
	print "\t\t'all' - All RBPmap stored motifs.\n";
	print "\t\t'all_human' - All RBPmap human/mouse stored motifs.\n";
	print "\t\t'all_drosophila' - All RBPmap drosophila stored motifs.\n";
	print "\t\t<protein1,protein2,...> - Select all the motifs of the mentioned proteins names (primary names as appear in RBPmap database, without spaces).\n";
	print "\t\t'none' - None of RBPmap stored motifs (requires providing user-defined motifs).\n";
	print "\n\t-consensus <protein1:motif1,protein2:motif2,...>\n\t\tAn option to provide one or more user-defined consensus motifs. Valid motifs are 4-10 characters long and contain IUPAC symbols only, no spaces between the motifs.\n\t\tFor example: -consensus PTBP1:cucucu,MBNL1:ygcuky\n";
	print "\n\t-pssm <probability matrix file path>\n\t\tAn option to provide one or more user-defined PSSM motifs. The matrix file should be written in MEME format.\n";
	
	print "\nADVANCED OPTIONS:\n";
	print "\n\t-stringency <'default'/'high'>\n\t\tControl the stringency level of the results.\n";
	print "\n\t-conservation <'on'/'off'>\n\t\tApply conservation filter on intergenic regions (default is: 'off').\n";

	print "\nGENERAL OPTIONS:\n";
	print "\n\t-job_name <name>\n\t\tA name for the output directory of the job, without spaces (default is timestamp).\n";
	print "\n\t-delete_fasta <'on'/'off'>\n\t\tDelete the intermediate fasta files which are created during RBPmap run (default is 'on').\n";
}

########################################################################################################################################
sub print_sequence_results() {

	my $seq = shift;

	my $seq_title = $$seq{title};

	print "\nResults for sequence $seq_title\n";

	## Sequence has an error and was not calculated
	if ($$seq{error} =~ /^Error/) {
		print "\n$$seq{error}\n";
	}

	## Sequence was calculated - print results
	else {

		### There are genomic coordinates
		if ($ori =~ /[+-]/) {
			print "Genomic position: $chr".":".$start."-".$end."   Strand: $ori\n"; 

			## No motifs found
			if ($found_motif_in_sequence == 0) {
				print "\nNo motifs found.\n";
			}

			## Print links
			else {
				print "\nOutput files for the sequence:\n";
				print "Binding sites predictions: $pred_file\n";
				print "RBPmap custom tracks in BedGraph format: $bedgraph_for_user_file\n";
			}
		}

		### No coordinats => no genome browser display
		else {

			## Genome is human/mouse => there was a problem with BLAT
			unless ($arg{genome} eq "other") {	
				print "Genomic position: N/A ($$seq{warning}).\n";
			}

			## No motifs found
			if ($found_motif_in_sequence == 0) {
				print "\nNo motifs found.\n";
			}

			## Print only the prediction summary link
			else {
				print "\nOutput files for the sequence:\n";
				print "Binding sites predictions: $pred_file\n";
			}
		}
	}
}


