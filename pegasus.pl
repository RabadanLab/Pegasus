#!/usr/bin/env perl

use strict;
use warnings FATAL => 'all';
use Getopt::Long;
use Scalar::Util qw(looks_like_number);
use File::Basename;

my $data_input_file;
my $config_file;
my $log_folder;
my $out_folder;
my $sge = 0;
my $cmd = "";
my $res = 0;
my $keepDB = 0;
my $clf_model = "gbc";
my %config = ();

GetOptions (
	'c=s' => \$config_file, 
	'd=s' => \$data_input_file,
	'l=s' => \$log_folder,
	'o=s' => \$out_folder, 
	'p=i' => \$sge, 
	'k=i' => \$keepDB, 
        'm=s' => \$clf_model
);  

my $p_fail = 0;

if(!defined($config_file))
{
        print STDERR "Error - $0: config_file not specified\n";
	$p_fail=1;
}

if(!defined($data_input_file))
{
        print STDERR "Error - $0: data_input_file not specified\n";
	$p_fail=1;
}

if(!defined($log_folder))
{
        print STDERR "Error - $0: log_folder not specified\n";
	$p_fail=1;
}

if(!defined($out_folder))
{
        print STDERR "Error - $0: output_folder option not specified\n";
	$p_fail=1;
}

if($sge<0 || $sge>1)
{
        print STDERR "Error - $0: parallelization_sge parameter allows only 0 and 1 values\n";
	$p_fail=1;
}

if($keepDB<0 || $keepDB>1)
{
        print STDERR "Error - $0: keepDB parameter allows only 0 and 1 values\n";
	$p_fail=1;
}

if(!($clf_model eq "gbc" || $clf_model eq "rfc"))
{
        print STDERR "Error - $0: model_classification parameter must be in {gbc, rfc}\n";
        $p_fail=1;        
}

if ($p_fail)
{
        printHelp();
        exit(255);
}

open LOG, ">>".$log_folder."/Pegasus.log" or die "INFO - $0: Error opening Pegasus log file. \n";

open CONF, "<".$config_file or die "INFO - $0: Error opening configuration file. \n";

while(my $line = <CONF>)
{
	chomp($line);
	if($line =~ m/^$/ || $line =~ m/^#/)
	{
		next;
	}

	my @fields = split(/\t/, $line);
	my $key = $fields[0];
	my $value = $fields[1];

	$config{$key} = $value;

}

close(CONF);


die "pegasus_folder parameter not specified in configuration file" if !defined($config{'pegasus_folder'});
die "hg parameter not specified in configuration file" if !defined($config{'hg'});
die "hg_fai parameter not specified in configuration file" if !defined($config{'hg_fai'});
die "gtf_file parameter not specified in configuration file" if !defined($config{'gtf_file'});
#die "script parameter not specified in configuration file" if !defined($config{'script'});
#die "jar parameter not specified in configuration file" if !defined($config{'jar'});
die "sample_type parameter not specified in configuration file" if !defined($config{'sample_type'});
die "sge_param parameter not specified in configuration file" if !defined($config{'sge_param'});

die "Human genome reference file not found at $config{'hg'}" if ! -e $config{'hg'};
die "Human genome reference file is empty" if -z $config{'hg'};

die "Human genome reference index file not found at $config{'hg_fai'}" if ! -e $config{'hg_fai'};
die "Human genome reference index file is empty" if -z $config{'hg_fai'};

die "Human genome reference GTF file not found at $config{'gtf_file'}" if ! -e $config{'gtf_file'};
die "Human genome reference GTF file is empty" if -z $config{'gtf_file'};

die "Log folder $log_folder does not exist" if ! -e $log_folder;
die "Output folder $out_folder does not exist" if ! -e $out_folder;

$config{'script'}=$config{'pegasus_folder'}."/source/scripts";
$config{'resources'}=$config{'pegasus_folder'}."/resources";
$config{'jar'}=$config{'pegasus_folder'}."/jars";

# Checking FASTA and GTF reference file
if(! -e $log_folder."/InputChecking.job")
{
	print STDERR "[".`date | tr '\n' ' '`."] Checking input reference files:\n";
	print STDERR "[".`date | tr '\n' ' '`."] $config{'gtf_file'} file:\n";

	open( GTF, "<".$config{'gtf_file'} ) or die "$!";

	while(my $line = <GTF>)
	{
		chomp($line);
		if( $line !~ m/^chr/ )
		{
			die "Error: the chromosome format is not in the form of chr[chromosome number] in the specified $config{'gtf_file'} file.";
		}	
	}

	close(GTF);

	print STDERR "[".`date | tr '\n' ' '`."] $config{'hg'} file:\n";
	open( HG, "<".$config{'hg'} ) or die "$!";

	while(my $line = <HG>)
	{
		chomp($line);	
		if( $line =~ m/^>/ && $line !~ m/chr/ )
		{
			die "Error: the chromosome format is not in the form of chr[chromosome number] in the specified $config{'hg'} file.";
		}	
	}

	close(HG);

	checkpoint($log_folder."/InputChecking.job");
}
else
{
	print STDERR "[".`date | tr '\n' ' '`."] Checking input reference files [SKIPPED]\n";
}








#if(! -e $log_folder."/QueryFusionDatabase_Delete.job" && $keepDB==0)
if($keepDB==0)
{
	$cmd = "java -Xmx1024m -jar "." $config{'jar'}"."/QueryFusionDatabase.jar -t FUSIONS_COMPLETE_ID -c deleteAll -d "." $config{'pegasus_folder'}"."/resources/hsqldb-2.2.7/hsqldb/mydb";
	$cmd .= " >> ".$log_folder."/QueryFusionDatabase_Delete.log 2>> ".$log_folder."/QueryFusionDatabase_Delete.log";
	print STDERR "[".`date | tr '\n' ' '`."] Cleaning fusion Database\n";
	print LOG "[".`date | tr '\n' ' '`."] $cmd\n\n";
	$res = system($cmd);
	die "Error running QueryFusionDatabase.jar. Exit code $res." if $res!=0;
	checkpoint($log_folder."/QueryFusionDatabase_Delete.job");
}
else
{
	print STDERR "[".`date | tr '\n' ' '`."] Cleaning fusion Database [SKIPPED]\n";
}








if(! -e $log_folder."/LoadFusionReport.job")
{
	open DATAIN, "<".$data_input_file or die "INFO - $0: Error opening data input file. \n";

	while(my $line = <DATAIN>)
	{
		chomp($line);
		if($line =~ m/^$/ || $line =~ m/^#/)
		{
			next;
		}

		my @fields = split(/\t/, $line);
		my $sample_path = $fields[0];
		my $sample_name = $fields[1];
		my $sample_type = $fields[2];
		my $fusion_program = $fields[3];

		if ($fusion_program eq "general" || $fusion_program eq "defuse" || $fusion_program eq "chimerascan" || $fusion_program eq "bellerophontes")
		{
			printf LOG "[".`date | tr '\n' ' '`."] $line - elaborated\n\n";
		}
		else
		{
			die "Error. The specified fusion program has not been supported yet.\n ===> $line\n" ;
		}

		$cmd  = "java -Xmx1024m -jar "."$config{'jar'}"."/LoadFusionReport.jar -f "."$config{'pegasus_folder'}"."/resources/hsqldb-2.2.7/hsqldb/mydb -i ";
		$cmd .= "$sample_path"." -t FUSIONS_COMPLETE_ID -d xdb -s \""."$sample_name"."|"."$sample_type"."\" -p ".$fusion_program;
		$cmd .= " >> ".$log_folder."/LoadFusionReport.log 2>> ".$log_folder."/LoadFusionReport.log";
		printf STDERR "[".`date | tr '\n' ' '`."] Loading fusion report %s\n", $sample_name;
		printf LOG "[".`date | tr '\n' ' '`."] $cmd\n\n";
		$res = system($cmd);
		die "Error running LoadFusionReport.jar. Exit code $res." if $res!=0;


	}

	close(DATAIN);
	checkpoint($log_folder."/LoadFusionReport.job");
}
else
{
	print STDERR "[".`date | tr '\n' ' '`."] Loading fusion report [SKIPPED]\n";
}











if(! -e $log_folder."/QueryFusionDatabase.job")
{
	$cmd  = "java -Xmx1024m -jar ".$config{'jar'}."/QueryFusionDatabase.jar -t FUSIONS_COMPLETE_ID -c all_with_kinases ";
	$cmd .= "-s ".$config{'sample_type'}." -d ".$config{'pegasus_folder'}."/resources/hsqldb-2.2.7/hsqldb/mydb";
	$cmd .= " > ".$out_folder."/samples_info.txt 2>> ".$log_folder."/QueryFusionDatabase.log";
	printf STDERR "[".`date | tr '\n' ' '`."] Unifying reports for sample type %s\n", $config{'sample_type'};
	printf LOG "[".`date | tr '\n' ' '`."] $cmd\n\n";
	$res = system($cmd);
	die "Error running QueryFusionDatabase.jar. Exit code $res." if $res!=0;
	die "Error Running QueryFusionDatabase.jar . samples_info.txt file not generated." if ! -e $out_folder."/samples_info.txt";
	die "Error Running QueryFusionDatabase.jar . samples_info.txt is empty." if -z $out_folder."/samples_info.txt";
	checkpoint($log_folder."/QueryFusionDatabase.job");
}
else
{
	printf STDERR "[".`date | tr '\n' ' '`."] Unifying reports for sample type %s [SKIPPED]\n", $config{'sample_type'};
}
















if(! -e $log_folder."/format_forInframeAnalysis.job")
{
	$cmd  = $config{'script'}."/format_forInframeAnalysis.sh ".$out_folder."/samples_info.txt ";
	$cmd .= "sort -u > ".$out_folder."/samples_info.formatted.sorted.txt";
	printf STDERR "[".`date | tr '\n' ' '`."] Preparing the environment for annotation\n";
	printf LOG "[".`date | tr '\n' ' '`."] $cmd\n\n";
	$res = system($cmd);
	die "Error running format_forInframeAnalysis.sh . Exit code $res." if $res!=0;
	die "Error Running format_forInframeAnalysis.sh . samples_info.formatted.sorted.txt file not generated." if ! -e $out_folder."/samples_info.formatted.sorted.txt";
	die "Error Running format_forInframeAnalysis.sh . samples_info.formatted.sorted.txt is empty." if -z $out_folder."/samples_info.formatted.sorted.txt";
	checkpoint($log_folder."/format_forInframeAnalysis.job");
}
else
{
	print STDERR "[".`date | tr '\n' ' '`."] Preparing the environment for annotation [SKIPPED]\n";
}












if(! -e $log_folder."/RunAnnotation.job")
{
	if($sge==1)
	{
		printf STDERR "[".`date | tr '\n' ' '`."] Preparing the environment for SGE running\n";
		
		$cmd  = "mkdir -p ".$out_folder."/split.folder.input ";
		printf LOG "[".`date | tr '\n' ' '`."] $cmd\n\n";
		system($cmd);

		$cmd  = "mkdir -p ".$out_folder."/split.folder.output ";
		printf LOG "[".`date | tr '\n' ' '`."] $cmd\n\n";
		system($cmd);

		$cmd   = $config{'script'}."/split.sh ";
		$cmd  .= $out_folder."/samples_info.formatted.sorted.txt ".$out_folder."/split.folder.input 10";
		printf LOG "[".`date | tr '\n' ' '`."] $cmd\n\n";
		system($cmd);

		printf STDERR "[".`date | tr '\n' ' '`."] Annotating fusions on SGE system\n";
		
		$cmd   = $config{'script'}."/run_annotate.sh ";
		$cmd  .= $out_folder."/split.folder.input ";
		$cmd  .= $out_folder."/split.folder.output ";
		$cmd  .= $config{'pegasus_folder'}." ";
		$cmd  .= $log_folder." ";
		$cmd  .= $config{'hg'}." ";
		$cmd  .= $config{'gtf_file'}." ";
		$cmd  .= " ".$config{'sge_param'}." " ;
		$cmd .= " > ".$log_folder."/RunAnnotation.log 2>> ".$log_folder."/RunAnnotation.log";
		printf LOG "[".`date | tr '\n' ' '`."] $cmd\n\n";
		system($cmd);

		$cmd   = "cat ".$out_folder."/split.folder.output/* >  ";
		$cmd  .= $out_folder."/samples_info.annotated.txt ";
		printf LOG "[".`date | tr '\n' ' '`."] $cmd\n\n";
		$res = system($cmd);
		die "Error Running RunAnnotation procedure . samples_info.annotated.txt file not generated." if ! -e $out_folder."/samples_info.annotated.txt";
		die "Error Running RunAnnotation procedure . samples_info.annotated.txt is empty." if -z $out_folder."/samples_info.annotated.txt";
		checkpoint($log_folder."/RunAnnotation.job");

	}
	else
	{
		printf STDERR "[".`date | tr '\n' ' '`."] Annotating fusions\n";
		$cmd   = $config{'script'}."/MultiFusionSequenceFromGTF_wrap.sh -r ";
		$cmd  .= $config{'hg'}." -g ";
		$cmd  .= $config{'gtf_file'}." -s ";
		$cmd  .= $config{'script'}." -j ";
		$cmd  .= $config{'jar'}." -o ";
		$cmd  .= $out_folder."/samples_info.annotated.txt -t ";
		$cmd  .= $out_folder."/tmp -i ";
		$cmd  .= $out_folder."/samples_info.formatted.sorted.txt ";
		$cmd  .= " > ".$log_folder."/RunAnnotation.log 2>> ".$log_folder."/RunAnnotation.log";
		printf LOG "[".`date | tr '\n' ' '`."] $cmd\n\n";
		$res = system($cmd);
		die "Error Running RunAnnotation procedure . samples_info.annotated.txt file not generated." if ! -e $out_folder."/samples_info.annotated.txt";
		die "Error Running RunAnnotation procedure . samples_info.annotated.txt is empty." if -z $out_folder."/samples_info.annotated.txt";
		checkpoint($log_folder."/RunAnnotation.job");		
	}
}
else
{
	print STDERR "[".`date | tr '\n' ' '`."] Annotating fusions on SGE system [SKIPPED]\n";
	print STDERR "[".`date | tr '\n' ' '`."] Preparing the environment for SGE running [SKIPPED]\n";
}









if(! -e $log_folder."/Wrapper_final_report.job")
{
	printf STDERR "[".`date | tr '\n' ' '`."] Preparing files for ML module\n";
		
	$cmd   = "perl ".$config{'script'}."/merge_res_commonFormat.pl -i ";
	$cmd  .= $out_folder."/samples_info.txt -r ";
	$cmd  .= $out_folder."/samples_info.annotated.txt |  cut -f1-23,28- | cut -f1-26,28- > ";
	$cmd  .= $out_folder."/final_results_forXLS.txt";
	$cmd .= " 2>> ".$log_folder."/Merge_format.log";
	print LOG "[".`date | tr '\n' ' '`."] $cmd\n\n";
	$res = system($cmd);
	die "Error running merge_res_commonFormat.pl. Exit code $res." if $res!=0;
	die "Error Running merge_res_commonFormat.pl . final_results_forXLS.txt file not generated." if ! -e $out_folder."/final_results_forXLS.txt";
	die "Error Running merge_res_commonFormat.pl . final_results_forXLS.txt is empty." if -z $out_folder."/final_results_forXLS.txt";
	checkpoint($log_folder."/Merge_format.job");


	$cmd   = "perl ".$config{'script'}."/wrapper_finalreport.pl -s ";
	$cmd  .= $config{'script'}." -d ";
	$cmd  .= $out_folder."/ -i ";
	$cmd  .= $out_folder."/final_results_forXLS.txt -r ";
	$cmd  .= $config{'resources'}." -l ";
	$cmd  .= $log_folder."/Wrapper_final_report.log -o ";
	$cmd  .= $out_folder."/final_results_forXLS.ML.input.txt ";
	$cmd  .= " >> ".$log_folder."/Wrapper_final_report.log 2>> ".$log_folder."/Wrapper_final_report.log";
	print LOG "[".`date | tr '\n' ' '`."] $cmd\n\n";
	$res = system($cmd);
	die "Error running wrapper_finalreport.pl . Exit code $res." if $res!=0;
	die "Error Running wrapper_finalreport.pl . final_results_forXLS.ML.input.txt file not generated." if ! -e $out_folder."/final_results_forXLS.ML.input.txt";
	die "Error Running wrapper_finalreport.pl . final_results_forXLS.ML.input.txt is empty." if -z $out_folder."/final_results_forXLS.ML.input.txt";
	checkpoint($log_folder."/Wrapper_final_report.job");
}
else
{
	print STDERR "[".`date | tr '\n' ' '`."] Preparing files for ML module [SKIPPED]\n";
}












if(! -e $log_folder."/ML_module.job")
{
	printf STDERR "[".`date | tr '\n' ' '`."] Running ML module\n";

        $cmd   = "python ".$config{'script'}."/classify.py -i ";
	$cmd  .= $out_folder."/final_results_forXLS.ML.input.txt -m ";
	$cmd  .= $config{'pegasus_folder'}."/learn/models/trained_model_".$clf_model.".pk -o ";
	$cmd  .= $out_folder."/pegasus.output.txt -l ";
	$cmd  .= $log_folder." ";
	$cmd  .= " >> ".$log_folder."/ML_module.log 2>> ".$log_folder."/ML_module.log";
	print LOG "[".`date | tr '\n' ' '`."] $cmd\n\n";
	$res = system($cmd);
	die "Error running classify.py . Exit code $res check ".$log_folder."/ML_module.log for further details." if $res!=0;
	die "Error Running classify.py . pegasus.output.txt file not generated." if ! -e $out_folder."/pegasus.output.txt";
	die "Error Running classify.py . pegasus.output.txt is empty." if -z $out_folder."/pegasus.output.txt";
	checkpoint($log_folder."/ML_module.job");
}
else
{
	print STDERR "[".`date | tr '\n' ' '`."] Running ML module [SKIPPED]\n";
}

printf STDERR "[".`date | tr '\n' ' '`."] Execution Terminated.\n";


close(LOG);

















sub printHelp
{
	print "Pegasus Usage:\n";
	print "\t-c config_file: path to Pegasus configuration file (mandatory)\n"; 
	print "\t-d data_input_file: path to data input configuration file (mandatory)\n"; 
	print "\t-l log_folder: path to a log folder (mandatory)\n"; 
	print "\t-o output_folder: path to a output folder (mandatory)\n"; 
	print "\t-p parallelization_sge: 1 run on a SGE system, 0 otherwise (optional, default 0)\n"; 
	print "\t-k keepDB: 1 the database is cleaned, 0 otherwise (optional, default 1)\n"; 
	print "\t-m model_classification: gbc for gradient boosting, rfc for random forest (optional, default gbc)\n"; 
}

sub checkpoint
{
	my $path = shift;
	open(FW,">".$path);
	print FW "";
	close(FW)
}

