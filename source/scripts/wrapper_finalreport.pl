#!/usr/bin/env perl

use strict;
use warnings FATAL => 'all';
use Getopt::Long;
use Scalar::Util qw(looks_like_number);
use File::Basename;

my $script_path;
my $working_dir;
my $input_file;
my $output_file;
my $logfile;
my $resources;

GetOptions (
	's=s' => \$script_path,  	# 
	'd=s' => \$working_dir,  	# 
	'i=s' => \$input_file,  	# 
	'l=s' => \$logfile,  		# 
	'r=s' => \$resources,  		# 
	'o=s' => \$output_file  	# 
);  

# 1 is true, 0 false
my $p_fail;

if(!defined($script_path))
{
        print STDERR "Error - $0: script_path not specified\n";
		$p_fail=1;
}

if(!defined($resources))
{
        print STDERR "Error - $0: resources not specified\n";
		$p_fail=1;
}

if(!defined($working_dir))
{
        print STDERR "Error - $0: working_dir not specified\n";
		$p_fail=1;
}

if(!defined($input_file))
{
        print STDERR "Error - $0: input_file not specified\n";
		$p_fail=1;
}

if(!defined($output_file))
{
        print STDERR "Error - $0: output_file file not specified\n";
		$p_fail=1;
}

if(!defined($logfile))
{
        print STDERR "Error - $0: logfile file not specified\n";
		$p_fail=1;
}

if ($p_fail)
{
	printHelp();
	exit(255);
}

open LOG, ">".$logfile or die "INFO - $0: Error opening $logfile file. \n";

my $cmd = "";

$cmd = "cat $input_file | cut -f36,38 | sed -e 's/\"//g' | tr '\t' '\@' | sed -e 's/\@//g' | sed -e 's/, /,/g' > $working_dir/Conserved_12.txt";
print LOG "[".`date | tr '\n' ' '`."] $cmd\n\n";
system($cmd);

$cmd = "cat $input_file | cut -f37,39 | sed -e 's/\"//g' | tr '\t' '\@' | sed -e 's/\@//g' | sed -e 's/, /,/g' > $working_dir/Lost_12.txt";
print LOG "[".`date | tr '\n' ' '`."] $cmd\n\n";
system($cmd);

$cmd = "perl $script_path/DomainsFormat.pl -r $resources/biological/Domains.Normal.ChimeraDB.txt -s $working_dir/Conserved_12.txt -l CD > $working_dir/results.CD.tsv";
print LOG "[".`date | tr '\n' ' '`."] $cmd\n\n";
system($cmd);

$cmd = "perl $script_path/DomainsFormat.pl -r $resources/biological/Domains.Normal.ChimeraDB.txt -s $working_dir/Lost_12.txt -l LD > $working_dir/results.LD.tsv";
print LOG "[".`date | tr '\n' ' '`."] $cmd\n\n";
system($cmd);

$cmd = "cat $input_file | cut -f14,15 > $working_dir/fusionsList.txt";
print LOG "[".`date | tr '\n' ' '`."] $cmd\n\n";
system($cmd);

$cmd = "perl $script_path/addCosmicToFusionList.pl -r $working_dir/fusionsList.txt -s $resources/biological/CosmicTable3.txt | cut -f1-4 > $working_dir/cosmic.txt";
print LOG "[".`date | tr '\n' ' '`."] $cmd\n\n";
system($cmd);

$cmd = "perl $script_path/CheckOverlappingDomains.pl -r $working_dir/Conserved_12.txt -s $resources/biological/Oncogene/ONCO_ALL.txt -l ONCO_CONSERVED > $working_dir/Onco_domains_CD.txt";
print LOG "[".`date | tr '\n' ' '`."] $cmd\n\n";
system($cmd);

$cmd = "perl $script_path/CheckOverlappingDDI.pl -r $working_dir/Conserved_12.txt -s $resources/biological/Oncogene/ONCO_DDI.txt -l ONCO_DDI_CONSERVED > $working_dir/Onco_DDI_CD.txt";
print LOG "[".`date | tr '\n' ' '`."] $cmd\n\n";
system($cmd);

$cmd = "perl $script_path/CheckOverlappingDomains.pl -r $working_dir/Lost_12.txt -s $resources/biological/Oncogene/ONCO_ALL.txt -l ONCO_LOST > $working_dir/Onco_domains_LD.txt";
print LOG "[".`date | tr '\n' ' '`."] $cmd\n\n";
system($cmd);

$cmd = "perl $script_path/CheckOverlappingDDI.pl -r $working_dir/Lost_12.txt -s $resources/biological/Oncogene/ONCO_DDI.txt -l ONCO_DDI_CONSERVED > $working_dir/Onco_DDI_LD.txt";
print LOG "[".`date | tr '\n' ' '`."] $cmd\n";
system($cmd);

$cmd = "perl $script_path/CheckOverlappingDomains.pl -r $working_dir/Conserved_12.txt -s $resources/biological/TumorSuppressors/TS_ALL.txt -l TS_CONSERVED > $working_dir/TS_domains_CD.txt";
print LOG "[".`date | tr '\n' ' '`."] $cmd\n\n";
system($cmd);

$cmd = "perl $script_path/CheckOverlappingDomains.pl -r $working_dir/Lost_12.txt -s $resources/biological/TumorSuppressors/TS_ALL.txt -l TS_LOST > $working_dir/TS_domains_LD.txt";
print LOG "[".`date | tr '\n' ' '`."] $cmd\n\n";
system($cmd);

$cmd = "perl $script_path/CheckOverlappingDDI.pl -r $working_dir/Conserved_12.txt -s $resources/biological/TumorSuppressors/TS_DDI.txt -l TS_DDI_CONSERVED > $working_dir/TS_DDI_CD.txt";
print LOG "[".`date | tr '\n' ' '`."] $cmd\n\n";
system($cmd);

$cmd = "perl $script_path/CheckOverlappingDDI.pl -r $working_dir/Lost_12.txt -s $resources/biological/TumorSuppressors/TS_DDI.txt -l TS_DDI_LOST > $working_dir/TS_DDI_LD.txt";
print LOG "[".`date | tr '\n' ' '`."] $cmd\n\n";
system($cmd);

$cmd = "paste $input_file $working_dir/Onco_domains_CD.txt $working_dir/Onco_DDI_CD.txt $working_dir/Onco_domains_LD.txt $working_dir/Onco_DDI_LD.txt $working_dir/TS_domains_CD.txt $working_dir/TS_domains_LD.txt $working_dir/results.CD.tsv $working_dir/results.LD.tsv | cut -f1-623,625- | cut -f1-1183 > $output_file";
print LOG "[".`date | tr '\n' ' '`."] $cmd\n\n";
system($cmd);

close(LOG);


sub printHelp
{
	print "\t-i [input file]\n";
	print "\t-o [output file]\n";
	print "\t-s [script_path]\n";
	print "\t-d [working_dir]\n"; 
	print "\t-l [logfile]\n"; 
}


